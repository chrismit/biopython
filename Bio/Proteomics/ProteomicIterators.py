"""
Author: Chris Mitchell (chris.mit7@gmail.com)
Copyright (C) 2012 Chris Mitchell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import re, os, ProteomicMasses, sqlite3, StringIO, zipfile, time, ProteomicObjects as PObj
import xml.etree.cElementTree as etree

class AnyIterator(object):
    """
    Tries to figure out what iterator we want and use it
    """
    def __init__(self, filename, **kwrds):
        ftype = os.path.splitext(filename)[1]
        iterTypes = ({'.xml':XTandemXMLIterator, '.msf':ThermoMSFIterator, '.mgf':MGFIterator})
        try:
            if kwrds:
                self.ipointer = iterTypes[ftype](filename,kwrds)
            else:
                self.ipointer = iterTypes[ftype](filename)
        except KeyError:
            return None
            
    def __iter__(self):
        return self.ipointer.__iter__()
    
    def next(self):
        return self.ipointer.next()
    
class XTandemXMLIterator(object):
    """
    Parser for X!Tandem XML Files.
    """
    def __init__(self, filename):
        #parse in our X!Tandem xml file
        dom1 = etree.parse(filename)
        self.scanSplit = re.compile(r'[\s\t]')
        self.group = dom1.findall("group")
        self.groupMap = {}
        self.db = None
        
    def __iter__(self):
        return self
    
    def parselxml(self, group):
        try:
            expect = group.attrib["expect"]
        except KeyError:
            if self.group:
                self.group.pop(0)
            self.next()
        subnote = list(group.iter("note"))
        for i in subnote:
            if (i.attrib["label"] == "Description"):
                experiment = i.text.strip()
        charge = group.attrib["z"]
        premass = group.attrib["mh"]
        rt = group.attrib["rt"]
        proteins = list(group.iter("protein"))
        fullProtein = ""
        for protein in proteins:
            scanObj = PObj.PeptideObject()
            scanObj.charge = charge
            scanObj.mass = premass*int(charge)
            scanObj.rt = rt
            sgroup = group.iter("group")
            for i in sgroup:
                #This is horridly inefficient...
                if "fragment ion mass spectrum" in i.attrib["label"]:
                    ab = i.iter('{http://www.bioml.com/gaml/}Xdata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
                            mz = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    ab = i.iter('{http://www.bioml.com/gaml/}Ydata')
                    for j in ab:
                        mzIter = j.iter('{http://www.bioml.com/gaml/}values')
                        for k in mzIter:
                            inten = [mval for mval in k.text.strip().replace('\n',' ').split(' ')]
                    for j,k in zip(mz,inten):
                        scanObj.addScan(j,k)
            domain = list(protein.iter("domain"))[0]#we only have one domain per protein instance
            note = list(protein.iter("note"))[0]#same here
            mods = list(protein.iter("aa"))#we can have multiple modifications
            if not self.db:
                files = list(protein.iter("file"))
                self.db = files[0].attrib["URL"]
            id = domain.attrib["id"]
            start = domain.attrib["start"]
            end = domain.attrib["end"]
            peptide = domain.attrib["seq"]
            pExpect = domain.attrib["expect"]
            for mod in mods:
                scanObj.addModification(mod.attrib["type"],int(mod.attrib["at"])-1,float(mod.attrib["modified"]), False)
            scanObj.peptide = peptide
            scanObj.expect = pExpect
            scanObj.id = id
            self.groupMap[id] = scanObj
            scanObj.title = id
            scanObj.accession = note.text
        return scanObj
    
    def next(self):
        try:
            group = self.group.pop(0)
        except IndexError:
            raise StopIteration
        return self.parselxml(group)
                
    def getScan(self, id):
        try:
            return self.groupMap[id]
        except KeyError:
            raise Exception(ValueError,"Id not found in scan index") 
    
class MGFIterator(object):
    def __init__(self, filename):
        #load our index      
        tFile = list(filename)
        tFile.reverse()
        tFile = tFile[tFile.index('.')+1:]
        tFile.reverse()
        self.indexFile=''.join(tFile)+'.mgfi'
        self.scanSplit = re.compile(r'[\s\t]')
        self.titleMap = {}
        self.distillerParse = re.compile(r'_DISTILLER_RAWFILE\[(\d+)\]=\(1\)(.+)')
        self.ra = {}
        if isinstance(filename,(str,unicode)):
            self.f = open(filename, 'rb')
        elif isinstance(filename,file):
            self.f = filename
        else:
            raise Exception(TypeError,"Unknown Type of filename -- must be a file handle or a file path")
        self.tparse = re.compile(r'TITLE=(\d+),(\d+): Scan (\d+) \(rt=(.+)\)')
        self.openMGFIndex()
        
    def openMGFIndex(self):
        path = self.indexFile
        try:
            f = open(path, 'rb')
            for row in f:
                entry = row.strip().split('\t')
                self.ra[entry[0]] = (int(entry[1]),int(entry[2]))
            try:
                self.epos = entry[2]
            except UnboundLocalError:
                self.epos=1
        except IOError:
            print 'no index available, call saveIndex() on iterator object to create one'
            
    def saveIndex(self):
        path = self.indexFile
        if os.path.exists(path):
            raise Exception("Index file path: %s found, but appears incomplete"%path)
        f = open(path, 'wb')
        while True:
            try:
                self.next()
            except StopIteration:
                break
        for i in self.ra:
            f.write('%s\t%d\t%d\n'%(i,self.ra[i][0],self.ra[i][1]))
        f.flush()
        self.f.seek(0)
            
    def getScan(self, title):
        """
        allows random lookup
        """
        if self.ra.has_key(title):
            self.f.seek(self.ra[title][0],0)
            toRead = self.ra[title][1]-self.ra[title][0]
            info = self.f.read(toRead)
            scan = self.parseScan(info)
        else:
            return None
        return scan
        
    def parseScan(self, scan):
        """
        All input follows the BEGIN IONS row and ends before END IONS
        """
        setupScan = True
        foundCharge = False
        foundMass = False
        foundTitle = False
        scanObj = PObj.ScanObject()
        for row in scan.split('\n'):
            if not row:
                continue
            entry = row.strip().split('=')
            if len(entry) >= 2:
                if entry[0] == 'PEPMASS':
                    scanObj.mass = entry[1]
                    foundMass = True
                elif entry[0] == 'CHARGE':
                    scanObj.charge = entry[1]
                    foundCharge = True
                elif entry[0] == 'TITLE':
                    title = '='.join(entry[1:])
                    foundTitle = True
                    scanObj.title = title
                elif entry[0] == 'RTINSECONDS':
                    scanObj.rt = entry[1]
                else:
                    scanObj.entry[0] = entry[1]
            else:
                mz,intensity = self.scanSplit.split(row.strip())
                scanObj.addScan(mz,intensity)
        if foundCharge and foundMass and foundTitle:
            return scanObj
        return None
            
    def __iter__(self):
        return self
    
    def next(self):
        row = self.f.readline()
        if row == '':
            raise StopIteration
        setupScan = False
        scanInfo = ""
        while row:
            if '_DISTILLER' in row:
                if row:
                    m = self.distillerParse.match(row)
                    if m:
                        self.titleMap[int(m.group(1))+1] = m.group(2)
            elif 'BEGIN IONS' in row:
                if self.rand:
                    pStart=self.f.tell()
                setupScan=True
            elif 'END IONS' in row:
                scan = self.parseScan(scanInfo)
                if scan:
                    self.ra[scan.getTitle()] = (pStart,pos)
                    return scan
                return None
            elif setupScan:
                scanInfo+=row
            pos = self.f.tell()
            row = self.f.readline()
            
class ThermoMSFIterator(object):
    """
    Thermo Scientific's Proteome Discoverer Parser
    Optional Keywords:
    full=True(False default) to parse all information in, much slower
    confidence=1(1 default) minimum confidence level
    rank=1(1 default) minimun search engine rank confidence
    """
    def __init__(self, filename, full=False, confidence=1, rank=1):
        lastSplit = re.compile(r'.+[/\\](.+)')
        if isinstance(filename,(str,unicode)):
            self.f = open(filename, 'rb')
        else:
            raise Exception(TypeError,"Unknown Type of filename -- must be a file path")
        self.conn = sqlite3.connect(filename, check_same_thread=False)
        self.full = full
        self.cur = self.conn.cursor()
        sql = 'select * from fileinfos'
        self.cur.execute(sql)
        self.fileMap = {}
        self.sFileMap = {}
        self.scans = []
        for i in self.cur.fetchall():
            self.fileMap[str(i[0])]=str(i[1])
            self.sFileMap[i[0]]=lastSplit.search(i[1]).group(1)
        #modification table
        sql = 'select a.AminoAcidModificationID,a.ModificationName, a.DeltaMass from aminoacidmodifications a'
        self.cur.execute(sql)
        self.modTable = {}
        for i in self.cur.fetchall():
            self.modTable[i[0]] = (i[1],i[2])
        #We fetch all modifications here for temporary storage because it is VERY expensive to query peptide by peptide (3 seconds per 100 on my 500 MB test file, with 300,000 scans that's horrid)
        sql = 'select pam.PeptideID, GROUP_CONCAT(pam.AminoAcidModificationID), GROUP_CONCAT(pam.Position) from peptidesaminoacidmodifications pam GROUP BY pam.PeptideID'
        self.cur.execute(sql)
        self.mods = {}
        for i in self.cur.fetchall():
            self.mods[i[0]] = i[1:]
        try:
            if self.full:
                sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.SearchEngineRank),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID, sp.Spectrum from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join spectra sp on (sp.UniqueSpectrumID=sh.UniqueSpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel <= ? and p.SearchEngineRank <= ? GROUP BY p.SpectrumID'
            else:
                sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.SearchEngineRank),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel <= ? and p.SearchEngineRank <= ? GROUP BY p.SpectrumID'                
            self.cur.execute(sql,(confidence,rank))
        except sqlite3.OperationalError:#older version of PD created this file
            if self.full:
                sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID, sp.Spectrum from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join spectra sp on (sp.UniqueSpectrumID=sh.UniqueSpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel = ? GROUP BY p.SpectrumID'
            else:
                sql = 'select GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.ConfidenceLevel),GROUP_CONCAT(p.Sequence),GROUP_CONCAT(p.PeptideID), GROUP_CONCAT(pp.ProteinID), p.SpectrumID, sh.Charge, sh.RetentionTime, sh.FirstScan, sh.LastScan, mp.FileID from peptides p join peptidesproteins pp on (p.PeptideID=pp.PeptideID) left join spectrumheaders sh on (sh.SpectrumID=p.SpectrumID) left join masspeaks mp on (sh.MassPeakID=mp.MassPeakID) where p.PeptideID IS NOT NULL and p.ConfidenceLevel = ? GROUP BY p.SpectrumID'
            self.cur.execute(sql,(confidence,))
        self.index = 0
            
    def getScan(self, specId, peptide):
        """
        get a random scan
        """
        sql = "select sp.Spectrum, p.Sequence, p.PeptideID from spectrumheaders sh left join spectra sp on (sp.UniqueSpectrumID=sh.UniqueSpectrumID) left join peptides p on (sh.SpectrumID=p.SpectrumID) where sh.SpectrumID = %d and p.Sequence = '%s'"%(int(specId),peptide)
        self.cur.execute(sql)
        i = self.cur.fetchone()
        if not i:
            return None
        return self.parseFullScan(i)
            
    def __iter__(self):
        return self
    
    def processZip(self, specZip, scanObj, pid):
        sInfo = specZip
        fp = StringIO.StringIO(sInfo)
        zf = zipfile.ZipFile(fp, 'r')
        sql = 'select aam.ModificationName,pam.Position,aam.DeltaMass from peptidesaminoacidmodifications pam left join aminoacidmodifications aam on (aam.AminoAcidModificationID=pam.AminoAcidModificationID) where pam.PeptideID=%s'%pid
        for row in self.conn.execute(sql):
            scanObj.addModification(scanObj.peptide[row[1]], str(row[1]), str(row[2]), row[0])
        for j in zf.namelist():
            msInfo = zf.read(j)
            msStr = msInfo.split('\n') 
            stage = 0
            #this is dirty, but unfortunately the fastest method at the moment -- lxml/celementree is far too slow
            for row in msStr:
                if stage == 0:
                    if 'FileID' in row:
                        stage=1
                elif stage == 1:
                    if 'PrecursorInfo' in row:
                        finfo = row.split('"')
                        smass = finfo[5]
                        scanObj.mass = smass
                        stage=2
                elif stage == 2:
                    if '<PeakCentroids>' in row:
                        stage = 3
                elif stage == 3:
                    #grab the ms/ms peaks
                    if 'Peak X' in row:
                        finfo = row.split('"')
                        scanObj.addScan(finfo[1],finfo[3])
                    elif '</PeakCentroids>' in row:
                        break
        if msInfo:
            return scanObj
        return None
    
    def parseScan(self, i):
        objs = []
        self.index+=1
        added = set([])#for some reason redundant scans appear
        for confidence, searchRank, sequence, pepId, proId in zip(i[0].split(','),i[1].split(','),i[2].split(','),i[3].split(','),i[4].split(',')):
            if (sequence,pepId,proId) in added:
                continue
            else:
                added.add((sequence,pepId,proId))
            scanObj = PObj.PeptideObject()
            try:
                mods = self.mods[int(pepId)]
                for modId, modPosition in zip(mods[0].split(','),mods[1].split(',')):
                    modEntry = self.modTable[int(modId)]
                    scanObj.addModification(sequence[int(modPosition)], modPosition, modEntry[1], modEntry[0])
            except KeyError:
                pass
            scanObj.peptide = sequence
            scanObj.rank = searchRank
            scanObj.confidence = confidence
            scanObj.accession = proId
            scanObj.charge = i[6]
            fName = self.sFileMap[i[10]]
            fScan = i[8]
            lScan = i[9]
            sid = '%s.%s.%s'%(fName, fScan,lScan)
            scanObj.title = sid
            scanObj.id = i[5]
            scanObj.spectrumId=i[5]
            if self.full:
                scanObj = self.processZip(i[11],scanObj, pepId)
            objs.append(scanObj)
        return objs
    
    def parseFullScan(self, i):
        """
        parses scan info, takes significantly longer since it has to unzip/parse xml
        """
        scanObj = PObj.PeptideObject()
        sInfo = i[0]
        fp = StringIO.StringIO(sInfo)
        zf = zipfile.ZipFile(fp, 'r')
        peptide = str(i[1])
        pid=i[2]
        sql = 'select aam.ModificationName,pam.Position,aam.DeltaMass from peptidesaminoacidmodifications pam left join aminoacidmodifications aam on (aam.AminoAcidModificationID=pam.AminoAcidModificationID) where pam.PeptideID=%s'%pid
        for row in self.conn.execute(sql):
            scanObj.addModification(peptide[row[1]], str(row[1]), str(row[2]), row[0])
        scanObj.peptide = peptide
        for j in zf.namelist():
            msInfo = zf.read(j)
            msStr = msInfo.split('\n') 
            stage = 0
            #this is dirty, but unfortunately the fastest method at the moment
            for row in msStr:
                if stage == 0:
                    if 'FileID' in row:
                        finfo = row.split('"')
                        fileName = int(finfo[1])
                        #msScanSum = finfo[5]
                        fName = self.sFileMap[fileName]
                        sid = '%s.%s.%s'%(fName, finfo[2],finfo[2])
                        scanObj.title = sid
                        scanObj.id = sid
                        stage=1
                elif stage == 1:
                    if 'PrecursorInfo' in row:
                        finfo = row.split('"')
                        charge = finfo[3]
                        smass = finfo[5]
                        scanObj.charge = charge
                        scanObj.mass = smass
                        stage=2
                elif stage == 2:
                    if '<PeakCentroids>' in row:
                        stage = 3
                elif stage == 3:
                    #we just grab the ms/ms peaks at the moment
                    if 'Peak X' in row:
                        finfo = row.split('"')
                        scanObj.addScan(finfo[1],finfo[3])
                    elif '</PeakCentroids>' in row:
                        break
        if msInfo:
            return scanObj
        else:
            return None
    
    def next(self):
        if not self.scans:
            i = self.cur.fetchone()
            #we go by groups
            if not i:
                raise StopIteration
            self.scans = self.parseScan(i)
            if not self.scans:
                raise StopIteration
        return self.scans.pop(0)