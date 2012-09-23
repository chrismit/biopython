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
import ProteomicMasses as PMass

class ScanObject(object):
    """
    A scan object to store peaklist information in
    attributes to populate: title, charge, mass, rt
    """
    def __init__(self):
        self.scans = []
        self.rt = ''#populated in init so if it's fetched, it won't throw an error
            
    def addScan(self, mz, inten):
        self.scans.append((mz,inten))
    
    def writeScan(self, o):
        o.write('BEGIN IONS\n')
        o.write('TITLE=%s\n'%self.title)
        if self.rt:
            o.write('RTINSECONDS=%s\n'%self.rt)
        o.write('PEPMASS=%s\n'%self.mass)
        o.write('CHARGE=%s\n'%self.charge)
        for i in self.scans:
            o.write('%s %s\n'%(i[0],i[1]))
        o.write('END IONS\n\n')
        
class PeptideObject(ScanObject):
    """
    An enhanced scan object that can store peptide information as well.
    attributes to populate: peptide, expect, id, accession 
    """
    def __init__(self):
        ScanObject.__init__(self)
        self.mods = set([])
        self.peptide = ""
        
    def addModification(self, aa,position, modMass, modType):
        """
        MODIFICATION POSITION IS 0 BASED
        Modifications are stored internally as a tuple with this format:
        (amino acid modified, index in peptide of amino acid, modification type, modification mass)
        ie (M, 7, Oxidation, 15.9...)
        """
        if not modType:
            #try to figure out what it is
            tmass = abs(modMass)
            smass = str(tmass)
            prec = len(str(tmass-int(tmass)))-2
            precFormat = '%'+'0.%df'%prec
            modType = ""
            for i in PMass.mod_weights:
                if tmass in PMass.mod_weights[i] or smass == precFormat%PMass.mod_weights[i][0]:
                    #found it
                    modType = i
            if not modType:
                print 'mod not found',modMass
        self.mods.add((aa,str(position),str(modMass),str(modType)))
        
    def getModifications(self):
        return '|'.join([','.join(i) for i in self.mods])