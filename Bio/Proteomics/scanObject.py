class scanObject(object):
    """
    A scan object to store peaklist information in
    """
    def __init__(self):
        self.scans = []
        pass
    
    def addTitle(self, title):
        self.title = title
        
    def addCharge(self, charge):
        self.charge = charge
        
    def addMass(self, mass):
        self.mass = mass
        
    def addScan(self, scan):
        s = scan.split(' ')
        if len(s) > 1 and float(s[1]) != 0.0:
            self.scans.append(scan)
        
    def getInfo(self):
        return self.scans
    
    def getCharge(self):
        return self.charge
    
    def addRT(self, rt):
        self.rt = rt
    
    def getTitle(self):
        return self.title
    
    def getRT(self):
        try:
            return self.rt
        except:
            return None
        
    def getMZ(self):
        out = []
        for i in self.scans:
            out.append(i.split(' '))
        return out
    
    def getPrecursor(self):
        return self.mass
    
    def writeScan(self, o):
        o.write('BEGIN IONS\n')
        o.write('TITLE=%s\n'%self.title)
        try:
            o.write('RTINSECONDS=%s\n'%self.rt)
        except:
            pass
        o.write('PEPMASS=%s\n'%self.mass)
        o.write('CHARGE=%s\n'%self.charge)
        for i in self.scans:
            o.write('%s\n'%i)
        o.write('END IONS\n\n')