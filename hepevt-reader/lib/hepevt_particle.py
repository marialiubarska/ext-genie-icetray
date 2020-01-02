from I3Tray import *

from icecube import dataclasses

from math import acos, sqrt, atan2

class HepEvtParticle:

    # the members
    isthep=-1
    idhep=-1
    jmohep1=-1
    jmohep2=-1
    jdahep1=-1
    jdahep2=-1
    phepx=0.0
    phepy=0.0
    phepz=0.0
    phepe=0.0
    phepm=0.0
    vhepx=0.0
    vhepy=0.0
    vhepz=0.0
    vhept=0.0

    def ReadLine(self,line) : 
        self.isthep = int(line.split()[0])
        self.idhep = int(line.split()[1])
        self.jmohep1 = int(line.split()[2])
        self.jmohep2 = int(line.split()[3])
        self.jdahep1 = int(line.split()[4])
        self.jdahep2 = int(line.split()[5])
        self.phepx = float(line.split()[6]) * I3Units.GeV
        self.phepy = float(line.split()[7]) * I3Units.GeV
        self.phepz = float(line.split()[8]) * I3Units.GeV
        self.phepe = float(line.split()[9]) * I3Units.GeV
        self.phepm = float(line.split()[10]) * I3Units.GeV
        # read in vertex position and time
        # position is in mm, time in mm/c 
        self.vhepx = float(line.split()[11]) * I3Units.mm
        self.vhepy = float(line.split()[12]) * I3Units.mm
        self.vhepz = float(line.split()[13]) * I3Units.mm
        self.vhept = float(line.split()[14]) * I3Units.mm / (2.998e8 * I3Units.m/I3Units.s)


    def Print(self) : 
        print "status code = ", self.isthep
        print "PDG ID = ", self.idhep
        print "Mother index = ",self.jmohep1, " ",self.jmohep2
        print "Daughter index = ",self.jdahep1, " ",self.jdahep2
        print "(px,py,pz) [GeV/c] = (",self.phepx,", ",self.phepy,", ",self.phepz,")"
        print "(E,m) [GeV,GeV/c**2] = (",self.phepe,", ",self.phepm,")"
        print "(x,y,z,t) [m,m,m,ns] =(",self.vhepx,", ",self.vhepy,", ",\
              self.vhepz,", ",self.vhept,")"
                

    def GetI3Particle(self) :
        ''' Translate HepEvt particle to I3 particle
        '''
        i3part = dataclasses.I3Particle()
        if self.idhep == -11 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.EPlus)
        elif self.idhep == 11 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.EMinus)
        elif self.idhep == -12 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuEBar)
        elif self.idhep == 12 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuE)
        elif self.idhep == -13 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.MuPlus)
        elif self.idhep == 13 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.MuMinus)
        elif self.idhep == -14 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuMuBar)
        elif self.idhep == 14 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuMu)
        elif self.idhep == -15 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.TauPlus)
        elif self.idhep == 15 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.TauMinus)
        elif self.idhep == -16 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuTauBar)
        elif self.idhep == 16 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.NuTau)
        elif self.idhep == 22 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.Gamma)
        elif self.idhep == 111 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.Pi0)
        elif self.idhep == 211 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.PiPlus)
        elif self.idhep == -211 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.PiMinus)
            # not really correct but "closest" match
        elif self.idhep == 311 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.K0_Long)
            # not really correct but "closest" match
        elif self.idhep == -311 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.K0_Short)
        elif self.idhep == 321 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.KPlus)
        elif self.idhep == -321 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.KMinus)
        elif self.idhep == 2112 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.Neutron)
        elif self.idhep == 2212 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.PPlus) 
        elif self.idhep == -2212 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.PMinus) 
        elif self.idhep == 1000080160 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.O16Nucleus)
        # used for the sum of the hadronic particles    
        elif self.idhep == 99 :
            i3part.SetType(dataclasses.I3Particle.ParticleType.Hadrons)
        else :
            print "Unknown particle with pdg ID: ", self.idhep
            i3part.SetType(dataclasses.I3Particle.ParticleType.unknown)
                
        i3part.SetEnergy(self.phepe)
                
        ptot = sqrt(self.phepx**2 + self.phepy**2 + self.phepz**2)
        theta = 0
        if abs(ptot)>0 : theta = acos(self.phepz/ptot) * I3Units.radian 
        phi = atan2(self.phepy,self.phepx) * I3Units.radian 
        
        i3part.SetThetaPhi(theta,phi)
        i3part.SetPos(self.vhepx,self.vhepy,self.vhepz)
        i3part.SetTime(self.vhept)

        i3part.SetLocationType(dataclasses.I3Particle.LocationType.InIce)
        return i3part

