#!/usr/bin/env python
from I3Tray import *

from icecube import icetray, dataclasses, dataio

from hepevt_particle import HepEvtParticle

from math import sqrt

from os.path import expandvars
import sys
import gzip

    
class HEPEvtFileReaderModule(icetray.I3Module):
    ''' An I3Module to read in a hepevt file and convert to I3 format. It was originally written
    to read in events from GENIE, put is not specific to this generator. Intermediate states are
    ignored.

    This module has 4 parameters:
        - Filename            Name of input file .hepevt.gz
        - AddHadronic         Add all particles from hadronic origin together as one system
                              Default: True
        - AttachAllToLepton   Attach all particles as child to 1st primary lepton
                              Default: True
        - MomentumThreshold   Momentum threshold for particles to be added
                              Default: 1.*I3Units.keV

      M. Dierckxsens 2010/01/20
    '''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")
        self.AddParameter('Filename', \
                          'Name of file to read',\
                          'output.hepevt.gz')
        self.AddParameter('AddHadronic', \
                          'Add all the particles coming from the nucleus',\
                          True)
        self.AddParameter('AttachAllToLepton', \
                          'Attach all the particles as a child to the first primary lepton/neutrino',\
                          True)
        self.AddParameter('MomentumThreshold', \
                           'Threshold for momentum for particles to be added',\
                           1.*I3Units.keV)

    def Configure(self):
        self.filename = self.GetParameter('Filename')
        self.addhadronic = self.GetParameter('AddHadronic')
        self.attachtolepton = self.GetParameter('AttachAllToLepton')
        self.pthres = self.GetParameter('MomentumThreshold')

        ###
        # Open the file and keep a pointer as a "member variable"
        ###        
        self.f = gzip.open(self.filename,'rb')

    def Physics(self, frame):
        
        mctree = dataclasses.I3MCTree()

        ###
        # Read a line from the file
        ###
        line = self.f.readline()

        ###
        # If an empty line is returned this means the end of the file is reached
        # In that case I3Module::RequestSuspension is called to stop the IceTray flow.
        ###
        if len(line) == 0 :
            self.RequestSuspension()
            return

        # ignore possible lines starting with #
        while line[0] == '#' :
            line = self.f.readline()


        ###
        # First line of an event is event number and number of particles
        ###
        eventNumber = int(line.split()[0])
        numberParticles = int(line.split()[1])
        #print eventNumber, numberParticles
        ###
        #Hadronic sum, if required
        ###
        hadsum = HepEvtParticle()
        hadsum.idhep=99
        nhadorigin=0
        firsthadind=-1
        firstleptonind=-1
        ###
        # The list that will hold all the HepEvtParticles 
        ###
        particlelist=[]
        ###
        # The list that will hold all the primary particles to be saved
        ###
        primI3PartList=[]
        ###
        # The list of indices to the primary particles
        ###
        ind_prims_added=[]
        ###
        # Read in numberParticles lines, get HEPEVtPartile and add to the list
        # If primary particle, add them to the mctree, depending on requested
        # readout scheme
        ###
        for ip in range(numberParticles):
            firstlepton=False
            firsthadron=False
            islepton=False
            line = self.f.readline()
            hep = HepEvtParticle()
            hep.ReadLine(line)
            #hep.Print()
            particlelist.append(hep)
            if hep.isthep==0:
                # save the index of the first lepton
                if (abs(hep.idhep)>10 and abs(hep.idhep)<17):
                    islepton=True
                    if firstleptonind<0:
                        firstleptonind=ip
                        firstlepton=True
                # save the index of the first non-leptonic particle
                else :
                    if firsthadind<0:
                        firsthadind=ip
                        firsthadron=True
                ###
                # Add the primaries and save them in a list 
                # This depends on the options for how to attach the final state
                # particles to the primaries
                ###
                addparticle=False
                #
                # if everything has to be attached to the lepton, only add this
                # to primary list if it is the first found lepton
                if self.attachtolepton:
                    if firstlepton:
                        addparticle=True
                #
                # if hadronic particles are to be added, save all leptons
                # + first hadronic particle
                elif self.addhadronic:
                    if islepton or firsthadron:
                        addparticle=True
                #
                # in other cases just add all the particles
                else :
                    addparticle=True
                    
                if addparticle:
                    ind_prims_added.append(ip)
                    i3part = hep.GetI3Particle()
                    primI3PartList.append(i3part)
                    mctree.AddPrimary(i3part)                
        ###
        # Now find all the final state particles
        ###
        for ind in range(len(particlelist)):
            heppart = particlelist[ind]
            # only consider final state particles
            if heppart.isthep != 1 or heppart.jmohep1 == -1: continue
            # ignore slow final state particles
            ptot = heppart.phepx*heppart.phepx
            ptot += (heppart.phepy*heppart.phepy)
            ptot += (heppart.phepz*heppart.phepz)
            ptot = sqrt(ptot)
            if ptot < self.pthres: continue
            # find primary mother particle
            misthep = 1
            mind = ind
            nstep = 0
            while misthep != 0 :
                mind = particlelist[mind].jmohep1
                try: misthep = particlelist[mind].isthep
                except IndexError :
                    print "Warning: Index out of particlelist range"
                    break
                nstep += 1
                # just a check when ending up in an infinite loop 
                if nstep>20:
                    print "Warning: more than 20 steps in finding primary mother particle"
                    break
            #print "prim mother index ",mind, " with pdg id ", particlelist[mind].idhep
            if misthep != 0 :
                print "Did not find mother of final state particle, will not be added to the tree!"
                continue

            # if the particle doesn't come from a lepton, it is said to be from hadronic
            # origin and will be added to the hadronic particle (shower) if requested
            hadronicOrigin=False
            if (abs(particlelist[mind].idhep)<10 or \
                abs(particlelist[mind].idhep)>17):
                hadronicOrigin=True
                nhadorigin += 1
            if self.addhadronic and hadronicOrigin:
                # add particle momenta & energy to hadsum
                hadsum.phepx += heppart.phepx
                hadsum.phepy += heppart.phepy
                hadsum.phepz += heppart.phepz                
                hadsum.phepe += heppart.phepe
                # this will result in the vertex and time set from the last particle
                hadsum.vhepx = heppart.vhepx
                hadsum.vhepy = heppart.vhepy
                hadsum.vhepz = heppart.vhepz
                hadsum.vhept = heppart.vhept
                #skip adding the particle to the tree, it will be done after the loop
                continue

            # get the index in the list that holds the indices of the primary particles
            # It will depend on the requested manner of attaching the final state
            # particles to their parent
            #
            # if all particle have to be attached to the lepton, the index will be of
            # the first lepton in the primary particle list
            # the case where hadronic particle are summed does not reach this step
            # the one remaining case just keeps the previous found value for mind
            if self.attachtolepton :
                if firstleptonind<0:
                    print "Requested to add all particles to the primary lepton " \
                          "but did not found one! Particle will be attached to " \
                          "the original primary"
                else:
                    mind = firstleptonind
            # find the index in the primary i3particle list of the mother particle 
            prind=-1            
            try: prind=ind_prims_added.index(mind) 
            except ValueError :
                print "Primary with index ",mind," not found in primary particle list" 
                continue

            fsI3Part = heppart.GetI3Particle()
            mctree.AppendChild(primI3PartList[prind],fsI3Part)

        # Add sum of hadronic particle if requested, and if particles where added
        if self.addhadronic and nhadorigin>0 :
            hadsum.phepm = hadsum.phepe**2
            hadsum.phepm -= hadsum.phepx**2
            hadsum.phepm -= hadsum.phepy**2
            hadsum.phepm -= hadsum.phepz**2
            if hadsum.phepm<0: hadsum.phepm=0

            hadsum.phepm = sqrt(hadsum.phepm)

            hadsum.phepe = sqrt( hadsum.phepx**2 + hadsum.phepy**2 + hadsum.phepz**2)

            fsI3parthadsum = hadsum.GetI3Particle()

            mind=-1
            if self.attachtolepton :
                if firstleptonind<0:
                    print "Requested to add summed hadronic particles to the primary "\
                          "lepton but did not found one! Particle will be attached " \
                          "to the first hadronic primary"
                else:
                    mind = firstleptonind
            if mind<0 :
                if firsthadind<0:
                    print "Did not find first hadronic primary!"
                else :
                    mind = firsthadind
            

            # find the index in the primary i3particle list of the mother particle 
            prind=-1            
            try: prind=ind_prims_added.index(mind) 
            except ValueError :
                if mind < 0:
                    print " Did not find primary particle "
                else:
                    print "Primary with index ",mind," not found in primary " \
                          "particle list" 
                print "Summed hadronic system will not be added to the MC tree"

            if prind>=0 :
                mctree.AppendChild(primI3PartList[prind],fsI3parthadsum)
                
        #sim_services.I3GeoShifter.ShiftTree(frame,mctree)
        frame['I3MCTree'] = mctree

        #add event id to the eventheader
        evthdr = frame.Get('I3EventHeader')
        evthdr.EventID=eventNumber
        
        self.PushFrame(frame)


    
            
