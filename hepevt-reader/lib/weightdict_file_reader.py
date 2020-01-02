#!/usr/bin/env python
from I3Tray import *

from icecube import icetray, dataclasses, dataio

#from hepevt_particle import HepEvtParticle

from math import sqrt

from os.path import expandvars
import sys
import gzip

    
class WDictFileReaderModule(icetray.I3Module):
    ''' An I3Module to read in a wdict file and convert to I3 format (I3WeightDict).

    This module has 1 parameter:
        - Filename            Name of input file .wdict.dat.gz

      M. Dierckxsens 2010/09/19
    '''
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")
        self.AddParameter('Filename', \
                          'Name of file to read',\
                          'output.hepevt.gz')

    def Configure(self):
        self.filename = self.GetParameter('Filename')
 
        ###
        # Open the file and keep a pointer as a "member variable"
        ###
        if self.filename.find(".gz") != '-1':
            self.f = gzip.open(self.filename,'rb')
        else:
            self.f = open(self.filename, 'r')


    def Physics(self, frame):
        
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
        # First line of an event is event number, number of particles
        ###
        #print "line is ",line
        eventNumber = int(line.split()[0])
        numberMapEntries = int(line.split()[1])
        numberParticles = int(line.split()[2])
        nuPdg = int(line.split()[3])
        nuE = float(line.split()[4])

        ## @@check if we are dealing with the same event!
        ## get the physics frame and compare with values 2,3 & 4
        #frame['I3MCTree'] = mctree

        #add event id to the eventheader
        #evthdr = frame.Get('I3EventHeader')
        #evthdr.EventID=eventNumber
        
        #print "evtnum etc is ",eventNumber, numberMapEntries, numberParticles,nuPdg, nuE

        weightdict = dataclasses.I3MapStringDouble()
        for ip in range(numberMapEntries):
            line = self.f.readline()
            wdkey = str(line.split()[0])
            wdvalue = float(line.split()[1])
            weightdict[wdkey]=wdvalue

        #for kk,vv in weightdict.iteritems():
        #    print "   >>> ", kk, vv

        frame['I3MCWeightDict']=weightdict
        
        self.PushFrame(frame)


