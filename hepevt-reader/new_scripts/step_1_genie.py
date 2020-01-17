#!/usr/bin/env python

# This file has very small modifications from the original example script. It adds MMC for NuMu files.

from optparse import OptionParser
from os.path import expandvars
import os, random

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="test_genie.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-l", "--filenr",type="int",default=1,
                   dest="FILENR", help="File number, stream of I3SPRNGRandomService")
parser.add_option("-r", "--runnumber", type="int", default=1,
                  dest="RUNNUMBER", help="The run number for this simulation, also the seed for random number simulations")
parser.add_option("-n", "--numevents", type="int", default=100,
                  dest="NUMEVENTS", help="The number of events per run")
parser.add_option("-f", "--flavor", default="NuMu", 
                  dest="FLAVOR", help="The flavor of the neutrino produced")
parser.add_option("", "--energy-range", default="P", 
                  dest="ENERGYRANGE", help = "A=lowest, B=medium, C=high, P=production")
# python step_1_genie.py -o ~/scratch/simulation/L2/test_alex.i3.gz -s 4332 -r 1 -n 1000 -f NuMu

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()
if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
                crap += a
                crap += " "
        parser.error(crap)

print 'Using RUNNUMBER: ', options.RUNNUMBER
if options.FLAVOR == "NuE":
    en_range = [0.0,0.0]
    spectralIndex = 2.5
    if options.ENERGYRANGE=="A":
        spectralIndex = 1.9
        en_range[0]  = 1.0; en_range[1] = 4.0
    elif options.ENERGYRANGE=="B":
        spectralIndex = 2.0
        en_range[0]  = 4.0; en_range[1] = 12.0
    elif options.ENERGYRANGE=="C":
        en_range[0]  = 12.0; en_range[1] = 9999.999999
    elif options.ENERGYRANGE=="Z":
        spectralIndex = 1.
        en_range[0]  = 0.5; en_range[1] = 1.0
    else:
        print "Unknown option, return!"
        exit()
    radius = 250.0 ; length = 500.0
elif options.FLAVOR == "NuMu":
    en_range = [0.0,0.0]
    spectralIndex = 2.
    if options.ENERGYRANGE=="A":
        en_range[0] = 1.0; en_range[1] = 5.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="B":
        en_range[0] = 5.0; en_range[1] = 80.0
        radius = 330.0 ; length = 900.0
    elif options.ENERGYRANGE=="C":
#        en_range[0] = 80.0; en_range[1] = 9999.999999
        en_range[0] = 9998.0; en_range[1] = 9999.999999
        radius = 330.0 ; length = 1500.0
    elif options.ENERGYRANGE=="P":
        en_range[0]  = 3.0; en_range[1] = 9999.999999
        radius = 330.0 ; length = 1200.0
        print "Setting P option!!! not divided on ranges"
    elif options.ENERGYRANGE=="Z":
        spectralIndex = 1.
        en_range[0] = 0.5; en_range[1] = 1.0
        radius = 250.0 ; length = 500.0
    else:
        print "Unknown option, return!"
        exit()
elif options.FLAVOR == "NuTau":
    en_range = [0.0,0.0]
    spectralIndex = 2.   
    if options.ENERGYRANGE=="A":
        spectralIndex = 1.5
        en_range[0]  = 4.0; en_range[1] = 10.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="B":
        spectralIndex = 2.00
        en_range[0]  = 10.0; en_range[1] = 30.0
        radius = 250.0 ; length = 500.0
    elif options.ENERGYRANGE=="C":
        spectralIndex = 3.5
        en_range[0]  = 30.0; en_range[1] = 9999.999999  
        radius = 250.0 ; length = 1000.0    
    else: 
        print "Unknown option! Return!"
        exit()
else:
    en_range = [0.0,0.0]
    spectralIndex = 2.
    en_range[0]  = 3.0; en_range[1] = 9999.999999
print "Radius: ", radius
print "Length: ", length
print "Energy range " + options.ENERGYRANGE+ ": [" + str(en_range[0]) + " , " + str(en_range[1]) + " ] GeV"
    

from I3Tray import *
import os
import sys

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import genie_icetray



load("libsim-services")


tray = I3Tray()

# Random number generator
randomService = phys_services.I3SPRNGRandomService(
        seed = options.RUNNUMBER, # +10 for files, which want to fail anyway (usually 1 file per a couple of thousands)
        nstreams = 50000,
        streamnum = options.FILENR)


tray.AddModule("I3InfiniteSource","streams",
	       Stream=icetray.I3Frame.DAQ)

tray.AddModule("I3MCEventHeaderGenerator","gen_header",
	       Year=2009,
	       DAQTime=158100000000000000,
	       RunNumber=1,
	       EventID=1,
	       IncrementEventID=True)

# Generate the neutrino interactions
tray.AddModule("I3GENIEGenerator","genie_generator",
    RandomService = randomService, 
#    SplineFilename = expandvars("/afs/ifh.de/user/t/terliuk/lustre/fs19/group/icecube/terliuk/GENIE_splines/GENIE_water_splines_10TeV_N500.xml"),
#    SplineFilename = expandvars("/cvmfs/icecube.opensciencegrid.org/data/genie-splines/GENIE_water_splines_10TeV_N500.xml"),
#    SplineFilename = expandvars("/home/hignight/project/hignight/genie_2_12_8_splines/GENIE_2_12_8_Water_splines.xml"),
    SplineFilename = expandvars("/home/mliubar/Software/genie_workspace/genie-generator/xsec_splines/GENIE_2_12_8_Water_splines.xml"),
    LHAPDFPath = expandvars("$I3_BUILD/genie-icetray/resources/PDFsets"),
    NuEnergyMin = en_range[0]*I3Units.GeV, #3, 195
    NuEnergyMax = en_range[1]*I3Units.GeV,
    PowerLawIndex = spectralIndex, # E^-2.5 spectrum
    GenVolRadius = radius*I3Units.m,
    GenVolLength = length*I3Units.m,
    GenVolDepth = 1950.*I3Units.m, # This option does nothing. Changed it manually in file.
    NeutrinoFlavor = options.FLAVOR, 
    NuFraction = 0.70, # to match lifetime
    MaterialDensity = 0.93*I3Units.g/I3Units.cm3, # ice density
    TargetMixIngredients = [1000080160,1000010010], # O16, H1
    TargetMixQuantities = [1,2], # H2O (O16->1x, H1->2x)
    ForceSingleProbScale = False,
    NEvents = options.NUMEVENTS,
    SystematicNames = ['MaCCRES','MaNCRES','MaCCQE','MaNCEL','MaCOHpi','AhtBY','BhtBY','CV1uBY','CV2uBY'],
    SystematicSteps = [-2,-1,1,2],
    OutputGST=True,
    PositionShift=dataclasses.I3Position(46.29,-34.88,-330.0),  # Changed from dataclasses.I3Position(40,-50,-40), 
    MCTreeName = "I3MCTree_GENIE"
)


# Set up the Driving Time
time = dataclasses.I3Time()
time.set_mod_julian_time(55697, 0, 0)
def DrivingTime( frame ):
	if "DrivingTime" in frame : 
		del frame["DrivingTime"]
	frame.Put("DrivingTime", time )
def NEvents( frame ):
    if "NEvPerFile" in frame:
        del frame['NEvPerFile']
    frame.Put('NEvPerFile', icetray.I3Int(options.NUMEVENTS))

tray.AddModule(DrivingTime, "dt",
	       Streams = [icetray.I3Frame.DAQ] )
tray.AddModule(NEvents, "ne",
	       Streams = [icetray.I3Frame.DAQ] )
if not options.FLAVOR == 'NuMu':
	# Pass the results into an I3MCTree
	tray.AddModule("I3GENIEResultDictToMCTree", "toMcTree",
                       MCTreeName="I3MCTree_preprop", WeightDictName="I3MCWeightDict_GENIE"
			)	

if options.FLAVOR == 'NuMu':
    ### ADDING PROPAGATOR ###
	tray.AddModule("I3GENIEResultDictToMCTree", "toMcTree", 
			MCTreeName = "I3MCTree_preprop", WeightDictName="I3MCWeightDict_GENIE")	

	# tray.AddModule("Rename", 
	# 		Keys = ["I3MCTree","I3MCTree_preprop"])	
	from icecube import PROPOSAL, sim_services
	propagators = sim_services.I3ParticleTypePropagatorServiceMap()
	
	mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/mediadef')
	
	muMinusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
			mediadef=mediadef,
			cylinderRadius=1200,
			cylinderHeight=1700,
			type=dataclasses.I3Particle.ParticleType.MuMinus)
	muPlusPropagator = PROPOSAL.I3PropagatorServicePROPOSAL(
			mediadef=mediadef,
			cylinderRadius=1200,
			cylinderHeight=1700,
			type=dataclasses.I3Particle.ParticleType.MuPlus)
		
	propagators[dataclasses.I3Particle.ParticleType.MuMinus] = muMinusPropagator
	propagators[dataclasses.I3Particle.ParticleType.MuPlus] = muPlusPropagator
	tray.AddModule('I3PropagatorModule', 'muon_propagator',
			PropagatorServices=propagators,
			RandomService=randomService,
			InputMCTreeName="I3MCTree_preprop",
#			OutputMCTreeName="I3MCTree)
                       OutputMCTreeName="I3MCTree_postprop")	

#######
### Testing: split and pass to root file for checks
#######
# from icecube.tableio import I3TableWriter
# from icecube.rootwriter import I3ROOTTableService
# rootout = options.OUTFILE[:-6] + '.root'
# root_service = I3ROOTTableService(rootout)
# tray.AddModule("I3NullSplitter", "nullsplit")
# tray.AddModule(I3TableWriter,'table_writer',
#                TableService    = [root_service],
#                SubEventStreams = ['nullsplit','in_ice'],
#                BookEverything = True
#                )
# tray.AddModule("Dump", "ShowMeTheMoney")

#SkipKeys = ["I3MCTree_preprop"]

tray.AddModule("I3Writer","writer",
#	SkipKeys= SkipKeys,
    Filename = options.OUTFILE)

tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()
