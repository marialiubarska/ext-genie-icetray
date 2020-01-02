#!/usr/bin/env python
from I3Tray import *

from os.path import expandvars,isfile

import os
import sys

if(os.getenv('I3_BUILD') == None):
    print "I3_BUILD not set."
    sys.exit()

#photonics = True
#if(os.getenv('PHOTONICS_TABLES') == None):
#    photonics = False

data_dir=sys.argv[1]
base_name=sys.argv[2]
job=sys.argv[3]
runnum=int(job)

events=5000
seed=runnum

fname_noext=base_name+"."+str(runnum)
infile=data_dir+"/"+fname_noext+".hepevt"
infile_gz=data_dir+"/"+fname_noext+".hepevt.dat.gz"
wd_infile_gz=data_dir+"/"+fname_noext+".wdict.dat.gz"

outfile=base_name+"_"+job+".p1.i3.gz"

print infile
print infile_gz
print wd_infile_gz
print outfile

## check for the existance of the wdict file

addWeights=True
if not isfile(wd_infile_gz):
    addWeights=False
    print "Could not find file containing weights, will not be added"

## variables for geometry selector
ic22 = "21,29,30,38,39,40,49,50,59,58,67,66,74,73,65,72,78,48,57,47,46,56"
ic2008= "63,64,55,71,70,76,77,75,69,60,68,61,62,52,44,53,54,45"
ic2009 = "2,3,4,5,6,10,11,12,13,17,18,19,20,26,27,28,36,37,83"
stationsToExclude = "1:80"

## gcd file
gcd_file="/user/mdier/scratch/icecube/genie-generator/data_prod/GCD/NEW_GeoCalibDetectorStatus_IC80_DC6.54655.i3.gz"


load("libdataclasses")
load("libphys-services")
load("libdataio")
load("libsim-services")
load("libmmc-icetray")

tray = I3Tray()

tray.AddService("I3JavaVMFactory","java")(
    ("Options",[expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar")])
    )

tray.AddService("I3FileOMKey2MBIDFactory","omkey2mbid")(
    ("INFILE",expandvars("$I3_BUILD/phys-services/resources/doms.txt"))
    )

tray.AddService("I3ReaderServiceFactory","gcd")(
    ("filename",gcd_file),
    ("OmitEvent",True)
    )

tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom")(
	("NStreams",2),
	("Seed",seed),
        ("StreamNum",1)
        )

tray.AddService("I3MCTimeGeneratorServiceFactory","time-gen")(
        ("Year",2008),
        ("DAQTime",163165000000000000),
        ("runnumber",runnum)
        )

tray.AddService("I3GeometrySelectorServiceFactory","geo-selector")(
     ("StringsToUse","1:86"),
    ("StationsToExclude",stationsToExclude),
    ("GeoSelectorName","IC80_DC6-Geo")
    )

tray.AddModule("I3Muxer","muxer")(
    ("GeometryService","IC80_DC6-Geo"),
    ("eventservice","I3EventService")
    )

tray.AddModule("I3EventCounter","counter3")(
    ("nevents",events),
    ("physicscountername","Generated Events"),
    )

from hepevt_file_reader import HEPEvtFileReaderModule

tray.AddModule(HEPEvtFileReaderModule,"hepevt-reader")(
    ("Filename",infile_gz),
    ("AttachAllToLepton",True),
    ("AddHadronic",True)
    )

if addWeights:
    from weightdict_file_reader import WDictFileReaderModule
    tray.AddModule(WDictFileReaderModule,"wdict-reader")(
        ("Filename",wd_infile_gz)
        )

# Skip this for now, generating low energy events
#
#tray.AddModule("I3CascadeMCModule","cmcprop")(
#    ("energythresholdsplit",1000),
#    ("inputmctree","I3MCTree"),
#    ("energythresholdsimulation",1000000),
#    ("splitwidth",3),
#    ("outputmctree","I3MCTreeCMC")
#    )

mmcmode = -3
mmcopts = '-seed=1 -radius=1200 -length=1700'
tray.AddModule("I3PropagatorMMC","mmc")(
    ('rerr',"/dev/null"),
    ('mode', mmcmode),
    ('opts', mmcopts)
    )

tray.AddModule("I3Writer","writer")(
    ("streams",[icetray.I3Frame.Physics]),
    ("filename", outfile)
    )

#tray.AddModule("Dump", "the dump")
tray.AddModule("TrashCan", "the can")

tray.Execute(events+4)
tray.Finish()
