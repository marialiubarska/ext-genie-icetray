#!/usr/bin/env python
from I3Tray import *

from os.path import expandvars

import os
import sys

if(os.getenv('I3_BUILD') == None):
    print "I3_BUILD not set."
    sys.exit()

photonics = True
if(os.getenv('PHOTONICS_TABLES') == None):
    photonics = False

data_dir=sys.argv[1]
base_name=sys.argv[2]
job=sys.argv[3]
runnum=int(job)

events=5000
seed=runnum

fname_noext=base_name+"_"+job
infile=fname_noext+".p1.i3.gz"
outfile=fname_noext+".p2.i3.gz"

## variables for geometry selector
#ic22 = "21,29,30,38,39,40,49,50,59,58,67,66,74,73,65,72,78,48,57,47,46,56"
#ic2008= "63,64,55,71,70,76,77,75,69,60,68,61,62,52,44,53,54,45"
#ic2009 = "2,3,4,5,6,10,11,12,13,17,18,19,20,26,27,28,36,37,83"
#stationsToExclude = "1:80"

## gcd file
gcd_file="/user/mdier/scratch/icecube/genie/data_prod/GCD/NEW_GeoCalibDetectorStatus_IC80_DC6.54655.i3.gz"

load("libdataclasses")
load("libphys-services")
load("libdataio")
load("libsim-services")
load("libicepick")
load("libphotonics-service")
load("libhit-maker")

tray = I3Tray()

tray.AddService("I3ReaderServiceFactory","reader")(
    ("filenamelist",
     [gcd_file,infile])
    )

tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom")(
	("NStreams",2),
	("Seed",seed),
        ("StreamNum",1)
        )


if(photonics):
    tablesDir = expandvars("$PHOTONICS_TABLES")
    driverFilePath = tablesDir + "/listfiles_AHA07v2h2/I3Coord_I3Span_z80_a60/"
    level1Driver = "level1_shower.list"
    level2Driver = "level2_muon.list"

    tray.AddService("I3PhotonicsServiceFactory","photonics")(
        ("PhotonicsTopLevelDirectory",tablesDir),
        ("DriverFileDirectory",driverFilePath),
        ("PhotonicsLevel1DriverFile",level1Driver),
        ("PhotonicsLevel2DriverFile",level2Driver)
        )
else:
    print "******************************************************"
    print "*** Environment variable PHOTONICS_TABLES not set. ***"
    print "******************************************************"
    sys.exit()

tray.AddModule("I3Muxer","muxer")

#tray.AddModule("I3IcePickModule<I3SkipNEventFilter>", "skip")(
#    ("discardEvents", True),
#    ("SkipNevents", 1108),
#    ("NeventStopick", 20)
#)

tray.AddModule("I3HitMakerModule","hit-maker")(
    ("EnableBinning",True),
    ("MaxPEs",500000),
    ("NBins",300),
    ("DeepCoreFactor",1.25)
    )

tray.AddModule("I3Writer","writer")(
    ("streams",[icetray.I3Frame.Physics]),
    ("filename", outfile)
    )

#tray.AddModule("Dump", "the dump")
tray.AddModule("TrashCan", "the can")

tray.Execute()
tray.Finish()


