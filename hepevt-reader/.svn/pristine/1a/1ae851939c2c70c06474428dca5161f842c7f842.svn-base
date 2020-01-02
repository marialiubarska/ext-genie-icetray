#!/usr/bin/env python
from I3Tray import *

from os.path import expandvars

import os
import sys

if(os.getenv('I3_BUILD') == None):
    print "I3_BUILD not set."
    sys.exit()

data_dir=sys.argv[1]
base_name=sys.argv[2]
job=sys.argv[3]
runnum=int(job)

events=5000
seed=runnum

fname_noext=base_name+"_"+job
infile=fname_noext+".p2.i3.gz"
outfile=fname_noext+".p3.i3.gz"

## gcd file
gcd_file="/user/mdier/scratch/icecube/genie/data_prod/GCD/NEW_GeoCalibDetectorStatus_IC80_DC6.54655.i3.gz"

load("libdataclasses")
load("libphys-services")
load("libdataio")
load("libsim-services")
load("libpmt-simulator")
load("libnoise-generator")
load("libDOMsimulator")
load("libphotonics-service")

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

tray.AddModule("I3Muxer","muxer")

tray.AddModule("I3NoiseGeneratorModule","noise")(
    ("InIce",True),
    ("StartWindow",10.*I3Units.microsecond),
    ("EndWindow",10.*I3Units.microsecond),
    ("Rate",650.*I3Units.hertz),
    ("TWR",False),
    ("DeepCoreFactor",1.54),
    ("Icetop",False),
    ("InputHitSeriesMapName","MCHitSeriesMap"),
    )

tray.AddModule("I3PMTSimulator","pmt")

tray.AddModule("I3DOMsimulator","domsimulator")

#tray.AddModule("I3LETrigger","trigger2")(
#   ("FilterMode",False),
#   ("triggersource",0),
#   ("TriggerName","le_trig"),
#   ("Threshold",2),
#   ("TimeWindow",5000*I3Units.ns),
#   )

#tray.AddModule("I3SMTrigger","smtrigger")(
# #  ("TriggerName","I3Triggers")
#   )

#tray.AddModule("I3GlobalTriggerSim","globaltrigger")(
##("I3TWRReadoutOffset",0.*I3Units.ns),
#("I3ReadoutWindow",20000.*I3Units.ns),
##("I3TWRReadoutWindow",50000.*I3Units.ns),
#)

tray.AddModule("I3EventCounter","counter")(
("physicscountername","IC86 Triggered Events"),
)

tray.AddModule("I3Writer","writerend")(
    ("streams",[icetray.I3Frame.Physics]),
    ("filename", outfile)
    )

#tray.AddModule("Dump","dump")

tray.AddModule("TrashCan","trashcan")

tray.Execute(events+4)
tray.Finish()


