#!/bin/bash

export GENIE=/home/hignight/work/simulation/GENIE/R-2_12_8
export PYTHIA6=/home/hignight/work/simulation/GENIE/pythia6_4_28/lib
export LHAPDF=/home/hignight/work/simulation/GENIE/lhapdf-5.9.1
export LHAPATH=$LHAPDF
export LOG4CPP=/cvmfs/icecube.opensciencegrid.org/py2-v3.1.0/RHEL_7_x86_64/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:$PYTHIA6/lib:$LHAPDF/lib:$LOG4CPP/lib:$GENIE/lib
export PATH=$PATH:$ROOTSYS/bin:$GENIE/bin:$LHAPDF/bin

export MYGENIE=/home/mliubar/Software/genie_workspace/genie-generator
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MYGENIE/lib

export PYTHONPATH=$PYTHONPATH:/home/mliubar/Software/genie_workspace/hepevt-reader/lib

i3env=~/my_icetray_Wcut0/full_build/env-shell.sh

testgen=/home/mliubar/Software/genie_workspace/testing/scripts/ext_prod/gen_test.sh
testgen2=/home/mliubar/Software/genie_workspace/testing/scripts/ext_prod/gen_test_notdelete.sh

i3plot=/home/mliubar/Software/genie_workspace/testing/scripts/plot_MCWeightDict.py
