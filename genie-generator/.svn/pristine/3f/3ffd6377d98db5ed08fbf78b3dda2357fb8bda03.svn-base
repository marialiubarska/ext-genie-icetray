#!/bin/bash

echo "GENIE dir is $GENIE"

export GSPLOAD=$MYGENIE/xsec_splines/gxspl-water-v2.6.0.xml

echo "Using xsec file: $GSPLOAD"

export GMSGCONF=$MYGENIE/scripts/Messenger.xml

#export GPRODMODE=YES

date=`date +%Y%m%d_%k%M`
runnum=1
nevts=5000

# generate events: numu x**(-2) flux, between 5 and 200 GeV
# check the range of the cross section xml file, if you go 
# beyond GENIE might crash for unobvious reasons!
$MYGENIE/genie/gicevgen  -r $runnum -n $nevts -e 5,200 -p 14 -f 2 -R 1200 -L 2000 


# convert to hepevt format
echo 
echo "> Converting to hepevt format"

infile=genie_ic.$runnum.ghep.root
outfile=genie_ic.$runnum.hepevt.dat
$MYGENIE/genie/hepevt_writer -f $infile -o $outfile

# gzip the two output files
gzip -f genie_ic.$runnum.hepevt.dat
gzip -f genie_ic.$runnum.wdict.dat

