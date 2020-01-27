#!/bin/bash

genie_app=$MYGENIE/genie/gicevgen
hepev_writer=$MYGENIE/genie/hepevt_writer
hepev_reader=/home/mliubar/Software/genie_workspace/hepevt-reader/new_scripts/step_1_external_genie.py

filename=$1
nu_type=$2
runnumber=$3
filename=$filename.$nu_type.$runnumber
nevents=$4
en_min=1
en_max=2

echo "------------------------"
echo "filename: $filename"
echo "nevents: $nevents"
echo "flavor: $nu_type"
echo "run number: $runnumber"
echo "energy range: [$en_min,$en_max]"
echo "------------------------"

# gen root file
echo "Generating root file"
$genie_app -r $runnumber -o $filename -n $nevents -e $en_min,$en_max -t 1000080160[0.3333],1000010010[0.6666] -z 0,180 -a 0,360 -R 400 -L 900 --nu-type $nu_type --nu-fraction 0.7 --gamma 2 --cross-sections /home/mliubar/Software/genie_workspace/genie-generator/xsec_splines/GENIE_2_12_8_Water_splines.xml --seed $runnumber #--force-singleprob-scale 

# conv to hepev
echo "Converting root to hepev"
$hepev_writer -f $filename -o $filename.hepev

# gzip
echo "Gzipping hepev"
gzip -f $filename.hepev
gzip -f $filename.wdict.dat

# conv to i3
echo "converting hepev.gz to i3"
$i3env python $hepev_reader -i $filename.hepev.gz -w $filename.wdict.dat.gz -o $filename.i3 -n $nevents -r $runnumber -l $runnumber -f $nu_type 

# deleting intermediate files
echo "Deleting $filename, $filename.$nu_type.status, $filename.hepev.gz and $filename.wdict.dat.gz"
rm -f $filename $1.$nu_type.status $filename.hepev.gz $filename.wdict.dat.gz

echo "Done!"
