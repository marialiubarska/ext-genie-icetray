#!/bin/bash

echo "GENIE dir is $GENIE"

export GSPLOAD=$MYGENIE/xsec_splines/gxspl-water-v2.6.0.xml

echo "Using xsec file: $GSPLOAD"

export GMSGCONF=$MYGENIE/scripts/Messenger.xml

#export GPRODMODE=YES

outfile=xsec.root
maxe=200

probes=( 12 -12 14 -14 16 -16 )
targets=( 1000010010 1000080160 )

for tgt in ${targets[@]}
do
  for prb in ${probes[@]}
  do 
    echo 
    echo ">>> Starting: "
    echo ">>>    gspl2root -f $GSPLOAD -o $outfile -e $maxe -p $prb -t $tgt "
    echo ">>> "
    echo 
    gspl2root -f $GSPLOAD -o $outfile -e $maxe -p $prb -t $tgt
  done
done



