#!/bin/bash
#
#  Usage: 
# subm_gicevgen.sh pdg_id runnum_start nr_runs emin emax expon
if [ $# -ne 6 ]
then
    echo "Usage: `basename $0` pdg_id runnum_start nr_runs Emin(GeV) Emax(GeV) exponent"
    exit 1
fi

export flavor=$1
runnum_start=$2
nr_runs=$3
export emin=$4
export emax=$5
export expon=$6
#if [ $emin > $emax ]
#then
#    echo "Minumun energy (${emin} GeV) larger than maximum energy (${emax} GeV)"
#    exit 1
#fi

runnum_end=`expr ${runnum_start} + ${nr_runs}`
runnum_end=`expr ${runnum_end} - 1`

flavorset=0
if [ $flavor -eq 12 ]; then
    flavorset=1
elif [ $flavor -eq -12 ]; then
    flavorset=1
elif [ $flavor -eq 14 ]; then
    flavorset=1
elif [ $flavor -eq -14 ]; then
    flavorset=1
elif [ $flavor -eq 16 ]; then
    flavorset=1
elif [ $flavor -eq -16 ]; then
    flavorset=1
fi

if [ $flavorset -eq 0 ]; then
    echo "Unknown neutrino flavor id: $flavor"
    exit
fi

genie_job=$MYGENIE/pbs/gicevgen_pbs.job

for i in `seq ${runnum_start} ${runnum_end}`; do
    export runnum=$i
    varlist="runnum,flavor,emin,emax,expon"

    qsub -d `pwd` -q localcluster -l walltime=2:00:00 -v $varlist $genie_job
done

