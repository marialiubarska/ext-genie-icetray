#!/bin/bash

module --force purge
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
source /home/mliubar/Software/genie_workspace/setup.sh

genie_app=$MYGENIE/genie/gicevgen
hepev_writer=$MYGENIE/genie/hepevt_writer
hepev_reader=/home/mliubar/Software/genie_workspace/hepevt-reader/new_scripts/step_1_external_genie.py

# filename=$4
# nu_type=$2
# runnumber=$1
# filename=$filename.$nu_type.$runnumber
# nevents=$3
# en_min=10
# en_max=100

# echo "------------------------"
# echo "filename: $filename"
# echo "nevents: $nevents"
# echo "flavor: $nu_type"
# echo "run number: $runnumber"
# echo "energy range: [$en_min,$en_max]"
# echo "------------------------"

OUTDIR=/home/mliubar/Software/genie_workspace/testing/output/ext_prod

# set run_number, neutrino flavor and select energy range
RUNNUM=$1
echo "Run number: " $RUNNUM
FLV=$2
echo "Flavor is: " ${FLV}
E=$3
echo "Energy Range is: " ${E} # A, B, C, D

case ${FLV} in
    NuE)
        NU=12
	RADIUS=250
	LENGTH=500
        case ${E} in
            A)
                NEVENTS=450000
		EMIN=1
		EMAX=4
                ;;
            B)
                NEVENTS=100000
		EMIN=4
		EMAX=12
                ;;
            C)
                NEVENTS=100000
		EMIN=12
		EMAX=9999.9999
                ;;
            D)
                NEVENTS=57500
		EMIN=0.5
		EMAX=1
                ;;
            *)
                echo ${E} " is not an acceptable neutrino energy range (A B C D)"
                exit 2
                ;;
	esac
	;;
    NuMu)
	NU=14
        case ${E} in
            A)
                NEVENTS=408000
		EMIN=1
		EMAX=5
		RADIUS=250
		LENGTH=500
                ;;
            B)
                NEVENTS=440000
		EMIN=5
		EMAX=80
		RADIUS=330
		LENGTH=900
		;;
            C)
                NEVENTS=57500
		EMIN=80
		EMAX=9999.9999
		RADIUS=330
		LENGTH=1500
                ;;
            D)
                NEVENTS=6700
                EMIN=0.5
                EMAX=1
                RADIUS=250
                LENGTH=500
                ;;
	    *)
		echo ${E} " is not an acceptable neutrino energy range (A B C D)"
                exit 2
                ;;
        esac
        ;;
    NuTau)
        NU=16
        case ${E} in
            A)
                NEVENTS=300000
		EMIN=4
		EMAX=10
		RADIUS=250
		LENGTH=500
                ;;
            B)
                NEVENTS=375000
                EMIN=10
                EMAX=30
                RADIUS=250
                LENGTH=500
                ;;
            C)
                NEVENTS=200000
                EMIN=30
                EMAX=9999.9999
                RADIUS=250
                LENGTH=1000
                ;;
            *)
                echo ${E} " is not an acceptable neutrino energy range (A B C D)"
                exit 2
                ;;
        esac
        ;;
	    *)
		echo ${FLV} " is not an acceptable neutrino type (NuE NuMu NuTau)"
		exit 3
		;;

esac

RUNNUM=${NU}${RUNNUM}
echo "Run Number        : "${RUNNUM}
echo "Flavor            : "${FLV}
echo "Energy Range      : "${E}
echo "Number of Events  : "${NEVENTS}
echo "Energy range      : "${EMIN},${EMAX}
echo "Radius            : "${RADIUS}
echo "Length            : "${LENGTH}


filename=${OUTDIR}/${FLV}_${E}_${RUNNUM}_${FILE_NR}_ext_step1
echo "filename          : "${filename}

# gen root file
echo "Generating root file"
$genie_app -r ${RUNNUM} -o ${filename} -n ${NEVENTS} -e ${EMIN},${EMAX} -t 1000080160[0.95],1000010010[0.05] -z 0,180 -a 0,360 -R ${RADIUS} -L ${LENGTH} --nu-type ${FLV} --nu-fraction 0.7 --gamma 2 --cross-sections /home/mliubar/Software/genie_workspace/genie-generator/xsec_splines/GENIE_2_12_8_Water_splines.xml --seed ${RUNNUM} #--force-singleprob-scale 

# conv to hepev
echo "Converting root to hepev"
$hepev_writer -f $filename -o $filename.hepev

# gzip
echo "Gzipping hepev"
gzip -f $filename.hepev
gzip -f $filename.wdict.dat

# conv to i3
echo "converting hepev.gz to i3"
$i3env python $hepev_reader -i $filename.hepev.gz -w $filename.wdict.dat.gz -o $filename.i3 -n $NEVENTS -r $RUNNUM -l $RUNNUM -f $FLV

# deleting intermediate files
echo "Deleting $filename, $filename.$nu_type.status and $filename.hepev.gz"
rm -f $filename $1.$FLV.status $filename.hepev.gz

echo "Done!"
