#!/bin/bash
date

echo "Starting the neutrino job"
echo "Argument line : " $@

# source /home/hignight/setup_oscnext.sh
echo "TASK ID " $SLURM_ARRAY_TASK_ID
FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
FILE_NR=`printf "%06d\n" $FILE_NR`
echo "Filename ID : " $FILE_NR

echo "Starting cvmfs " 
module --force purge
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
source /home/mliubar/Software/genie_workspace/setup.sh
echo "GENIE: "$GENIE

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
# i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/mliubar/Software/genie_workspace/hepevt-reader/new_scripts/step_1_genie.py
echo "Will use script: " $script
echo  $I3_BUILD " ; ->/genie-icetray/resources/PDFsets"

RUNNUM=$1
echo "Run number: " $RUNNUM
FLV=$2
echo "Flavor is: " ${FLV}
E=$3
echo "Energy Range is: " ${E}

OUTDIR=/home/mliubar/Software/genie_workspace/testing/output/off_prod

GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
#Get set variables 
case ${RUNNUM} in
    0000)
	
	;;
    *)
	echo "No configuration for " $RUNNUM "... exiting"
	exit 1
	;;
esac

#Get FLV setup
case ${FLV} in
    NuE)
        NU=12
        case ${E} in
            A)
                NEVENTS=450000
                ;;
            B)
                NEVENTS=100000
                ;;
            C)
                NEVENTS=100000
                ;;
            D)
                NEVENTS=57500
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
                ;;
            B)
                NEVENTS=440000
                ;;
            C)
                NEVENTS=57500
                ;;
            D)
                NEVENTS=6700
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
                ;;
            B)
                NEVENTS=375000
                ;;
            C)
                NEVENTS=200000
                ;;
            D)
                NEVENTS=26000
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


OUTNAME=${FLV}_${E}_${RUNNUM}_${FILE_NR}_step1.zst

echo "OUTFILE NAME : " ${OUTNAME}
$i3env python $script -o ${OUTDIR}/${FLV}/${OUTNAME} -l ${FILE_NR}  -r ${RUNNUM} -n ${NEVENTS} -f ${FLV} --energy-range ${E} 

date


###########################
####     Exit codes    ####
#### 1 - No set number ####
#### 2 - No E code     ####
#### 3 - No FLV code   ####
###########################

