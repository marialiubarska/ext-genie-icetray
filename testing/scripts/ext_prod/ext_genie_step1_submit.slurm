#!/bin/bash
#SBATCH --time=09:59:00
#SBATCH --mem=4G
#SBATCH --output=/home/mliubar/scratch/ext_genie_logs/ext_log_job_%A_%a
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=genie_step1_ext
## #SBATCH --array=0-10
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))

startsecond=$(date +%s)
echo "Start second: " $startsecond 

echo "This is step 1 of neutrino simulations! "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
# generating name of the output file
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="

echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "Starting the singularity job"
singularity exec --bind /cvmfs --bind /project/6008051/mliubar --bind /project/6008051/hignight --bind /home/hignight --bind /home/mliubar --bind /scratch/mliubar --nv /project/6008051/hignight/singularity_images/centos7.img /home/mliubar/Software/genie_workspace/testing/scripts/ext_prod/ext_step1_genie.sh $@
date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
