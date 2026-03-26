#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --output=Output_PCA/output_thresh_%A_%a.out
#SBATCH --array=1-800

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"

idx=$SLURM_ARRAY_TASK_ID

# Map array index to distribution, SNR, replicate
dist_idx=$(( (idx-1)/400 + 1 ))
within_dist=$(( (idx-1)%400 + 1 ))

snr=$(( (within_dist-1)/50 + 1 ))
rep=$(( (within_dist-1)%50 + 1 ))

if [ "$dist_idx" -eq 1 ]; then
  dist="beta_sym"
elif [ "$dist_idx" -eq 2 ]; then
  dist="beta_u"
else
  echo "Invalid distribution index: $dist_idx"
  exit 1
fi

file="/storage/work/krf5429/GNG-Inhibition-Task/Simulation/Data_PCA/${dist}_thresh_roi_snr${snr}_${rep}.RData"

echo "Array ID: $idx"
echo "Distribution: $dist"
echo "SNR: $snr  Replicate: $rep"
echo "Processing $file"

Rscript PCA_single_sub_DCM.R "$file"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "