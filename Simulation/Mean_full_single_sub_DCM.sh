#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --output=Output_Mean/output_full_%A_%a.out
#SBATCH --array=1-400

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

snr=$(( (idx-1)/50 + 1 ))
rep=$(( (idx-1)%50 + 1 ))

file="/storage/work/krf5429/GNG-Inhibition-Task/Simulation/Data_Mean/full_roi_snr${snr}_${rep}.RData"

echo "Array ID: $idx"
echo "SNR: $snr  Replicate: $rep"
echo "Processing $file"

Rscript Mean_single_sub_DCM.R "$file"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "






