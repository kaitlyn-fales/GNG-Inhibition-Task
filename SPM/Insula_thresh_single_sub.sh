#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=open
#SBATCH --output=Output/Insula_thresh_%A_%a.out
#SBATCH --array=1-99

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load matlab/R2023a 

# --- Get subject ID from file ---
sub=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sub_run.txt)
echo "Running subject: $sub"

matlab -batch "sub='$sub'; run('Insula_thresh_single_sub.m')"


# Finish up
echo " "
echo "Job Ended at `date`"
echo " "






