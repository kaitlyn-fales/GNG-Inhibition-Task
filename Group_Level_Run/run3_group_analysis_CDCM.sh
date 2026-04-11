#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --account=statsresearch_sc_default
#SBATCH --output=output_%A_%a.out
#SBATCH --array=1-4

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"
export RUN=3

Rscript group_analysis_CDCM.R 

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "






