#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --account=statsresearch_sc_default
#SBATCH --output=Output/output_Insula_thresh_%A_%a.out
#SBATCH --array=1-99

# Get started
echo " "
echo "Job started on `hostname` at `date`"
echo " "

# Environment setup
module purge
module use /gpfs/group/RISE/sw7/modules
module load r/4.2.1

export R_LIBS="~/R_packages"

# Path to file list
FILE_LIST="/storage/work/krf5429/GNG-Inhibition-Task/file_list_thresh_roi.txt"

# Pick the file corresponding to this array task
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

echo "Processing $FILE"

Rscript Insula_single_sub_DCM.R "$FILE"

# Finish up
echo " "
echo "Job Ended at `date`"
echo " "






