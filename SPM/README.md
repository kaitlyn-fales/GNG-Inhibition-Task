# SPM

This folder contains the SPM-based DCM analyses used as a comparison framework and sanity check against the CDCM analyses. The workflows include single-subject SPM DCM estimation, preprocessing utilities, and group-level analyses based on the SPM outputs.

The analyses in this folder were implemented using MATLAB and SPM12.

## Files

### Single-Subject SPM DCM Scripts

- `MFG_full_single_sub.m`  
  Runs single-subject SPM DCM estimation for the MFG-centered network using the full ROI approach.

- `MFG_thresh_single_sub.m`  
  Runs single-subject SPM DCM estimation for the MFG-centered network using the thresholded ROI approach.

- `Insula_full_single_sub.m`  
  Runs single-subject SPM DCM estimation for the Insula-centered network using the full ROI approach.

- `Insula_thresh_single_sub.m`  
  Runs single-subject SPM DCM estimation for the Insula-centered network using the thresholded ROI approach.

### HPC Submission Scripts

- `*.sh`  
  SLURM submission scripts for running the SPM analyses on an HPC environment.

### Group-Level Analysis Scripts

- `group_analysis_SPM.R`  
  Runs the group-level analyses based on the estimated SPM DCM parameters.

- `group_results_SPM.R`  
  Compiles and summarizes the group-level SPM analysis outputs.

### Preprocessing and Supporting Files

- `preprocess_SPM.R`  
  Prepares and organizes SPM output files for downstream analyses.

- `sub_run.txt`  
  Supporting batch-processing file used for SPM workflow execution.

## Subfolders

### `Output/`

Placeholder directory for SPM model outputs, estimated DCM objects, compiled group-level summaries, and diagnostic results.
