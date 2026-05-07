# Group_Level_Run

This folder contains the run-specific CDCM group-level analyses. Unlike the multirun hierarchical model implemented elsewhere in the repository, the analyses in this folder model each run independently using posterior summaries obtained from the single-subject CDCM analyses.

These analyses use the `cdcm` package group-level meta-analysis framework to estimate run-level population effects separately for each task run.

## Files

- `group_analysis_CDCM.R`  
  Runs the CDCM group-level meta-analysis model for a single run using single-subject posterior summaries as inputs.

- `group_results_CDCM.R`  
  Compiles and summarizes posterior results from the run-level CDCM group analyses.

- `run1_group_analysis_CDCM.sh`  
  SLURM submission script for the Run 1 CDCM group-level analysis.

- `run2_group_analysis_CDCM.sh`  
  SLURM submission script for the Run 2 CDCM group-level analysis.

- `run3_group_analysis_CDCM.sh`  
  SLURM submission script for the Run 3 CDCM group-level analysis.

## Subfolders

### `Results/`

Placeholder directory for run-specific group-level model outputs, posterior summaries, diagnostics, and compiled result tables.
