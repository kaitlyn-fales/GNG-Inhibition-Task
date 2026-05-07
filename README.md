# GNG-Inhibition-Task

This repository contains the code, simulation framework, and analysis pipeline for the smoking abstinence Go/No-Go (GNG) inhibitory control fMRI analyses used in my dissertation work on Canonical Dynamic Causal Modeling (CDCM). The repository includes:

* Single-subject CDCM analyses
* Group-level CDCM meta-analysis models
* SPM-based comparison analyses
* Simulation studies on practical identifiability and ROI construction
* Preprocessing and ROI generation scripts

A major goal of this repository is to provide a reproducible implementation of the CDCM framework for task-based effective connectivity analysis.

---

# Repository Overview

## Main Analysis Components

### `Analysis/`

Follow-up CDCM diagnostics, posterior analyses, and downstream group-level analyses. This folder contains scripts used after the primary single-subject CDCM estimation pipeline, including analyses examining parameter behavior, sensitivity analyses, and additional model summaries.

### `Group_Level_Run/`

Run-specific CDCM group-level analyses using the `cdcm` package meta-analysis framework. Each run is modeled separately using posterior summaries from the single-subject CDCM analyses.

### `Group_Mdl/`

Stan implementations of the hierarchical multirun CDCM meta-analysis models. These models jointly incorporate multiple runs within a unified hierarchical framework.

### `SPM/`

SPM-based DCM analyses used as a comparison and sanity check against CDCM. This folder contains:

* Single-subject SPM DCM estimation scripts
* Group-level SPM analyses
* Multirun SPM models
* Supporting preprocessing and submission scripts

### `Simulation/`

Simulation studies used throughout Chapter 3 of my dissertation work, including:

1. Practical identifiability simulations for CDCM
2. ROI construction simulations examining full versus thresholded ROI approaches

## Key Scripts and Files

### Main Single-Subject CDCM Scripts

- `MFG_single_sub_DCM.R`  
  Runs single-subject CDCM estimation for the MFG-centered network configuration.

- `Insula_single_sub_DCM.R`  
  Runs single-subject CDCM estimation for the Insula-centered network configuration.

### Preprocessing and ROI Scripts

- `preprocess_data_ses1_run_full_roi.R`  
  Preprocesses BOLD time series using the full ROI approach.

- `preprocess_data_ses1_run_thresh_roi.R`  
  Preprocesses BOLD time series using the thresholded ROI approach.

- `VOI_mask_generate.R`  
  Generates masks and supporting files for ROI extraction workflows.

- `mask.nii`  
  NIfTI mask file used for ROI generation and preprocessing.

### Supporting Workflow Files

- `*.sh`  
  SLURM submission scripts for HPC-based model estimation.

- `*.txt`  
  Batch-processing file lists for SLURM jobs.

---

# Additional Repository Structure

### `Data/`

Placeholder directory for raw and processed data inputs.

Data are not distributed with this repository due to size and/or data-sharing restrictions.

### `Output/`

Placeholder directory for generated model outputs, posterior draws, figures, diagnostics, and summaries.

---

# CDCM Package Integration

This repository relies heavily on the companion `cdcm` R package, which implements the Canonical Dynamic Causal Modeling framework.

The `cdcm` package repository is available here: [https://github.com/kaitlyn-fales/cdcm](https://github.com/kaitlyn-fales/cdcm)

Before running the analyses in this repository, install the `cdcm` package and its dependencies.

---

# Software Requirements

The primary analyses were developed using:

* R (>= 4.2)
* CmdStan / CmdStanR
* Stan
* SPM12 (MATLAB)

Some analyses additionally require:

* MATLAB
* SPM12
* SLURM-based HPC environments for large-scale estimation

---

# Data Notes

This repository was developed for analyses of smoking abstinence Go/No-Go fMRI data. Because the underlying imaging data are not publicly distributed through this repository, users must supply their own data inputs and adapt paths accordingly.

Several scripts contain HPC-specific paths and submission configurations that may need to be modified for local environments.

---

# Reproducibility Notes

The repository preserves the original organizational structure used during dissertation development, including:

* Separate CDCM and SPM workflows
* Standalone simulation pipelines
* SLURM submission scripts
* Group-level hierarchical model implementations

Many scripts are designed to run independently once the required data objects and directory structure are in place.

---

# Citation

If you use this repository or the CDCM framework in your own work, please cite the associated dissertation and/or methodological manuscripts when available.



