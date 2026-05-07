# Simulation

This folder contains the simulation studies used throughout Chapter 3 of my dissertation work. The simulations focus on two primary objectives:

1. Practical identifiability of CDCM under different hypothesis structures
2. ROI construction and signal summarization strategies in task-based fMRI analyses

The scripts in this folder generate their own simulation-specific data and output directories during execution.

## Practical Identifiability Simulations

### Simplified Hypothesis Simulation

- `NoAgg_single_sub_DCM.R`  
  Runs the practical identifiability simulation under the simplified hypothesis structure.

- `NoAgg_single_sub_DCM.sh`  
  SLURM submission script for the simplified hypothesis simulation.

- `simulate_data_NoAgg.R`  
  Generates simulated datasets for the simplified hypothesis simulation.

### Original Hypothesis Simulation

- `OldHyp_single_sub_DCM.R`  
  Runs the practical identifiability simulation under the original hypothesis structure.

- `OldHyp_single_sub_DCM.sh`  
  SLURM submission script for the original hypothesis simulation.

- `simulate_data_OldHyp.R`  
  Generates simulated datasets for the original hypothesis simulation.

## ROI Construction Simulations

These simulations compare ROI summarization approaches, including PCA- and mean-based summaries under both full and thresholded ROI construction strategies.

### PCA-Based ROI Summaries

- `PCA_single_sub_DCM.R`  
  Runs the ROI construction simulation using PCA-based ROI summaries.

- `PCA_full_single_sub_DCM.sh`  
  SLURM submission script for the full ROI PCA simulation.

- `PCA_thresh_single_sub_DCM.sh`  
  SLURM submission script for the thresholded ROI PCA simulation.

### Mean-Based ROI Summaries

- `Mean_single_sub_DCM.R`  
  Runs the ROI construction simulation using mean-based ROI summaries.

- `Mean_full_single_sub_DCM.sh`  
  SLURM submission script for the full ROI mean-summary simulation.

- `Mean_thresh_single_sub_DCM.sh`  
  SLURM submission script for the thresholded ROI mean-summary simulation.

### Data Generation

- `simulate_data_Agg.R`  
  Generates simulated datasets for the ROI construction simulations.

## Additional Supporting Scripts

- `hypothesis_assumptions_check.R`  
  Performs assumption and identifiability checks for the simulation study designs.

- `simulation_results.R`  
  Compiles, summarizes, and visualizes simulation outputs across simulation settings.

## Generated Directories

The simulation workflows generate additional directories during execution, including:

- `Data_*`  
  Simulation-specific generated datasets.

- `Output_*`  
  Simulation-specific model outputs and posterior summaries.

- `Metrics/`  
  Additional simulation metrics, summaries, and compiled evaluation results.
