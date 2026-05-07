# Group_Mdl

This folder contains the Stan implementations of the hierarchical multirun CDCM group-level models. These models jointly analyze multiple runs within a unified hierarchical framework rather than modeling each run independently.

The multirun model extends the run-specific meta-analysis framework by incorporating run-level structure directly within a single hierarchical model.

## Files

- `meta_analysis_multirun.stan`  
  Stan implementation of the hierarchical multirun CDCM meta-analysis model.

- `meta_analysis_multirun`  
  Compiled CmdStan executable corresponding to the multirun Stan model.

## Notes

The multirun model is designed to jointly estimate population-level effects across all task runs while accounting for within-subject dependence across runs.
