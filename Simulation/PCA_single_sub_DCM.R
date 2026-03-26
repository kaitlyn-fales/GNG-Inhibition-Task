remove(list=ls())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No data file provided.")
}

data_file <- args[1]

# Load the data
load(data_file)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mcmcse))
suppressPackageStartupMessages(library(momentLS))

######## Change specifications here ###########

# Output specs
file_base <- basename(data_file)

dist <- sub("_(full|thresh)_roi_snr[0-9]+_[0-9]+\\.RData$", "", file_base)

output_dir <- file.path("Output_PCA", dist)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

basename <- tools::file_path_sans_ext(file_base)

# Source functions
source("../Canonical-DCM-Method/canonical_dcm_functions.R")

###############################################

########### Get data ready ####################
# Compile stan program
canonical_dcm = cmdstanr::cmdstan_model("../Canonical-DCM-Method/canonical_dcm.stan")

# Indices of parameters in hypothesis
A_idxs <- matrix(c(1,1,
                   2,1,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence check specs - 95% intervals with 5% tolerance
num_param <- get_num_param(stan_dat)
ess_check <- as.numeric(minESS(num_param, alpha = 0.05, eps = 0.05))
###############################################

########### Initialize sampler ################

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat)

# Run sampler until convergence
dcm_sample(mod = canonical_dcm, 
           data = stan_dat, 
           inits_list = inits_list, 
           output_dir = output_dir,
           basename = basename,
           metric = "dense_e",
           refresh = 100,
           warmup_iter = 5000,
           n_iter_chunk = 1000,
           max_iter = 100000,
           adapt_delta = 0.9,
           seed = 1234,
           chains = 1)
###############################################

