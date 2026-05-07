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
suppressPackageStartupMessages(library(cdcm))

######## Change specifications here ###########

# Output specs
output_dir <- "Output_OldHyp"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

basename <- tools::file_path_sans_ext(basename(data_file))

###############################################

########### Get data ready ####################
# Compile stan program
canonical_dcm = compile_cdcm()

# Indices of parameters in hypothesis
A_idxs <- matrix(c(1,1,
                   2,1,
                   1,2,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1,
                   2,1,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1), byrow = T, ncol = 2)

idxs <- list(A_idxs = A_idxs,
             B_idxs = B_idxs,
             C_idxs = C_idxs)

# Put data into form for sampler
stan_dat <- get_stan_dat(dat, idxs)

# Convergence check specs - 95% intervals with 5% tolerance
ess_check <- minESS_criterion(stan_dat, alpha = 0.05, eps = 0.05)
###############################################

########### Initialize sampler ################

# Use pathfinder to get good initial values
inits_list <- get_initial_vals(canonical_dcm, stan_dat)

# Run sampler until convergence
results <- dcm_sample(mod = canonical_dcm, 
                      data = stan_dat, 
                      inits_list = inits_list, 
                      output_dir = output_dir,
                      basename = basename,
                      ess_check = ess_check,
                      seed = 1234)
###############################################

