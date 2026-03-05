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
# Denote which activation region we are using - either MFG or Insula
region <- "MFG"
dat$y_obs <- dat$y_obs[,c(region,"Precuneus/PCC")]

# Subject specs
file_name <- basename(data_file)

# Match sub ID and run number, ignoring anything between run number and .RData
matches <- str_match(file_name, "sub-(\\d+)_ses1_run(\\d+)_(full|thresh)_roi\\.RData")

sub <- matches[2]        # subject number
run <- as.numeric(matches[3])  # run number
type <- matches[4]       # "full" or "thresh"

ses <- 1

# Output specs
output_dir <- "Output"
basename <- paste0("sub-",sub,"_ses",ses,"_run",run,"_",region,"_",type)

# Source functions
source("Canonical-DCM-Method/canonical_dcm_functions.R")

###############################################

paste0("DCM for Subject ",sub,", Session ", ses, ", Run ", run, ", ",type," ROI, using ",region)

########### Get data ready ####################
# Compile stan program
canonical_dcm = cmdstanr::cmdstan_model("Canonical-DCM-Method/canonical_dcm.stan")

# Indices of parameters in hypothesis
A_idxs <- matrix(c(1,1,
                   2,1,
                   1,2,
                   2,2), byrow = T, ncol = 2)
B_idxs <- matrix(c(2,1,1,
                   2,2,1,
                   2,1,2), byrow = T, ncol = 3)
C_idxs <- matrix(c(1,1,
                   2,1), byrow = T, ncol = 2)

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

