# PEB for combining subjects into group DCM
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(tidyverse))

# Get environment variables from Slurm
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Define all combinations of region and VOI type
region <- c("MFG","Insula")
VOI_type <- c("full","thresh")

# Expand grid
conditions <- expand.grid(region = region, type = VOI_type, stringsAsFactors = FALSE)

# Pick the corresponding row
region_condition <- conditions$region[task_id]
type_condition  <- conditions$type[task_id]

cat("Running group analysis for Region =", region_condition, "and ROI =", type_condition, "\n")

# Load diagnostics df to use indices for looping through
load("diagnostics_compilation.RData")

# Filter for one phase and mask condition - do for all four conditions
df <- filter(diagnostics_df, region == region_condition & VOI_type == type_condition)

# Filter out subject 110 and 111 - no covariate info right now
df <- filter(df, !(subject %in% c(110,111)))

# Add new subject index for stan
df$subj_idx_run <- rep(1:length(unique(df$subject)),each = 3)

# Set up MCMC posterior outputs for meta-analysis
###############################################
y_list <- list()
S_list <- list()

for (i in 1:nrow(df)){
  
  # Load your R posterior draws
  load(paste0("../Output/sub-", df$subject[i],
              "_ses1_run", df$run[i], "_", df$region[i],"_", df$VOI_type[i], "_draws.RData"))
  
  # Get rid of unnecessary columns
  draws <- suppressWarnings(draws_df[,c(4:9)])
  
  # Transform the diag(A) draws according to our reparameterization
  # nu_A[1] and nu_A[3]
  draws[,c(1,3)] <- -0.5*exp(draws[,c(1,3)])
  
  # Extract posterior means
  post_means <- summarize_draws(draws, mean)[,2]
  y_list[[i]] <- as.numeric(as.matrix(post_means))
  
  # Extract posterior covariance
  S_list[[i]] <- cov(draws)
  
}

# Subject-level covariates
X_subj1 <- read.csv("../Data/participants.csv")
X_subj2 <- read.csv("../Data/participants_ftnd.csv")

X_subj <- merge(X_subj1,X_subj2,by = "SubjectID") %>% filter(SubjectID != 176)
X_subj$subj_idx <- c(1:nrow(X_subj))

X_subj <- X_subj %>%
  arrange(subj_idx) %>%
  mutate(
    age_std = scale(Age),
    sex_code = ifelse(Sex %in% c("M", 1), 0.5, -0.5),  # effect coding
    ftnd_std = scale(FTND)
  ) %>%
  select(age_std, sex_code, ftnd_std) %>%
  as.matrix()

q_subj <- ncol(X_subj)

# Run-level covariates 
run_cov <- read.csv("../Data/participants_acc.csv")
run_cov <- run_cov %>% filter(SubjectID != 176)
run_cov$subj_idx <- rep(1:length(unique(df$subject)),each = 3)

X_run <- run_cov %>%
  arrange(subj_idx, Run) %>%
  mutate(accuracy_std = scale(accuracy_avg)) %>%
  select(accuracy_std) %>%
  as.matrix()

q_run <- ncol(X_run)

# Get data in proper format
K <- nrow(X_subj)         # number of subjects
p <- 6                     # number of parameters per subject
subj_idx_run <- run_cov$subj_idx  # length N
N <- nrow(X_run) 

meta_data <- list(
  N      = N,         # total runs
  K      = K,         # total subjects
  p      = p,  # number of parameters
  q_subj = q_subj,
  q_run  = q_run,
  y      = y_list,
  S      = S_list,
  X_subj = X_subj,
  X_run  = X_run,
  subj   = subj_idx_run         
)

# Compile stan model
mod <- cmdstan_model("../Group_Mdl/meta_analysis_multirun.stan")  

init_fun <- function() {
  list(
    alpha   = rep(0, p),
    B_subj  = matrix(0, p, q_subj),   # subject-level slopes
    B_run   = matrix(0, p, q_run),    # run-level slopes
    Ltau    = diag(p),
    tau     = rep(0.1, p),
    Lcorr   = diag(p)
  )
}

# Run pathfinder to get good starting MCMC values
pf <- mod$pathfinder(data = meta_data, init = init_fun, num_paths = 1)

# Number of chains
num_chains <- 5

# Column names
param_names <- colnames(pf$draws())

# Extract a single Pathfinder draw (e.g., the last one)
draw_i <- tail(pf$draws(), 1)

# alpha: vector of length p
alpha_vals <- as.numeric(draw_i[grep("^alpha\\[", param_names)])

# B_subj
B_subj_vals <- as.numeric(draw_i[grep("^B_subj\\[", param_names)])
B_subj_mat <- matrix(B_subj_vals, nrow = p, ncol = q_subj, byrow = TRUE)

# B_run
B_run_vals <- as.numeric(draw_i[grep("^B_run\\[", param_names)])
B_run_mat <- matrix(B_run_vals, nrow = p, ncol = q_run, byrow = TRUE)

# tau: vector of length p
tau_vals <- as.numeric(draw_i[grep("^tau\\[", param_names)])

# Lcorr: matrix p x p
Lcorr_vals <- as.numeric(draw_i[grep("^Lcorr\\[", param_names)])
Lcorr_mat <- matrix(Lcorr_vals, nrow = p, ncol = p, byrow = TRUE)

# Create a single chain init list
single_init <- list(
  alpha   = alpha_vals,
  B_subj  = B_subj_mat,
  B_run   = B_run_mat,
  tau     = tau_vals,
  Lcorr   = Lcorr_mat
)

# Replicate the same init for all chains
init_list <- rep(list(single_init), num_chains)

# Fit model
fit <- mod$sample(
  data = meta_data,
  chains = 5,
  parallel_chains = 5,
  init = init_list,
  iter_warmup = 1000,
  iter_sampling = 5000,
  adapt_delta = 0.8,
  max_treedepth = 10,
  seed = 1234
)

# Extract draws
draws <- as_draws_df(fit$draws())

# Save
save(draws, file = paste0("Results/",region_condition,"_",type_condition,".RData"))


