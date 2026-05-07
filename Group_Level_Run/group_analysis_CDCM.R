# PEB for combining subjects into group DCM
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(tidyverse))

# Get environment variables from Slurm
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
run_id <- as.integer(Sys.getenv("RUN"))

# Define all combinations of region and VOI type
region <- c("MFG","Insula")
VOI_type <- c("full","thresh")

# Expand grid
conditions <- expand.grid(region = region, type = VOI_type, stringsAsFactors = FALSE)

# Pick the corresponding row
region_condition <- conditions$region[task_id]
type_condition  <- conditions$type[task_id]

cat("Running group analysis for Region =", region_condition, "and ROI =", type_condition, "and Run =", run_id, "\n")

# Load diagnostics df to use indices for looping through
load("../Analysis/diagnostics_compilation.RData")

# Filter for one phase and mask condition - do for all four conditions
df <- filter(diagnostics_df, region == region_condition & VOI_type == type_condition & run == run_id)

# Filter out subject 110 and 111 - no covariate info right now
df <- filter(df, !(subject %in% c(110,111)))

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

# Covariates
X_subj1 <- read.csv("../Data/participants.csv")
X_subj2 <- read.csv("../Data/participants_ftnd.csv")

# Run-level covariates 
run_cov <- read.csv("../Data/participants_acc.csv")
run_cov <- run_cov %>% filter(Run == run_id)

X_subj <- merge(X_subj1,X_subj2,by = "SubjectID") %>% filter(SubjectID != 176)
X_subj <- merge(X_subj,run_cov,by = "SubjectID")
X_subj$subj_idx <- c(1:nrow(X_subj))

X_subj <- X_subj %>%
  arrange(subj_idx) %>%
  mutate(
    age_std = scale(Age),
    sex_code = ifelse(Sex %in% c("M", 1), 0.5, -0.5),  # effect coding
    ftnd_std = scale(FTND),
    accuracy_std = scale(accuracy_avg)
  ) %>%
  select(age_std, sex_code, ftnd_std, accuracy_std) %>%
  as.matrix()

# Get data in proper format
K <- length(y_list)         # number of studies
p <- 6                     # number of parameters per subject
q <- ncol(X_subj)                      # number of covariates

meta_data <- list(
  K = K,
  p = p,
  q = q,
  y = y_list,
  S = S_list,
  X = as.matrix(X_subj)          
)

# Compile stan model
mod <- compile_meta_analysis()

init_fun <- function() {
  list(
    alpha = rep(0, p),
    B     = matrix(0, p, q),
    Ltau  = diag(p),
    tau   = rep(0.1, p),          
    Lcorr = diag(p)              
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

# B: matrix p x q
B_vals <- as.numeric(draw_i[grep("^B\\[", param_names)])
B_mat <- matrix(B_vals, nrow = p, ncol = q, byrow = TRUE)

# tau: vector of length p
tau_vals <- as.numeric(draw_i[grep("^tau\\[", param_names)])

# Lcorr: matrix p x p
Lcorr_vals <- as.numeric(draw_i[grep("^Lcorr\\[", param_names)])
Lcorr_mat <- matrix(Lcorr_vals, nrow = p, ncol = p, byrow = TRUE)

# Create a single chain init list
single_init <- list(
  alpha = alpha_vals,
  B     = B_mat,
  tau   = tau_vals,
  Lcorr = Lcorr_mat
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
save(draws, file = paste0("Results/",region_condition,"_",type_condition,"_run",run_id,".RData"))


