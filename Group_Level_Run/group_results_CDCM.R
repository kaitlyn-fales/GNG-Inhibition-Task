# Extract and summarize posterior draws from group level model for MCMC
library(posterior)
library(dplyr)
library(tibble)

# Number of parameters in group level model of interest
p <- 6

# Names from original notation
par_names <- c("nu_A[1]", "nu_A[2]", "nu_A[3]", "nu_B[1]", "nu_B[2]", "nu_C[1]")

# Covariate names in the exact order used to build X_subj and X_run
subj_cov_names <- c("age_std", "sex_code", "ftnd_std","accuracy_std")

# Conditions and runs
conditions <- c("MFG_full", "Insula_full", "MFG_thresh", "Insula_thresh")
run_id <- c(1:3)

# Robust helper function to summarize draws
summ_fun <- function(draw_obj) {
  mat <- as_draws_matrix(draw_obj)
  
  out <- tibble(
    variable = colnames(mat),
    mean = apply(mat, 2, mean),
    sd = apply(mat, 2, sd),
    q2.5 = apply(mat, 2, quantile, probs = 0.025),
    q97.5 = apply(mat, 2, quantile, probs = 0.975),
    prob_gt0 = apply(mat, 2, function(x) mean(x > 0)),
    prob_lt0 = apply(mat, 2, function(x) mean(x < 0)),
    prob_sign = apply(mat, 2, function(x) max(mean(x > 0), mean(x < 0)))
  )
  
  return(out)
}

# Loop through runs
for (r in 1:length(run_id)){
  # Loop through conditions
  for (i in 1:length(conditions)) {
    
    load(paste0("Group_Level_Run/", conditions[i], "_run", run_id[r], ".RData"))
    
    # -----------------------------
    # alpha
    # -----------------------------
    alpha_vars <- paste0("alpha[", 1:p, "]")
    alpha_draws <- subset_draws(draws, variable = alpha_vars)
    alpha_summary <- summ_fun(alpha_draws)
    
    # -----------------------------
    # tau
    # -----------------------------
    tau_vars <- paste0("tau[", 1:p, "]")
    tau_draws <- subset_draws(draws, variable = tau_vars)
    tau_summary <- summ_fun(tau_draws)
    
    # -----------------------------
    # subject-level covariates: B[p, q]
    # order:
    # 1 = age_std
    # 2 = sex_code
    # 3 = ftnd_std
    # 4 = accuracy_std
    # -----------------------------
    B_subj_vars <- c()
    for (j in 1:p) {
      for (k in 1:length(subj_cov_names)) {
        B_subj_vars <- c(B_subj_vars, paste0("B[", j, ",", k, "]"))
      }
    }
    
    B_subj_draws <- subset_draws(draws, variable = B_subj_vars)
    B_subj_summary <- summ_fun(B_subj_draws)
    
    B_subj_index <- expand.grid(
      parameter = par_names,
      covariate = subj_cov_names,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    B_subj_table <- tibble(
      parameter = B_subj_index$parameter,
      covariate = B_subj_index$covariate,
      mean = B_subj_summary$mean,
      sd = B_subj_summary$sd,
      q2.5 = B_subj_summary$q2.5,
      q97.5 = B_subj_summary$q97.5,
      prob_gt0 = B_subj_summary$prob_gt0,
      prob_lt0 = B_subj_summary$prob_lt0,
      prob_sign = B_subj_summary$prob_sign
    )
    
    # -----------------------------
    # alpha + tau summary table
    # -----------------------------
    alpha_tau_table <- tibble(
      parameter = par_names,
      mean_alpha = alpha_summary$mean,
      sd_alpha = alpha_summary$sd,
      q2.5_alpha = alpha_summary$q2.5,
      q97.5_alpha = alpha_summary$q97.5,
      alpha_prob_gt0 = alpha_summary$prob_gt0,
      alpha_prob_lt0 = alpha_summary$prob_lt0,
      alpha_prob_sign = alpha_summary$prob_sign,
      between_subject_sd = tau_summary$mean,
      between_subject_sd_sd = tau_summary$sd,
      between_subject_q2.5 = tau_summary$q2.5,
      between_subject_q97.5 = tau_summary$q97.5
    )
    
    # -----------------------------
    # Save all summaries for this condition
    # -----------------------------
    condition_results <- list(
      alpha_tau = alpha_tau_table,
      B_subj = B_subj_table
    )
    
    assign(paste0(conditions[i],"_run",run_id[r]), condition_results)
  }
  
}

# Combine results into one list per condition across runs
CDCM_MFG_full <- list(
  run1 = MFG_full_run1,
  run2 = MFG_full_run2,
  run3 = MFG_full_run3
)

CDCM_MFG_thresh <- list(
  run1 = MFG_thresh_run1,
  run2 = MFG_thresh_run2,
  run3 = MFG_thresh_run3
)

CDCM_Insula_full <- list(
  run1 = Insula_full_run1,
  run2 = Insula_full_run2,
  run3 = Insula_full_run3
)

CDCM_Insula_thresh <- list(
  run1 = Insula_thresh_run1,
  run2 = Insula_thresh_run2,
  run3 = Insula_thresh_run3
)

