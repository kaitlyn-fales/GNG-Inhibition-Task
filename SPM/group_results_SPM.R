# Extract and summarize posterior draws from group level model for MCMC
library(posterior)
library(dplyr)
library(tibble)

# Number of parameters in group level model of interest
p <- 6

# Names from original notation
par_names <- c("nu_A[1]", "nu_A[2]", "nu_A[3]", "nu_B[1]", "nu_B[2]", "nu_C[1]")

# Covariate names in the exact order used to build X_subj and X_run
subj_cov_names <- c("age_std", "sex_code", "ftnd_std")
run_cov_names  <- c("accuracy_std")

# Conditions
conditions <- c("MFG_full", "Insula_full", "MFG_thresh", "Insula_thresh")

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

# Loop through conditions
for (i in seq_along(conditions)) {
  
  load(paste0("SPM/", conditions[i], ".RData"))
  
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
  # subject-level covariates: B_subj[p, q_subj]
  # order:
  # 1 = age_std
  # 2 = sex_code
  # 3 = ftnd_std
  # -----------------------------
  B_subj_vars <- c()
  for (j in 1:p) {
    for (k in 1:length(subj_cov_names)) {
      B_subj_vars <- c(B_subj_vars, paste0("B_subj[", j, ",", k, "]"))
    }
  }
  
  B_subj_draws <- subset_draws(draws, variable = B_subj_vars)
  B_subj_summary <- summ_fun(B_subj_draws)
  
  B_subj_index <- expand.grid(
    covariate = subj_cov_names,
    parameter = par_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::select(parameter, covariate)
  
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
  # run-level covariates: B_run[p, q_run]
  # order:
  # 1 = accuracy_std
  # -----------------------------
  B_run_vars <- c()
  for (j in 1:p) {
    for (k in 1:length(run_cov_names)) {
      B_run_vars <- c(B_run_vars, paste0("B_run[", j, ",", k, "]"))
    }
  }
  
  B_run_draws <- subset_draws(draws, variable = B_run_vars)
  B_run_summary <- summ_fun(B_run_draws)
  
  B_run_index <- expand.grid(
    covariate = run_cov_names,
    parameter = par_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::select(parameter, covariate)
  
  B_run_table <- tibble(
    parameter = B_run_index$parameter,
    covariate = B_run_index$covariate,
    mean = B_run_summary$mean,
    sd = B_run_summary$sd,
    q2.5 = B_run_summary$q2.5,
    q97.5 = B_run_summary$q97.5,
    prob_gt0 = B_run_summary$prob_gt0,
    prob_lt0 = B_run_summary$prob_lt0,
    prob_sign = B_run_summary$prob_sign
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
    B_subj = B_subj_table,
    B_run = B_run_table
  )
  
  assign(conditions[i], condition_results)
}

# Clear environment except for results
rm(list = setdiff(ls(), c("MFG_full", "Insula_full", "MFG_thresh", "Insula_thresh")))

# Combine results into one list once done and export
SPM_results <- list(
  MFG_full = MFG_full,
  Insula_full = Insula_full,
  MFG_thresh = MFG_thresh,
  Insula_thresh = Insula_thresh
)

save(SPM_results, file = "SPM/SPM_results.RData")
