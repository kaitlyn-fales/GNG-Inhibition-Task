# Extract and summarize posterior draws from group level model for MCMC
library(posterior)
library(dplyr)
library(R.matlab)

# Number of parameters in group level model of interest
p <- 6

# Names from original notation
par_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_B[1]","nu_B[2]","nu_C[1]")

# Conditions
conditions <- c("MFG_full","Insula_full","MFG_thresh","Insula_thresh")

# Loop through conditions
for (i in 1:length(conditions)){
  load(paste0("Analysis/",conditions[i],".RData"))
  
  # Extract draws for alpha[1] ... alpha[p]
  alpha_vars <- paste0("alpha[", 1:p, "]")
  alpha_draws <- subset_draws(draws, variable = alpha_vars)
  alpha_summary <- summarize_draws(alpha_draws, mean, sd, ~quantile(.x, probs = c(0.025, 0.975))) 
  
  # Extract and summarize tau draws
  tau_vars <- paste0("tau[", 1:p, "]")
  tau_draws <- subset_draws(draws, variable = tau_vars)
  tau_summary <- summarize_draws(tau_draws, mean, sd, ~quantile(.x, probs = c(0.025, 0.975))) 
  
  # General summary table
  summary_table <- tibble(
    parameter = par_names, # original parameter names
    mean_alpha = alpha_summary$mean,
    sd_alpha = alpha_summary$sd,    
    q2.5_alpha = alpha_summary$`2.5%`,
    q97.5_alpha = alpha_summary$`97.5%`,
    between_subject_sd = tau_summary$mean,
    between_subject_q2.5 = tau_summary$`2.5%`,
    between_subject_q97.5 = tau_summary$`97.5%`
  )
  
  assign(conditions[i], summary_table)
}

# Clear environment except for results
rm(list = setdiff(ls(), c("MFG_full","Insula_full","MFG_thresh","Insula_thresh")))

# Combine results into one list once done and export
CDCM_results <- list(MFG_full = MFG_full,
                     Insula_full = Insula_full,
                     MFG_thresh = MFG_thresh,
                     Insula_thresh = Insula_thresh)
save(CDCM_results, file = "Analysis/Results/CDCM_results.RData")


