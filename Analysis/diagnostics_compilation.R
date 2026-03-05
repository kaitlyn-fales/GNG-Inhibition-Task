library(dplyr)
library(posterior)
library(mcmcse)
library(momentLS)

path <- paste0(getwd(),"/Output")

# Vector of subjects to loop through
subjects <- c(110,111,120,121,124,128,134,143,152,160,171,172,173,181,
              184,196,199,214,215,223,227,230,247,252,256,258,265,266,268,
              271,275,276,277)

# Run
run <- c(1:3)

# Type of ROI/VOI approach
type <- c("full","thresh")

# Activating region used
region <- c("MFG","Insula")

# Combinations of all phases and mask types
combos <- expand.grid(run = run, VOI_type = type, region = region)

# Empty list to store results
diagnostic_list <- list()

# For loop to make master list of diagnostics across all combinations and subjects
for (i in 1:length(subjects)){
  
  # Temporary empty list for individual subject results as they are compiled into master list
  sub_results_list <- list()

    for (j in 1:nrow(combos)){
    
      # Load output diagnostic file
      load(paste0(path,"/sub-",subjects[i],"_ses1_run",combos$run[j],"_",combos$region[j],"_",combos$VOI_type[j],"_diagnostics.RData"))
      
      # Put in large df
      all_diagnostics <- do.call(rbind, do.call(c, lapply(all_diagnostics, function(x) if (is.data.frame(x)) list(x) else x)))
      
      # Record number of iterations
      iterations <- nrow(all_diagnostics)/6
      
      # Collect max treedepths
      max_treedepth <- all_diagnostics %>% filter(Parameter == "treedepth__" & Value == 10) %>% nrow()/iterations
      
      # Collect divergences
      divergence <- all_diagnostics %>% filter(Parameter == "divergent__" & Value == 1) %>% nrow()/iterations
      
      # Collect E-BFMI
      energy <- all_diagnostics %>% filter(Parameter == "energy__")
      e_bfmi <- (sum(diff(energy$Value)^2) / (length(energy$Value)-2)) / var(energy$Value)
      
      # Collect average acceptance prob
      accept_stat <- as.numeric(all_diagnostics %>% filter(Parameter == "accept_stat__") %>% summarise(mean(Value)))
      
      # Collect stepsize
      stepsize <- as.numeric(all_diagnostics %>% filter(Parameter == "stepsize__") %>% head(n = 1) %>% select(Value))
      
      # Load in output draws file
      load(paste0(path,"/sub-",subjects[i],"_ses1_run",combos$run[j],"_",combos$region[j],"_",combos$VOI_type[j],"_draws.RData"))
      
      # Get summary of param except lp_
      summary <- summarize_draws(draws_df, mcse_mean, ess_bulk, ess_tail)[-1,-1]
      avg_summary <- colMeans(summary)
      names(avg_summary) <- c("mean_mcse","mean_ess_bulk","mean_ess_tail")
      diag_cutoffs <- c(max_mcse = max(summary$mcse_mean),
                        min_ess_bulk = min(summary$ess_bulk),
                        min_ess_tail = min(summary$ess_tail))
      
      # Get multivariate ESS
      param_draws <- suppressWarnings(as.matrix(draws_df[,2:(ncol(draws_df)-3)]))
      avar <- momentLS::mtvMLSE(param_draws)$cov
      ess_multi <- multiESS(param_draws, covmat = avar)
      
      # Convergence indicators
      convergence <- c(mcse_ok = ifelse(diag_cutoffs[1]<0.01,T,F),
                       ess_bulk_ok = ifelse(diag_cutoffs[2]>100,T,F),
                       ess_tail_ok = ifelse(diag_cutoffs[3]>100,T,F),
                       ess_multi_good = ifelse(ess_multi>8793,T,F), # number for minESS(15, 0.05, 0.05)
                       ess_multi_ok = ifelse(ess_multi>2198,T,F)) #   number for minESS(15, 0.05, 0.1)
      
      # Combine into dataframe
      results <- as.data.frame.list(
        c(
          subject = subjects[i],
          run = combos$run[j],
          region = as.character(combos$region[j]),
          VOI_type = as.character(combos$VOI_type[j]),
          iterations = iterations,
          prop_max_treedepth = max_treedepth,
          prop_divergence = divergence,
          e_bfmi = e_bfmi,
          mean_accept_stat = accept_stat,
          stepsize = stepsize,
          avg_summary,
          diag_cutoffs,
          ess_multi = ess_multi,
          convergence
        ),
        stringsAsFactors = FALSE
      )
      
      # Format dataframe properly
      names(results) <- sub("\\..*$", "", names(results))
      results <- type.convert(results, as.is = TRUE)
      
      sub_results_list[[j]] <- results
    }
  
  # Bind rows of combination results for subject i
  sub_results <- bind_rows(sub_results_list)

  # Add to result diagnostic list
  diagnostic_list[[i]] <- sub_results

}

diagnostics_df <- bind_rows(diagnostic_list)

save(diagnostics_df, file = paste0(path,"/../Analysis/diagnostics_compilation.RData"))

# Summarize
summary(diagnostics_df)

