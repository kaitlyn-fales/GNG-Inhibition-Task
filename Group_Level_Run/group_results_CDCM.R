# Extract and summarize posterior draws from group level model for MCMC
library(posterior)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)

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
    
    load(paste0("Group_Level_Run/Results/", conditions[i], "_run", run_id[r], ".RData"))
    
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

rm(list = setdiff(ls(), c("CDCM_MFG_full","CDCM_MFG_thresh","CDCM_Insula_full","CDCM_Insula_thresh")))

# Plotting
all_models <- list(
  MFG_Full = CDCM_MFG_full,
  MFG_Thresh = CDCM_MFG_thresh,
  Insula_Full = CDCM_Insula_full,
  Insula_Thresh = CDCM_Insula_thresh
)

# Parameters to keep for 3x2 plots
params_keep <- c("nu_A[2]", "nu_B[2]", "nu_C[1]")

# Functions to extract baseline connectivity and FTND
extract_baseline <- function(model_obj, roi_name) {
  map_dfr(c("run1", "run2", "run3"), function(run_name) {
    model_obj[[run_name]]$alpha_tau %>%
      filter(parameter %in% params_keep) %>%
      transmute(
        roi = roi_name,
        run = factor(
          run_name,
          levels = c("run1", "run2", "run3"),
          labels = c("Run 1", "Run 2", "Run 3")
        ),
        parameter,
        mean  = mean_alpha,
        lower = q2.5_alpha,
        upper = q97.5_alpha
      )
  })
}

extract_ftnd <- function(model_obj, roi_name) {
  map_dfr(c("run1", "run2", "run3"), function(run_name) {
    model_obj[[run_name]]$B_subj %>%
      filter(covariate == "ftnd_std", parameter %in% params_keep) %>%
      transmute(
        roi = roi_name,
        run = factor(
          run_name,
          levels = c("run1", "run2", "run3"),
          labels = c("Run 1", "Run 2", "Run 3")
        ),
        parameter,
        mean  = mean,
        lower = q2.5,
        upper = q97.5
      )
  })
}

# Build datasets
baseline_df <- imap_dfr(all_models, extract_baseline)
ftnd_df     <- imap_dfr(all_models, extract_ftnd)

add_plot_metadata <- function(df) {
  df %>%
    mutate(
      region = case_when(
        str_detect(roi, "^MFG")    ~ "ROI = MFG",
        str_detect(roi, "^Insula") ~ "ROI = Insula",
        TRUE ~ NA_character_
      ),
      roi_type = case_when(
        str_detect(roi, "Full")   ~ "Full ROI",
        str_detect(roi, "Thresh") ~ "Thresholded ROI",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      parameter_row = case_when(
        parameter == "nu_A[2]" ~ "ROI \u2192 PCC",
        parameter == "nu_B[2]" ~ "Inhibition \u2192 (ROI \u2192 PCC)",
        parameter == "nu_C[1]" ~ "GNG Task \u2192 ROI",
        TRUE ~ parameter
      ),
      parameter_row = factor(
        parameter_row,
        levels = c(
          "ROI \u2192 PCC",
          "Inhibition \u2192 (ROI \u2192 PCC)",
          "GNG Task \u2192 ROI"
        )
      ),
      region = factor(region, levels = c("ROI = MFG", "ROI = Insula"))
    ) %>%
    mutate(
      parameter_label = case_when(
        parameter == "nu_A[2]" & region == "MFG"    ~ "MFG \u2192 PCC",
        parameter == "nu_A[2]" & region == "Insula" ~ "Insula \u2192 PCC",
        parameter == "nu_B[2]" & region == "MFG"    ~ "Inhibition \u2192 (MFG \u2192 PCC)",
        parameter == "nu_B[2]" & region == "Insula" ~ "Inhibition \u2192 (Insula \u2192 PCC)",
        parameter == "nu_C[1]" & region == "MFG"    ~ "GNG Task \u2192 MFG",
        parameter == "nu_C[1]" & region == "Insula" ~ "GNG Task \u2192 Insula",
        TRUE ~ parameter
      )
    )
}

baseline_df2 <- add_plot_metadata(baseline_df) %>%
  mutate(estimate_type = "Run")
ftnd_df2     <- add_plot_metadata(ftnd_df) %>%
  mutate(estimate_type = "Run")

# Load in run-pooled results and format accordingly
load("Analysis/Results/CDCM_results.RData")
all_pooled_models <- list(
  MFG_Full         = CDCM_results$MFG_full,
  MFG_Threshold    = CDCM_results$MFG_thresh,
  Insula_Full      = CDCM_results$Insula_full,
  Insula_Threshold = CDCM_results$Insula_thresh
)

extract_pooled_baseline <- function(model_obj, roi_name) {
  model_obj$alpha_tau %>%
    filter(parameter %in% params_keep) %>%
    transmute(
      roi = roi_name,
      parameter,
      mean = mean_alpha
    )
}

extract_pooled_ftnd <- function(model_obj, roi_name) {
  model_obj$B_subj %>%
    filter(covariate == "ftnd_std", parameter %in% params_keep) %>%
    transmute(
      roi = roi_name,
      parameter,
      mean = mean
    )
}

pooled_baseline_df <- imap_dfr(all_pooled_models, extract_pooled_baseline) %>%
  add_plot_metadata() %>%
  mutate(estimate_type = "Pooled")

pooled_ftnd_df <- imap_dfr(all_pooled_models, extract_pooled_ftnd) %>%
  add_plot_metadata() %>%
  mutate(estimate_type = "Pooled")

# Plotting function
make_run_plot_roi_compare <- function(df, title_text, ylab_text, pooled_df = NULL) {
  dodge <- position_dodge(width = 0.15)
  
  p <- ggplot(
    df,
    aes(
      x = run,
      y = mean,
      color = roi_type,
      shape = roi_type,
      group = interaction(roi_type, estimate_type)
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.5
    )
  
  if (!is.null(pooled_df)) {
    p <- p +
      geom_hline(
        data = pooled_df,
        aes(
          yintercept = mean,
          color = roi_type,
          linetype = estimate_type
        ),
        linewidth = 0.6,
        inherit.aes = FALSE
      )
  }
  
  p +
    geom_line(
      aes(linetype = estimate_type),
      linewidth = 0.8,
      position = dodge
    ) +
    geom_point(
      size = 3.0,
      position = dodge
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.08,
      linewidth = 0.7,
      position = dodge
    ) +
    facet_grid(
      rows = vars(parameter_row),
      cols = vars(region),
      scales = "fixed"
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(
      "Full ROI" = 16,
      "Thresholded ROI" = 17
    )) +
    scale_linetype_manual(values = c(
      "Run" = "solid",
      "Pooled" = "dotted"
    )) +
    labs(
      title = title_text,
      x = "Run",
      y = ylab_text,
      color = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(face = "bold", size = 11),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey85"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(face = "bold", size = 11)
    )
}

row_levels <- c(
  "MFG \u2192 PCC",
  "Inhibition \u2192 (MFG \u2192 PCC)",
  "GNG Task \u2192 MFG",
  "Insula \u2192 PCC",
  "Inhibition \u2192 (Insula \u2192 PCC)",
  "GNG Task \u2192 Insula"
)

baseline_df2 <- baseline_df2 %>%
  mutate(parameter_label = factor(parameter_label, levels = row_levels))

ftnd_df2 <- ftnd_df2 %>%
  mutate(parameter_label = factor(parameter_label, levels = row_levels))

pooled_baseline_df <- pooled_baseline_df %>%
  mutate(parameter_label = factor(parameter_label, levels = row_levels))

pooled_ftnd_df <- pooled_ftnd_df %>%
  mutate(parameter_label = factor(parameter_label, levels = row_levels))

p_baseline <- make_run_plot_roi_compare(
  baseline_df2,
  title_text = "Baseline Connectivity by Run",
  ylab_text = "Posterior mean (95% HDI)",
  pooled_df = pooled_baseline_df
)

p_ftnd <- make_run_plot_roi_compare(
  ftnd_df2,
  title_text = "FTND Effects by Run",
  ylab_text = "FTND effect (95% HDI)",
  pooled_df = pooled_ftnd_df
)

p_baseline
p_ftnd
