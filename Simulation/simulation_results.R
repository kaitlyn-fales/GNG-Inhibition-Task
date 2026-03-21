library(cmdstanr)
library(posterior)  
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)

# Function for summarizing simulation results - NoAgg, OldHyp
summarize_results <- function(dir,transform_idx,truth,par_names,
                              file_pattern = c("Type","NoType"),
                              type = c("full","thresh")){
  results_list <- list()
  
  for (snr in 1:8){
    for (rep in 1:50){
      
      if (file_pattern == "Type"){
        file_path <- paste0(dir,"/", type, "_roi_snr", snr, "_", rep, "_draws.RData")
      } else {
        file_path <- paste0(dir,"/snr", snr, "_", rep, "_draws.RData")
      }
      
      if(!file.exists(file_path)) next
      load(file_path) 
      
      # Transform diagonal of A draws
      draws_df[,transform_idx] <- suppressWarnings(-0.5*exp(draws_df[,transform_idx]))
      
      # Keep only the parameters of interest
      draws_subset <- suppressWarnings(draws_df[, par_names])
      
      # Summarize all parameters at once
      summ <- summarise_draws(draws_subset,
                              mean, median, sd,
                              ~quantile(.x, probs = c(0.025,0.975)))
      
      # Loop over parameters to compute coverage and proportion correct sign
      df <- tibble()
      for (i in seq_along(par_names)){
        draws_param <- draws_subset[[ i ]]
        lower <- summ$`2.5%`[i]
        upper <- summ$`97.5%`[i]
        mean_val <- summ$mean[i]
        median_val <- summ$median[i]
        
        df <- bind_rows(df,
                        tibble(
                          snr = snr,
                          rep = rep,
                          param = par_names[i],
                          hpd_lower = lower,
                          hpd_upper = upper,
                          interval_length = upper - lower,
                          coverage = as.numeric(truth[i] >= lower & truth[i] <= upper),
                          bias_mean = mean_val - truth[i],
                          bias_median = median_val - truth[i],
                          prop_correct_sign = mean(sign(draws_param) == sign(truth[i]))
                        ))
      }
      
      results_list[[length(results_list)+1]] <- df
    }
  }
  # Bind rows of list
  sim_results <- bind_rows(results_list)
  
  # Summarize into tidy df
  summary_results <- sim_results %>%
    group_by(snr, param) %>%
    summarise(
      coverage = mean(coverage),
      mean_interval = mean(interval_length),
      median_interval = median(interval_length),
      mean_bias = mean(bias_mean),
      median_bias = median(bias_median),
      prop_correct_sign = mean(prop_correct_sign),
      .groups = "drop"
    )
  
  # Add real SNR values
  SNR_vals <- seq(0.1, 2, length.out = 8)
  
  # Add a numeric column with actual SNR values
  summary_results <- summary_results %>%
    mutate(snr_val = SNR_vals[snr])
  
  return(summary_results)
}

# Function to summarize pre-estimation metrics from PCA/Mean aggregation
read_metrics <- function(file_pattern, method_label, metrics_dir = "Simulation/Metrics") {
  
  files <- list.files(
    path = metrics_dir,
    pattern = file_pattern,
    full.names = TRUE
  )
  
  df <- map_dfr(files, load_metrics_file) %>%
    mutate(method = method_label)
  
  long_df <- df %>%
    pivot_longer(
      cols = -c(rep, snr_index, SNR, sub, run, method),
      names_to = "metric_name",
      values_to = "value"
    ) %>%
    mutate(
      roi = case_when(
        str_detect(metric_name, "MFG") ~ "MFG",
        str_detect(metric_name, "PCC") ~ "PCC",
        TRUE ~ NA_character_
      ),
      scope = case_when(
        str_detect(metric_name, "full") ~ "Full",
        str_detect(metric_name, "thresh") ~ "Thresholded",
        str_detect(metric_name, "selected") ~ "Thresholded",
        TRUE ~ NA_character_
      ),
      metric = case_when(
        str_detect(metric_name, "^cor_") ~ "Signal recovery",
        str_detect(metric_name, "^snr_eff_") ~ "Effective ROI SNR",
        str_detect(metric_name, "^prop_selected") ~ "Proportion selected",
        TRUE ~ metric_name
      )
    )
  
  summary_df <- long_df %>%
    group_by(method, SNR, metric, scope, roi) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      se_value = sd_value / sqrt(n),
      .groups = "drop"
    )
  
  list(raw = df, long = long_df, summary = summary_df)
}

######### Results for original hypothesis (OldHyp) #########
# Parameter names to extract
param_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_A[4]","nu_B[1]","nu_B[2]","nu_B[3]","nu_C[1]")

# Set of true values from simulation
true_vals <- c(-0.5*exp(0.1),-0.3,-0.2,-0.5*exp(0.05),0.15,-0.1,-0.1,0.9) 

# Summarize simulation results
summary_OldHyp <- summarize_results(dir = "Simulation/Output_OldHyp",
                                    transform_idx = c(4,7),
                                    truth = true_vals,
                                    par_names = param_names,
                                    file_pattern = "NoType")

desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_B[3]", "nu_A[4]", "nu_C[1]"
)

summary_OldHyp <- summary_OldHyp %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(8, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:6], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,17,16,15)

# Coverage
p1 <- summary_OldHyp %>%
  ggplot(aes(x = snr_val, y = coverage, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    name = "Parameter",
    values = param_colors,
    labels = c(
      "MFG (Self)","Inhibition -> MFG (Self)", "MFG -> PCC", "Inhibition -> (MFG -> PCC)",
      "PCC -> MFG","Inhibition -> (PCC -> MFG)","PCC (Self)","GNG Task -> MFG"
    )
  ) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Coverage",
       x = "",
       y = "") +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = point_shapes,
        size = 3,
        linewidth = 1
      )
    ),
    linetype = "none",
    shape = "none"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Interval Length
p2 <- summary_OldHyp %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0, 3.2) +
  labs(title = "Interval Length",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Bias
p3 <- summary_OldHyp %>%
  ggplot(aes(x = snr_val, y = mean_bias, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(-1,1) +
  labs(title = "Bias",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Bias
p4 <- summary_OldHyp %>%
  ggplot(aes(x = snr_val, y = prop_correct_sign, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Proportion of Correct Sign",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Extract legend and plot all 
common_legend <- cowplot::get_legend(
  p1 +
    theme(legend.position = "right") +
    guides(
      color = guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(
          shape = point_shapes,
          size = 3,
          linewidth = 1
        )
      ),
      linetype = "none",
      shape = "none"
    )
)
legend_grob <- arrangeGrob(common_legend)
plots_without_legends <- list(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  p4 + theme(legend.position = "none")
)
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("Original State Space Hypothesis",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)

############################################################

rm(list = setdiff(ls(), c("summarize_results","read_metrics")))

##### Results for revised hyp, no aggregation (NoAgg) ######
# Parameter names to extract
param_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_B[1]","nu_B[2]","nu_C[1]")

# Set of true values from simulation
true_vals <- c(-0.5*exp(0.1),-0.3,-0.5*exp(0.05),0.15,-0.1,0.9)

# Summarize simulation results
summary_NoAgg <- summarize_results(dir = "Simulation/Output_NoAgg",
                                   transform_idx = c(4,6),
                                   truth = true_vals,
                                   par_names = param_names,
                                   file_pattern = "NoType")

# Fixed color palette for parameters
desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_C[1]"
)

summary_NoAgg <- summary_NoAgg %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(6, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:4], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,15)

# Coverage
p1 <- summary_NoAgg %>%
  ggplot(aes(x = snr_val, y = coverage, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    name = "Parameter",
    values = param_colors,
    labels = c(
      "MFG (Self)","Inhibition -> MFG (Self)", "MFG -> PCC", "Inhibition -> (MFG -> PCC)",
      "PCC (Self)","GNG Task -> MFG"
    )
  ) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Coverage",
       x = "",
       y = "") +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = point_shapes,
        size = 3,
        linewidth = 1
      )
    ),
    linetype = "none",
    shape = "none"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Interval Length
p2 <- summary_NoAgg %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0, 3.2) +
  labs(title = "Interval Length",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Bias
p3 <- summary_NoAgg %>%
  ggplot(aes(x = snr_val, y = mean_bias, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(-1,1) +
  labs(title = "Bias",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Bias
p4 <- summary_NoAgg %>%
  ggplot(aes(x = snr_val, y = prop_correct_sign, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Proportion of Correct Sign",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Extract legend and plot all 
common_legend <- cowplot::get_legend(
  p1 +
    theme(legend.position = "right") +
    guides(
      color = guide_legend(
        nrow = 2, byrow = TRUE,
        override.aes = list(
          shape = point_shapes,
          size = 3,
          linewidth = 1
        )
      ),
      linetype = "none",
      shape = "none"
    )
)
legend_grob <- arrangeGrob(common_legend)
plots_without_legends <- list(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  p4 + theme(legend.position = "none")
)
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("Simplified Hypothesis",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)

############################################################

rm(list = setdiff(ls(), "read_metrics"))

############## Results for PCA/Mean aggregation #################
# Analyze pre-estimation metrics
pca_out  <- read_metrics("^pca_metrics_snr[0-9]+_[0-9]+\\.RData$", "PCA")
mean_out <- read_metrics("^mean_metrics_snr[0-9]+_[0-9]+\\.RData$", "Mean")

all_metrics_long <- bind_rows(pca_out$long, mean_out$long)
all_metrics_summary <- bind_rows(pca_out$summary, mean_out$summary)

# Parameter names to extract
param_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_B[1]","nu_B[2]","nu_C[1]")

# Set of true values from simulation
true_vals <- c(-0.5*exp(0.1),-0.3,-0.5*exp(0.05),0.15,-0.1,0.9)

# add plotting and analysis of post-estimation metrics

############################################################


