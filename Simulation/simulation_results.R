library(cmdstanr)
library(posterior)  
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(purrr)
library(RColorBrewer)

# Function for summarizing simulation results
summarize_results <- function(dir, transform_idx, truth, par_names,
                              file_pattern = c("Type", "NoType"),
                              type = c("full", "thresh"),
                              dist = c("beta_sym","beta_u")) {
  
  results_list <- list()
  
  for (snr in 1:8) {
    for (rep in 1:50) {
      
      if (file_pattern == "Type") {
        file_path <- paste0(dir, "/", dist, "_", type, "_roi_snr", snr, "_", rep, "_draws.RData")
      } else {
        file_path <- paste0(dir, "/snr", snr, "_", rep, "_draws.RData")
      }
      
      if (!file.exists(file_path)) next
      load(file_path)
      
      # Transform diagonal of A draws
      draws_df[, transform_idx] <- suppressWarnings(-0.5 * exp(draws_df[, transform_idx]))
      
      # Keep only the parameters of interest
      draws_subset <- suppressWarnings(draws_df[, par_names])
      
      # Summarize all parameters at once
      summ <- summarise_draws(
        draws_subset,
        mean, median, sd,
        ~quantile(.x, probs = c(0.025, 0.975))
      )
      
      # Loop over parameters
      df <- tibble()
      for (i in seq_along(par_names)) {
        draws_param <- draws_subset[[i]]
        lower <- summ$`2.5%`[i]
        upper <- summ$`97.5%`[i]
        mean_val <- summ$mean[i]
        median_val <- summ$median[i]
        
        # Posterior probability of being on the correct side of 0
        post_prob_correct_side <- if (truth[i] > 0) {
          mean(draws_param > 0)
        } else if (truth[i] < 0) {
          mean(draws_param < 0)
        } else {
          NA_real_
        }
        
        df <- bind_rows(
          df,
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
            prop_correct_sign = mean(sign(draws_param) == sign(truth[i])),
            post_prob_correct_side = post_prob_correct_side
          )
        )
      }
      
      results_list[[length(results_list) + 1]] <- df
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
      median_bias = mean(bias_median),
      prop_correct_sign = mean(prop_correct_sign),
      mean_post_prob_correct_side = mean(post_prob_correct_side, na.rm = TRUE),
      median_post_prob_correct_side = median(post_prob_correct_side, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add real SNR values
  SNR_vals <- seq(0.1, 2, length.out = 8)
  
  summary_results <- summary_results %>%
    mutate(snr_val = SNR_vals[snr])
  
  return(summary_results)
}

# Wrapper function to load metrics file
load_metrics_file <- function(file_path) {
  e <- new.env()
  load(file_path, envir = e)
  
  obj_name <- ls(e)
  
  if (length(obj_name) != 1) {
    stop(paste("File", file_path, "contains multiple objects."))
  }
  
  df <- e[[obj_name]]
  df <- as.data.frame(df)
  
  return(df)
}

# Function to summarize pre-estimation metrics from PCA/Mean aggregation
read_metrics <- function(file_pattern, method_label, metrics_dir = "Simulation/Metrics") {
  
  files <- list.files(
    path = metrics_dir,
    pattern = file_pattern,
    full.names = TRUE
  )
  
  df <- map_dfr(files, function(f) {
    out <- load_metrics_file(f)
    
    # Extract full distribution label from filename
    out$dist <- sub("_(pca|mean)_metrics_snr[0-9]+_[0-9]+\\.RData$", "", basename(f))
    out
  }) %>%
    mutate(method = method_label)
  
  long_df <- df %>%
    pivot_longer(
      cols = -c(rep, snr_index, SNR, sub, run, method, dist),
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
    group_by(dist, method, SNR, metric, scope, roi) %>%
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"))

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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Bias
p4 <- summary_OldHyp %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Posterior Probability of Correct Sign",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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

rm(list = setdiff(ls(), c("summarize_results","read_metrics","load_metrics_file")))

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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"))

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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Bias
p4 <- summary_NoAgg %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Posterior Probability of Correct Sign",
       x = "",
       y = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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

rm(list = setdiff(ls(), c("summarize_results","read_metrics","load_metrics_file")))

############## Results for PCA/Mean aggregation #################
# Analyze pre-estimation metrics
pca_out  <- read_metrics("^beta_.*_pca_metrics_snr[0-9]+_[0-9]+\\.RData$", "PCA")
mean_out <- read_metrics("^beta_.*_mean_metrics_snr[0-9]+_[0-9]+\\.RData$", "Mean")

all_metrics_long <- bind_rows(pca_out$long, mean_out$long)
all_metrics_summary <- bind_rows(pca_out$summary, mean_out$summary)

dist_labels <- c(
  beta_sym = "Beta(3,3)",
  beta_u   = "Beta(0.5,0.5)"
)

paired_cols <- brewer.pal(10, "Paired")[c(4, 10)]
names(paired_cols) <- c("PCA", "Mean")

roi_cols <- brewer.pal(6, "Paired")[c(2, 6)]
names(roi_cols) <- c("MFG", "PCC")

make_main_plot <- function(dat, dist_label) {
  dat %>%
    filter(
      metric %in% c("Signal recovery", "Effective ROI SNR"),
      dist == dist_label
    ) %>%
    mutate(
      metric = factor(metric, levels = c("Signal recovery", "Effective ROI SNR")),
      scope  = factor(scope, levels = c("Full", "Thresholded")),
      method = factor(method, levels = c("PCA", "Mean")),
      roi    = factor(roi, levels = c("MFG", "PCC"))
    ) %>%
    ggplot(aes(
      x = SNR,
      y = mean_value,
      color = method,
      linetype = scope,
      group = interaction(method, scope)
    )) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    geom_errorbar(
      aes(ymin = mean_value - se_value,
          ymax = mean_value + se_value),
      width = 0.05,
      linewidth = 0.5,
      linetype = "solid"
    ) +
    facet_grid(metric ~ roi, scales = "free_y") +
    scale_color_manual(values = paired_cols) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(
      title = paste("Pre-Estimation Performance by ROI Approach and Summary Method -",
                    dist_labels[dist_label]),
      x = NULL,
      y = NULL,
      color = "Summary",
      linetype = "ROI approach"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      axis.text = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
      linetype = guide_legend(nrow = 1, byrow = TRUE, order = 2)
    )
}

make_prop_plot <- function(dat, dist_label) {
  dat %>%
    filter(
      metric == "Proportion selected",
      scope == "Thresholded",
      method == "PCA",
      dist == dist_label
    ) %>%
    mutate(
      roi = factor(roi, levels = c("MFG", "PCC"))
    ) %>%
    ggplot(aes(
      x = SNR,
      y = mean_value,
      color = roi,
      group = roi
    )) +
    geom_line(linewidth = 0.7, linetype = "dashed") +
    geom_point(size = 1.5) +
    geom_errorbar(
      aes(ymin = mean_value - se_value,
          ymax = mean_value + se_value),
      width = 0.05,
      linewidth = 0.5,
      linetype = "solid"
    ) +
    scale_color_manual(values = roi_cols) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = "Voxelwise Signal-to-Noise Ratio (SNR)",
      y = NULL,
      title = paste("Proportion of Voxels Selected -", dist_labels[dist_label]),
      color = "ROI"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "white"),
      legend.title = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    guides(
      color = guide_legend(nrow = 1)
    )
}

# Symmetric scaling plot
plot_main_sym <- make_main_plot(all_metrics_summary, "beta_sym")
plot_prop_sym <- make_prop_plot(all_metrics_summary, "beta_sym")

main_legend_sym <- get_legend(
  plot_main_sym + theme(legend.position = "right")
)

prop_legend_sym <- get_legend(
  plot_prop_sym + theme(legend.position = "right")
)

plot_main_sym_nolegend <- plot_main_sym + theme(legend.position = "none")
plot_prop_sym_nolegend <- plot_prop_sym + theme(legend.position = "none")

combined_legend_sym <- plot_grid(
  main_legend_sym,
  prop_legend_sym,
  nrow = 1,
  rel_widths = c(1.4, 1)
)

final_pre_plot_sym <- plot_grid(
  plot_main_sym_nolegend,
  plot_prop_sym_nolegend,
  combined_legend_sym,
  ncol = 1,
  rel_heights = c(2.8, 1.4, 0.4)
)

final_pre_plot_sym


# U-shaped scaling plot
plot_main_u <- make_main_plot(all_metrics_summary, "beta_u")
plot_prop_u <- make_prop_plot(all_metrics_summary, "beta_u")

main_legend_u <- get_legend(
  plot_main_u + theme(legend.position = "right")
)

prop_legend_u <- get_legend(
  plot_prop_u + theme(legend.position = "right")
)

plot_main_u_nolegend <- plot_main_u + theme(legend.position = "none")
plot_prop_u_nolegend <- plot_prop_u + theme(legend.position = "none")

combined_legend_u <- plot_grid(
  main_legend_u,
  prop_legend_u,
  nrow = 1,
  rel_widths = c(1.4, 1)
)

final_pre_plot_u <- plot_grid(
  plot_main_u_nolegend,
  plot_prop_u_nolegend,
  combined_legend_u,
  ncol = 1,
  rel_heights = c(2.8, 1.4, 0.4)
)

final_pre_plot_u

# Clear environment
rm(list = setdiff(ls(), c("summarize_results","read_metrics","load_metrics_file")))

#### Beta(3,3) symmetric

# Parameter names to extract
param_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_B[1]","nu_B[2]","nu_C[1]")

# Set of true values from simulation
true_vals <- c(-0.5*exp(0.1),-0.3,-0.5*exp(0.05),0.15,-0.1,0.9)

# Summarize simulation results
summary_PCA <- summarize_results(dir = "Simulation/Output_PCA/beta_sym",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "full",
                                 dist = "beta_sym")

# Fixed color palette for parameters
desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_C[1]"
)

summary_PCA_full <- summary_PCA %>%
  mutate(param = factor(param, levels = desired_order))

# Summarize simulation results
summary_PCA <- summarize_results(dir = "Simulation/Output_PCA/beta_sym",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "thresh",
                                 dist = "beta_sym")

summary_PCA_thresh <- summary_PCA %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(6, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:4], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,15)

# Interval Length
p1 <- summary_PCA_full %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
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
  ylim(0,1.2) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold")
        )

# Post. prob of correct sign
p2 <- summary_PCA_full %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Interval Length
p3 <- summary_PCA_thresh %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1.2) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Post. prob of correct sign
p4 <- summary_PCA_thresh %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  p4 + theme(legend.position = "none"))
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("PCA Region Aggregation: Beta(3,3)",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Voxelwise Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)


# Summarize simulation results
summary_mean <- summarize_results(dir = "Simulation/Output_Mean/beta_sym",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "full",
                                 dist = "beta_sym")

# Fixed color palette for parameters
desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_C[1]"
)

summary_mean_full <- summary_mean %>%
  mutate(param = factor(param, levels = desired_order))

# Summarize simulation results
summary_mean <- summarize_results(dir = "Simulation/Output_Mean/beta_sym",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "thresh",
                                 dist = "beta_sym")

summary_mean_thresh <- summary_mean %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(6, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:4], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,15)

# Interval Length
p1 <- summary_mean_full %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
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
  ylim(0,1.2) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold")
  )

# Post. prob of correct sign
p2 <- summary_mean_full %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Interval Length
p3 <- summary_mean_thresh %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1.2) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Post. prob of correct sign
p4 <- summary_mean_thresh %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  p4 + theme(legend.position = "none"))
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("Mean Region Aggregation: Beta(3,3)",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Voxelwise Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)

# Clear environment
rm(list = setdiff(ls(), c("summarize_results","read_metrics","load_metrics_file")))

#### Beta(0.5,0.5) U-shape

# Parameter names to extract
param_names <- c("nu_A[1]","nu_A[2]","nu_A[3]","nu_B[1]","nu_B[2]","nu_C[1]")

# Set of true values from simulation
true_vals <- c(-0.5*exp(0.1),-0.3,-0.5*exp(0.05),0.15,-0.1,0.9)

# Summarize simulation results
summary_PCA <- summarize_results(dir = "Simulation/Output_PCA/beta_u",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "full",
                                 dist = "beta_u")

# Fixed color palette for parameters
desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_C[1]"
)

summary_PCA_full <- summary_PCA %>%
  mutate(param = factor(param, levels = desired_order))

# Summarize simulation results
summary_PCA <- summarize_results(dir = "Simulation/Output_PCA/beta_u",
                                 transform_idx = c(4,6),
                                 truth = true_vals,
                                 par_names = param_names,
                                 file_pattern = "Type",
                                 type = "thresh",
                                 dist = "beta_u")

summary_PCA_thresh <- summary_PCA %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(6, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:4], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,15)

# Interval Length
p1 <- summary_PCA_full %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
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
  ylim(0,1.2) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold")
  )

# Post. prob of correct sign
p2 <- summary_PCA_full %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Interval Length
p3 <- summary_PCA_thresh %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1.2) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Post. prob of correct sign
p4 <- summary_PCA_thresh %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  p4 + theme(legend.position = "none"))
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("PCA Region Aggregation: Beta(0.5,0.5)",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Voxelwise Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)


# Summarize simulation results
summary_mean <- summarize_results(dir = "Simulation/Output_Mean/beta_u",
                                  transform_idx = c(4,6),
                                  truth = true_vals,
                                  par_names = param_names,
                                  file_pattern = "Type",
                                  type = "full",
                                  dist = "beta_u")

# Fixed color palette for parameters
desired_order <- c(
  "nu_A[1]", "nu_B[1]", "nu_A[2]", "nu_B[2]",
  "nu_A[3]", "nu_C[1]"
)

summary_mean_full <- summary_mean %>%
  mutate(param = factor(param, levels = desired_order))

# Summarize simulation results
summary_mean <- summarize_results(dir = "Simulation/Output_Mean/beta_u",
                                  transform_idx = c(4,6),
                                  truth = true_vals,
                                  par_names = param_names,
                                  file_pattern = "Type",
                                  type = "thresh",
                                  dist = "beta_u")

summary_mean_thresh <- summary_mean %>%
  mutate(param = factor(param, levels = desired_order))

# Fixed color palette for parameters
param_colors <- RColorBrewer::brewer.pal(6, "Paired")  # your first 6 colors
# pick extra colors from larger palette
extra_colors <- RColorBrewer::brewer.pal(12, "Paired")[c(8,10)]
param_colors <- c(param_colors[1:4], extra_colors)
names(param_colors) <- desired_order

point_shapes <- c(16,17,16,17,16,15)

# Interval Length
p1 <- summary_mean_full %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
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
  ylim(0,1.2) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold")
  )

# Post. prob of correct sign
p2 <- summary_mean_full %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "Full ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Interval Length
p3 <- summary_mean_thresh %>%
  ggplot(aes(x = snr_val, y = mean_interval, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(
    values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1.2) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Interval Length") +
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
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Post. prob of correct sign
p4 <- summary_mean_thresh %>%
  ggplot(aes(x = snr_val, y = mean_post_prob_correct_side, color = param, shape = param)) +
  geom_line(linewidth = 1) +
  geom_point(size=3) +
  scale_color_manual(values = param_colors) +
  scale_shape_manual(values = point_shapes) +
  ylim(0,1) +
  labs(title = "GLM Thresholded ROI",
       x = NULL,
       y = "Posterior Prob. of Correct Sign") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

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
  p4 + theme(legend.position = "none"))
plot_grid <- arrangeGrob(
  grobs = plots_without_legends,
  ncol = 2, nrow = 2,
  top = textGrob("Mean Region Aggregation: Beta(0.5,0.5)",
                 gp = gpar(fontface = "bold", fontsize = 16)),
  bottom = textGrob("Voxelwise Signal-to-Noise Ratio (SNR)", gp = gpar(fontface = "bold", fontsize = 12))
)
grid.arrange(
  plot_grid,
  legend_grob,
  ncol = 1,
  heights = c(10, 1.2)
)

############################################################


