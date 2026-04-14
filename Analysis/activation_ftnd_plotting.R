library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Load files
load("Data/voxel_count_summary_df.RData")
ftnd_df <- read.csv("Data/participants_ftnd.csv", stringsAsFactors = FALSE)

# SubjectID linking column
ftnd_df <- ftnd_df %>%
  rename(subject = SubjectID)   

# Join for plotting
plot_df <- voxel_count_summary_df %>%
  left_join(ftnd_df, by = "subject")

# Reshape for plotting
plot_df_long <- plot_df %>%
  select(subject, run, roi, FTND, n_activated, n_deactivated) %>%
  pivot_longer(
    cols = c(n_activated, n_deactivated),
    names_to = "signal_type",
    values_to = "n_voxels"
  ) %>%
  mutate(
    signal_type = recode(
      signal_type,
      n_activated = "Activated",
      n_deactivated = "Deactivated"
    ),
    roi = factor(roi, levels = c("MFG", "Insula", "Precuneus/PCC")),
    signal_type = factor(signal_type, levels = c("Activated", "Deactivated"))
  )


# Plot: 3 x 2 grid
# rows = ROI
# cols = Activated / Deactivated
ggplot(plot_df_long,
       aes(x = FTND,
           y = n_voxels,
           shape = factor(run),
           color = factor(run))) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Run") +
  scale_shape_discrete(name = "Run") +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  facet_grid(roi ~ signal_type) +
  labs(
    x = "FTND",
    y = "Number of voxels",
    title = "ROI Task Related Voxel Counts by FTND"
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold")
  )

