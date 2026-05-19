
# 05_visualization.R
# Purpose:
# Generate visualization outputs for sequence-level oxidative stress residue analysis.
#
# IMPORTANT:
# These visualizations reflect sequence composition and clustering only.
# They do not imply biochemical modification or structural accessibility.

library(tidyverse)
library(patchwork)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load results
residue_data <- read.csv(
  "results/tables/residue_mapping_results.csv"
)

window_data <- read.csv(
  "results/tables/sliding_window_density.csv"
)

null_data <- read.csv(
  "results/tables/null_model_results.csv"
)

# -----------------------------
# Plot 1: Residue Frequencies
# -----------------------------

frequency_plot <- ggplot(
  residue_data,
  aes(
    x = residue,
    y = frequency,
    fill = residue
  )
) +
  geom_bar(
    stat = "identity",
    position = "dodge"
  ) +
  facet_wrap(~ protein) +
  labs(
    title = "Frequency of Oxidative Stress-Associated Residues",
    x = "Residue",
    y = "Frequency"
  ) +
  theme_bw()

# -----------------------------
# Plot 2: Sliding Window Density
# -----------------------------

window_plot <- ggplot(
  window_data,
  aes(
    x = window_start,
    y = target_residue_density,
    color = protein
  )
) +
  geom_line(linewidth = 1) +
  labs(
    title = "Sliding-Window Residue Density",
    x = "Sequence Position",
    y = "Residue Density"
  ) +
  theme_bw()

# -----------------------------
# Plot 3: Null Model Comparison
# -----------------------------

null_plot <- ggplot(
  null_data,
  aes(
    x = protein,
    y = observed_max_density
  )
) +
  geom_col(fill = "steelblue") +
  geom_point(
    aes(y = mean_randomized_max_density),
    color = "red",
    size = 3
  ) +
  labs(
    title = "Observed vs Randomized Maximum Density",
    x = "Protein",
    y = "Maximum Density"
  ) +
  theme_bw()

# Save figures

ggsave(
  "results/figures/residue_frequency_plot.png",
  frequency_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "results/figures/sliding_window_density_plot.png",
  window_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "results/figures/null_model_comparison_plot.png",
  null_plot,
  width = 8,
  height = 5,
  dpi = 300
)

message("Visualization pipeline completed.")
message("Figures saved in results/figures/")

