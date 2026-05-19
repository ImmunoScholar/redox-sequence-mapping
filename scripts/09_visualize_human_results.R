# 09_visualize_human_results.R
# Purpose:
# Generate visual summaries for the UniProt-derived human protein sequence analysis.
#
# Scientific scope:
# These figures visualise sequence-level residue frequencies and null-model
# comparisons. They do not demonstrate biochemical reactivity, oxidative
# modification, electrophile adduction, structural exposure, or functional
# perturbation.

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

residue_mapping <- read_csv(
  "results/tables/human_residue_mapping_results.csv",
  show_col_types = FALSE
)

null_results <- read_csv(
  "results/tables/human_null_model_results.csv",
  show_col_types = FALSE
)

extract_gene <- function(header) {
  str_match(header, "GN=([^ ]+)")[, 2]
}

residue_mapping <- residue_mapping %>%
  mutate(gene = extract_gene(protein))

null_plot_data <- null_results %>%
  mutate(
    gene = extract_gene(protein),
    density_difference = observed_max_density - mean_randomized_max_density
  )

residue_frequency_plot <- ggplot(
  residue_mapping,
  aes(x = gene, y = frequency, fill = residue)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Cys/Lys/His Frequencies in UniProt-Derived Human Proteins",
    x = "Gene",
    y = "Residue frequency",
    fill = "Residue"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "results/figures/human_residue_frequency_plot.png",
  plot = residue_frequency_plot,
  width = 9,
  height = 5,
  dpi = 300
)

null_model_plot <- ggplot(
  null_plot_data,
  aes(x = reorder(gene, density_difference), y = density_difference)
) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Observed vs Randomized Maximum Residue Density in Human Proteins",
    subtitle = "Positive values indicate observed maximum density above the composition-preserving randomized mean",
    x = "Gene",
    y = "Observed max density - randomized mean max density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "results/figures/human_null_model_comparison_plot.png",
  plot = null_model_plot,
  width = 8,
  height = 5,
  dpi = 300
)

message("Human residue frequency plot written to results/figures/human_residue_frequency_plot.png")
message("Human null-model comparison plot written to results/figures/human_null_model_comparison_plot.png")