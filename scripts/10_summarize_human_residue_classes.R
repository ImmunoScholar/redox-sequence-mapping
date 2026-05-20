# 10_summarize_human_residue_classes.R
# Purpose:
# Summarise Cys/Lys/His residue-class counts and frequencies across the UniProt-derived human protein set.
#
# Scientific scope:
# This script summarises sequence-level residue-class abundance. It does not infer oxidation, electrophile adduction, biochemical reactivity, structural exposure, or functional perturbation.

library(readr)
library(dplyr)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

human_mapping <- read_csv(
  "results/tables/human_residue_mapping_results.csv",
  show_col_types = FALSE
)

residue_rationale <- read_csv(
  "metadata/residue_class_rationale.csv",
  show_col_types = FALSE
)

summary_table <- human_mapping %>%
  group_by(residue) %>%
  summarise(
    total_count = sum(count),
    mean_frequency = mean(frequency),
    median_frequency = median(frequency),
    min_frequency = min(frequency),
    max_frequency = max(frequency),
    proteins_analysed = n_distinct(protein),
    .groups = "drop"
  ) %>%
  left_join(residue_rationale, by = "residue") %>%
  select(
    residue,
    residue_name,
    chemical_group,
    total_count,
    mean_frequency,
    median_frequency,
    min_frequency,
    max_frequency,
    proteins_analysed,
    rationale,
    interpretation_boundary
  ) %>%
  arrange(residue)

write_csv(
  summary_table,
  "results/tables/human_residue_class_summary.csv"
)

message("Human residue-class summary written to results/tables/human_residue_class_summary.csv")
