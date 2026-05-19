# 08_summarize_human_results.R
# Purpose:
# Summarise UniProt-derived human protein sequence analysis results.
#
# Scientific scope:
# This summary ranks sequence-level residue clustering outputs. It does not
# infer biochemical modification, redox reactivity, structural exposure, or
# functional perturbation.

library(readr)
library(dplyr)
library(stringr)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

null_results <- read_csv(
  "results/tables/human_null_model_results.csv",
  show_col_types = FALSE
)

extract_uniprot_accession <- function(header) {
  str_match(header, "^sp\\|([^|]+)\\|")[, 2]
}

extract_uniprot_name <- function(header) {
  str_match(header, "^sp\\|[^|]+\\|([^ ]+)")[, 2]
}

extract_gene <- function(header) {
  str_match(header, "GN=([^ ]+)")[, 2]
}

summary_table <- null_results %>%
  mutate(
    accession = extract_uniprot_accession(protein),
    uniprot_name = extract_uniprot_name(protein),
    gene = extract_gene(protein),
    density_difference = observed_max_density - mean_randomized_max_density,
    evidence_interpretation = case_when(
      empirical_p_value < 0.05 & density_difference > 0 ~
        "Observed clustering exceeds randomized expectation at nominal p < 0.05",
      density_difference > 0 ~
        "Observed clustering is higher than randomized expectation, but not statistically strong",
      TRUE ~
        "Observed clustering does not exceed randomized expectation"
    )
  ) %>%
  select(
    accession,
    uniprot_name,
    gene,
    observed_max_density,
    mean_randomized_max_density,
    density_difference,
    empirical_p_value,
    n_permutations,
    evidence_interpretation
  ) %>%
  arrange(empirical_p_value, desc(density_difference))

write_csv(
  summary_table,
  "results/tables/human_sequence_clustering_summary.csv"
)

message("Human sequence clustering summary written to results/tables/human_sequence_clustering_summary.csv")