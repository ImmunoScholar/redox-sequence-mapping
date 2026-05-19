
# 04_null_model_analysis.R
# Purpose:
# Compare observed target-residue density against randomized sequence controls.
#
# Biological interpretation:
# This tests whether local clustering of Cys/Lys/His exceeds what would be expected
# from the same amino acid composition randomly rearranged.
#
# Limitation:
# This is a sequence-level null model only. It does not model structure,
# solvent accessibility, pKa, or true biochemical reactivity.

library(Biostrings)
library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fasta_file <- "data/raw/example_proteins.fasta"
protein_sequences <- readAAStringSet(fasta_file)

target_residues <- c("C", "K", "H")
window_size <- 10
n_permutations <- 100
set.seed(123)

calculate_max_density <- function(sequence_chars, window_size, residues) {

  sequence_length <- length(sequence_chars)

  if (sequence_length < window_size) {
    return(NA_real_)
  }

  densities <- numeric(sequence_length - window_size + 1)

  for (start_position in seq_along(densities)) {
    end_position <- start_position + window_size - 1
    window_chars <- sequence_chars[start_position:end_position]
    densities[start_position] <- sum(window_chars %in% residues) / window_size
  }

  max(densities)
}

null_results <- list()

for (i in seq_along(protein_sequences)) {

  protein_name <- names(protein_sequences)[i]
  sequence_chars <- strsplit(as.character(protein_sequences[[i]]), "")[[1]]

  observed_max_density <- calculate_max_density(
    sequence_chars,
    window_size,
    target_residues
  )

  randomized_max_densities <- replicate(
    n_permutations,
    calculate_max_density(
      sample(sequence_chars),
      window_size,
      target_residues
    )
  )

  empirical_p_value <- mean(randomized_max_densities >= observed_max_density)

  null_results[[i]] <- data.frame(
    protein = protein_name,
    observed_max_density = observed_max_density,
    mean_randomized_max_density = mean(randomized_max_densities),
    sd_randomized_max_density = sd(randomized_max_densities),
    empirical_p_value = empirical_p_value,
    n_permutations = n_permutations
  )
}

null_model_results <- bind_rows(null_results)

write.csv(
  null_model_results,
  "results/tables/null_model_results.csv",
  row.names = FALSE
)

message("Null-model analysis completed.")
message("Results saved to results/tables/null_model_results.csv")

