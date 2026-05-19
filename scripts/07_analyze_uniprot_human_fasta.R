# 07_analyze_uniprot_human_fasta.R
# Purpose:
# Run sequence-level residue mapping, sliding-window density analysis, and
# composition-preserving null-model analysis on UniProt-derived human proteins.
#
# Scientific scope:
# This script performs primary-sequence analysis only. It does not infer
# biochemical reactivity, oxidative modification, electrophile adduction,
# structural exposure, or functional perturbation.

library(Biostrings)
library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fasta_file <- "data/raw/human_redox_relevant_proteins_uniprot.fasta"

if (!file.exists(fasta_file)) {
  stop("Human UniProt FASTA file not found. Run scripts/06_download_uniprot_human_fasta.R first.")
}

protein_sequences <- readAAStringSet(fasta_file)

target_residues <- c("C", "K", "H")
window_size <- 10
n_permutations <- 1000
set.seed(123)

count_residues <- function(sequence, residues) {
  sequence_chars <- strsplit(as.character(sequence), "")[[1]]
  total_length <- length(sequence_chars)

  residue_counts <- sapply(residues, function(residue) {
    sum(sequence_chars == residue)
  })

  residue_frequencies <- residue_counts / total_length

  data.frame(
    residue = residues,
    count = as.integer(residue_counts),
    frequency = as.numeric(residue_frequencies)
  )
}

calculate_window_density <- function(sequence, protein_name, window_size, residues) {
  sequence_chars <- strsplit(as.character(sequence), "")[[1]]
  sequence_length <- length(sequence_chars)

  if (sequence_length < window_size) {
    return(NULL)
  }

  window_results <- vector("list", sequence_length - window_size + 1)

  for (start_position in seq_along(window_results)) {
    end_position <- start_position + window_size - 1
    window_chars <- sequence_chars[start_position:end_position]

    residue_count <- sum(window_chars %in% residues)
    residue_density <- residue_count / window_size

    window_results[[start_position]] <- data.frame(
      protein = protein_name,
      window_start = start_position,
      window_end = end_position,
      window_size = window_size,
      target_residue_count = residue_count,
      target_residue_density = residue_density
    )
  }

  bind_rows(window_results)
}

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

residue_results <- list()
window_results <- list()
null_results <- list()

for (i in seq_along(protein_sequences)) {
  protein_name <- names(protein_sequences)[i]
  sequence <- protein_sequences[[i]]
  sequence_chars <- strsplit(as.character(sequence), "")[[1]]

  residue_data <- count_residues(sequence, target_residues)
  residue_data$protein <- protein_name
  residue_data$length <- length(sequence_chars)
  residue_results[[i]] <- residue_data

  window_results[[i]] <- calculate_window_density(
    sequence = sequence,
    protein_name = protein_name,
    window_size = window_size,
    residues = target_residues
  )

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

human_residue_mapping <- bind_rows(residue_results)
human_sliding_window_density <- bind_rows(window_results)
human_null_model_results <- bind_rows(null_results)

write.csv(
  human_residue_mapping,
  "results/tables/human_residue_mapping_results.csv",
  row.names = FALSE
)

write.csv(
  human_sliding_window_density,
  "results/tables/human_sliding_window_density.csv",
  row.names = FALSE
)

write.csv(
  human_null_model_results,
  "results/tables/human_null_model_results.csv",
  row.names = FALSE
)

message("Human UniProt sequence analysis completed.")
message("Proteins analysed: ", length(protein_sequences))
message("Null-model permutations per protein: ", n_permutations)
