
# 03_sliding_window_analysis.R
# Purpose:
# Calculate local sequence-level density of Cys, Lys, and His residues
# using a sliding window approach.
#
# Interpretation:
# High-density windows indicate residue clustering in primary sequence only.
# They do not prove oxidative modification or biochemical reactivity.

library(Biostrings)
library(tidyverse)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fasta_file <- "data/raw/example_proteins.fasta"
protein_sequences <- readAAStringSet(fasta_file)

target_residues <- c("C", "K", "H")
window_size <- 10

calculate_window_density <- function(sequence, protein_name, window_size, residues) {

  sequence_string <- as.character(sequence)
  sequence_chars <- strsplit(sequence_string, "")[[1]]
  sequence_length <- length(sequence_chars)

  if (sequence_length < window_size) {
    return(NULL)
  }

  window_results <- list()

  for (start_position in 1:(sequence_length - window_size + 1)) {

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

all_windows <- list()

for (i in seq_along(protein_sequences)) {
  all_windows[[i]] <- calculate_window_density(
    sequence = protein_sequences[[i]],
    protein_name = names(protein_sequences)[i],
    window_size = window_size,
    residues = target_residues
  )
}

window_density_results <- bind_rows(all_windows)

write.csv(
  window_density_results,
  "results/tables/sliding_window_density.csv",
  row.names = FALSE
)

message("Sliding-window density analysis completed.")
message("Results saved to results/tables/sliding_window_density.csv")

