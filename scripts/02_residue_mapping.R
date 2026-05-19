
# 02_residue_mapping.R
# Purpose:
# Quantify oxidative stress-associated residues (Cys, Lys, His)
# across protein sequences.
#
# IMPORTANT:
# This analysis is sequence-level only and does NOT infer
# biochemical reactivity or structural accessibility.

library(Biostrings)
library(tidyverse)
library(data.table)

# Create output directory
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load FASTA
fasta_file <- "data/raw/example_proteins.fasta"

protein_sequences <- readAAStringSet(fasta_file)

# Define residues of interest
target_residues <- c("C", "K", "H")

# Function to count residues
count_residues <- function(sequence, residues) {

  sequence_chars <- strsplit(as.character(sequence), "")[[1]]

  residue_counts <- sapply(residues, function(residue) {
    sum(sequence_chars == residue)
  })

  total_length <- nchar(as.character(sequence))

  residue_frequencies <- residue_counts / total_length

  return(data.frame(
    residue = residues,
    count = residue_counts,
    frequency = residue_frequencies
  ))
}

# Analyze all proteins
results_list <- list()

for (i in seq_along(protein_sequences)) {

  protein_name <- names(protein_sequences)[i]

  residue_data <- count_residues(
    protein_sequences[[i]],
    target_residues
  )

  residue_data$protein <- protein_name

  residue_data$length <- nchar(
    as.character(protein_sequences[[i]])
  )

  results_list[[i]] <- residue_data
}

# Combine results
final_results <- bind_rows(results_list)

# Save results
write.csv(
  final_results,
  "results/tables/residue_mapping_results.csv",
  row.names = FALSE
)

message("Residue mapping completed.")
message("Results saved to results/tables/residue_mapping_results.csv")

