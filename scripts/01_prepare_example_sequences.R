
# 01_prepare_example_sequences.R
# Purpose: Create a small example protein FASTA file for testing the residue mapping pipeline.
# Biological note: These are synthetic toy protein sequences for pipeline validation only.

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

example_fasta <- c(
  ">Protein_A_synthetic",
  "MKTCCGHHKLLACDEFGHIKLMNPQRSTVWY",
  ">Protein_B_synthetic",
  "GGGGKKKKHHHHCCCCAAAAVVVVLLLLMMMM",
  ">Protein_C_synthetic",
  "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGV"
)

writeLines(example_fasta, "data/raw/example_proteins.fasta")

message("Example FASTA written to data/raw/example_proteins.fasta")

