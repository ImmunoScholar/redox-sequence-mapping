# 06_download_uniprot_human_fasta.R
# Purpose:
# Download reviewed human protein FASTA entries from UniProt accessions listed in
# metadata/human_protein_accessions.csv.
#
# Data provenance:
# Protein sequences are retrieved from UniProt by accession using the UniProt REST FASTA endpoint.
# No protein sequences are manually typed into this repository.

library(readr)
library(stringr)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

accession_file <- "metadata/human_protein_accessions.csv"
output_fasta <- "data/raw/human_redox_relevant_proteins_uniprot.fasta"

accessions <- read_csv(
  accession_file,
  show_col_types = FALSE
)

required_columns <- c("accession", "gene", "protein_name", "rationale")

missing_columns <- setdiff(required_columns, colnames(accessions))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

download_fasta <- function(accession) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/",
    accession,
    ".fasta"
  )

  fasta_text <- readLines(url, warn = FALSE)

  if (length(fasta_text) == 0) {
    stop(paste("No FASTA returned for accession:", accession))
  }

  if (!startsWith(fasta_text[1], ">")) {
    stop(paste("Invalid FASTA returned for accession:", accession))
  }

  fasta_text
}

all_fasta <- list()

for (acc in accessions$accession) {
  message("Downloading UniProt FASTA: ", acc)
  all_fasta[[acc]] <- download_fasta(acc)
}

fasta_lines <- unlist(
  lapply(
    all_fasta,
    function(x) c(x, "")
  )
)

writeLines(
  fasta_lines,
  output_fasta
)

message("Downloaded FASTA written to: ", output_fasta)
message("Number of accessions processed: ", nrow(accessions))