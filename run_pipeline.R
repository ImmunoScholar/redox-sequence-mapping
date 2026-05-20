# run_pipeline.R
# Purpose:
# Execute the full sequence-level redox residue mapping workflow.
#
# Scientific scope:
# This pipeline performs primary-sequence analysis only. It does not infer
# biochemical reactivity, oxidative modification, electrophile adduction,
# structural accessibility, or functional perturbation.

message("Running sequence-level redox residue mapping pipeline...")

message("Running example-sequence workflow...")
source("scripts/01_prepare_example_sequences.R")
source("scripts/02_residue_mapping.R")
source("scripts/03_sliding_window_analysis.R")
source("scripts/04_null_model_analysis.R")
source("scripts/05_visualization.R")

message("Running UniProt-derived human-protein workflow...")
source("scripts/06_download_uniprot_human_fasta.R")
source("scripts/07_analyze_uniprot_human_fasta.R")
source("scripts/08_summarize_human_results.R")
source("scripts/09_visualize_human_results.R")
source("scripts/10_summarize_human_residue_classes.R")

message("Pipeline completed successfully.")