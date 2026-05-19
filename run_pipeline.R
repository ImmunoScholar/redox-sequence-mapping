# run_pipeline.R
# Purpose:
# Execute the full sequence-level redox residue mapping workflow.
#
# Scientific scope:
# This pipeline performs primary-sequence analysis only. It does not infer
# biochemical reactivity, oxidative modification, electrophile adduction,
# or structural accessibility.

message("Running sequence-level redox residue mapping pipeline...")

source("scripts/01_prepare_example_sequences.R")
source("scripts/02_residue_mapping.R")
source("scripts/03_sliding_window_analysis.R")
source("scripts/04_null_model_analysis.R")
source("scripts/05_visualization.R")

message("Pipeline completed successfully.")
