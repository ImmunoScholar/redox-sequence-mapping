# Computational Mapping of Oxidative Stress-Sensitive Residues in Human Proteins

## Overview

This repository contains an R-based sequence-level analytical pipeline for mapping amino acid residues associated with oxidative stress-mediated protein modification.

The pipeline focuses on:

- cysteine (Cys, C)

- lysine (Lys, K)

- histidine (His, H)

These residues are biologically relevant because they can participate in oxidative stress-associated or electrophile-associated protein modification depending on residue chemistry, structural exposure, local microenvironment, and cellular context.

## Scientific Positioning

This project is deliberately constrained to **primary sequence-level analysis**.

It does not attempt to predict biochemical reactivity. Instead, it quantifies residue composition and local sequence clustering as a reproducible first-pass strategy for hypothesis prioritisation.

The central question addressed is:

> Are oxidative stress-associated residues locally enriched in protein primary sequences beyond what would be expected from amino acid composition alone?

## Interpretation Framework: Inference vs Evidence

This project separates computational inference from biochemical evidence.

### What this pipeline provides

The pipeline provides evidence for:

- sequence-level abundance of Cys, Lys, and His residues

- local clustering of these residues along primary protein sequences

- comparison of observed clustering against randomized sequence controls

- prioritisation of sequence regions for downstream structural or experimental investigation

### What this pipeline does not provide

The pipeline does **not** provide evidence for:

- actual oxidative modification

- electrophile adduction

- redox reactivity

- covalent modification probability

- solvent accessibility

- residue pKa effects

- protein functional perturbation

- disease relevance without additional biological validation

Therefore, outputs should be interpreted as **sequence-derived prioritisation signals**, not as direct biochemical predictions.

## Interpretation and Evidence Boundary

This pipeline provides a reproducible sequence-level framework for identifying local enrichment of oxidative stress-associated residues in protein primary sequences.

In the current example dataset, the analysis shows that selected proteins contain local Cys/Lys/His clustering patterns that exceed composition-preserving randomized sequence expectations. This supports the use of the pipeline as a hypothesis-prioritisation tool for identifying sequence regions that may merit downstream structural or experimental evaluation.

Importantly, the outputs should be interpreted as residue-enrichment and clustering signals, not as direct evidence of oxidative modification or biochemical reactivity.

Claims about actual oxidation, electrophile adduction, solvent accessibility, residue pKa, or functional protein perturbation would require orthogonal evidence from structural biology, redox proteomics, targeted mutagenesis, mass spectrometry, or biochemical assays.

## Pipeline Modules

1. Generate example protein FASTA input

2. Quantify Cys/Lys/His residue counts and frequencies

3. Perform sliding-window residue density analysis

4. Compare observed clustering against composition-preserving randomized sequence null models

5. Generate interpretable visual outputs

## Repository Structure

```text

redox-sequence-mapping/

├── data/

│   ├── raw/

│   └── processed/

├── docs/

├── environment/

├── metadata/

├── notebooks/

├── results/

│   ├── figures/

│   └── tables/

├── scripts/

└── README.md

