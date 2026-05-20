# Scientific Rationale

## Purpose

This document explains the scientific rationale for selecting cysteine, lysine, and histidine as residue classes for sequence-level mapping in the context of oxidative stress-mediated and electrophile-associated protein modification biology.

## Residue-class selection

The pipeline focuses on Cys, Lys, and His because these residues contain side-chain chemistries that are relevant to oxidative stress and lipid electrophile-associated protein modification contexts.

- Cysteine contains a thiol side chain and is central to many redox-regulated protein mechanisms.
- Lysine contains a primary amine side chain and is relevant to carbonyl and lipid electrophile-adduction contexts under suitable biochemical conditions.
- Histidine contains an imidazole side chain and can participate in catalysis, metal coordination, and context-dependent oxidative or electrophile-associated modification biology.

## What the sequence-level workflow measures

The workflow measures residue abundance, residue frequency, local sliding-window density, and composition-preserving randomized expectations in protein primary sequences.

The UniProt-derived human protein workflow applies the same analysis to selected reviewed Homo sapiens proteins retrieved by accession.

## Interpretation boundary

The pipeline does not infer that an individual residue is oxidized, electrophile-adducted, solvent-exposed, catalytically active, structurally reactive, or functionally perturbed.

Residue clustering should be interpreted as a sequence-level prioritisation feature only.

Biochemical claims would require orthogonal evidence such as structural analysis, redox proteomics, targeted mass spectrometry, mutagenesis, biochemical assays, or cellular functional readouts.

## Defensible interpretation

The strongest defensible interpretation is that the pipeline maps oxidative stress-associated residue classes in protein sequences and prioritises local residue-enriched regions for downstream structural or experimental evaluation.
