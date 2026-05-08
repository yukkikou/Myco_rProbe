# About

Myco_rProbe is a bioinformatics pipeline for designing and evaluating fungal rRNA-targeted hybridization probes. It uses a combination of computational prediction, sequence alignment, and phylogenetic analysis to generate probe sets with optimal coverage and specificity.

## Overview

The pipeline takes fungal genome assemblies as input, predicts rRNA sequences using Barrnap, annotates them against the Silva reference database, and generates probe candidates via sliding-window fragmentation and CD-HIT-EST clustering. The multi-species workflow additionally evaluates probe coverage across taxonomic levels and automatically selects optimal probe configurations.

## Features

- **Single-species probe design** — Design rRNA probes for a single fungal species with dimer filtering
- **Multi-species (pan-fungal) probe design** — Design probes targeting conserved rRNA regions across multiple fungal species with phylogenetic core/plus grouping
- **Circular RNA (circRNA) identification** — Identify and quantify circRNAs using CIRIquant
- **In silico probe evaluation** — Evaluate probe specificity and coverage using RNA-seq mapping data
- **Taxonomic decision tree** — Generate hierarchical reports of probe coverage across taxonomic levels
- **Optimized probe selection** — Automatically select best probe configurations based on coverage and count scores
- **Core vs. plus probe grouping** — Differentiate probes targeting conserved (core) vs. variable (plus) regions

## Workflows

### Single-species design
Genome → Barrnap (rRNA prediction) → Silva DB annotation → Deduplication → Sliding-window probe fragments → CD-HIT-EST clustering → Primer3 dimer filter → Reverse-complement final probes

### Multi-species design
Genome list → Per-strain Barrnap + Silva → rRNA type separation (5S/5.8S/18S/28S) → Parameter sweep over probe lengths and thresholds → Sliding-window → CD-HIT-EST → BLAST coverage evaluation → Phylogenetic grouping → Config scoring → Decision tree report

### Evaluation
RNA-seq reads → fastp QC → Bowtie2 mapping (rRNA + genome) → Depth statistics → Coverage visualization

## License

Apache License 2.0
