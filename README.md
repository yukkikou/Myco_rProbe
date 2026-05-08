# Myco_rProbe

**Fungal rRNA Probe Design and Evaluation Pipeline**

Myco_rProbe is a comprehensive bioinformatics pipeline for designing, evaluating, and optimizing species-specific and pan-fungal rRNA-targeted hybridization probes. It supports both single-species and multi-species probe design, circular RNA (circRNA) identification, and in silico probe evaluation using RNA-seq data.

---

## Features

- **Single-species probe design**: Design rRNA probes for a single fungal species
- **Multi-species (pan-fungal) probe design**: Design probes targeting conserved rRNA regions across multiple fungal species with phylogenetic grouping
- **Circular RNA (circRNA) identification**: Identify and quantify circRNAs using CIRIquant
- **In silico probe evaluation**: Evaluate probe specificity and coverage using RNA-seq mapping data
- **Taxonomic decision tree**: Generate hierarchical reports of probe coverage across taxonomic levels
- **Optimized probe selection**: Automatically select optimal probe configurations based on coverage and count scores
- **Core vs. plus probe grouping**: Differentiate probes targeting conserved (core) vs. variable (plus) regions

---

## Pipeline Overview

The pipeline consists of three major workflows. Below is a step-by-step breakdown of each.

### A. Single-species Probe Design (`myco_rprobe_single.sh`)

```
Genome FASTA ──► Barrnap (rRNA prediction) ──► Silva DB (species annotation via seqkit grep)
    │
    ▼
rRNA sequences merged & deduplicated ──► Sliding window (probe fragment generation)
    │
    ▼
CD-HIT-EST (redundancy clustering, -c threshold) ──► Dimer filter (primer3, ΔG cutoff) ──► Final probes
```

1. **rRNA prediction** — Barrnap identifies 5S, 5.8S, 18S, and 28S rRNA sequences from the input genome.
2. **Silva DB integration** — Annotation-aware seqkit grep matches predicted rRNA to the Silva reference database, adding species lineage info to FASTA headers.
3. **Deduplication** — Redundant sequences are removed by sequence identity (seqkit rmdup).
4. **Sliding window** — Overlapping probe fragments of length `-l` are generated at step `-w`.
5. **CD-HIT-EST clustering** — Redundant probe sequences are clustered at identity threshold `-c`.
6. **Dimer filter** — Primer3 evaluates dimerization free energy (ΔG); probes below the threshold are removed via `remove_dimer_para.py`.
7. **Reverse complement** — Final probe sequences are reverse-complemented for experimental use.

### B. Multi-species (Pan-fungal) Probe Design (`myco_rprobe_multi.sh`)

```
Genome list (TSV) ──► Per-strain: Barrnap + Silva DB ──► All rRNA merged
    │
    ▼
rRNA separation by type (5S / 5.8S / 18S / 28S) via seqkit grep
    │
    ▼
For each rRNA config line (type, step, identity, width):
    │
    ├── Sliding window (step = configured step, width = configured width)
    ├── CD-HIT-EST clustering (threshold = configured identity)
    ├── Internal probe removal (remove_inner.py — strips overlapping probes)
    ├── BLAST against Silva reference (coverage evaluation)
    ├── Coverage parsing (parse_coverage.py → rRNA_coverage_summary.txt)
    ├── Probe phylogenetic grouping (group_probes.py / group_probes_plus.py)
    │   └── Core (conserved) vs Plus (variable) assignment
    └── Coverage statistics by taxonomic level (grouped_coverage.py)
    │
    ▼
Score = (1 / probe_num) × q3 ──► Best config auto-selected by rank
    │
    ▼
Decision tree (decision_tree.py) — hierarchical report across all ranks:
    subkingdom → phylum → subphylum → class → subclass → order → suborder → family → subfamily
```

1. **Multi-strain rRNA extraction** — Each strain in the genome list is processed through Barrnap and Silva-annotated, then all rRNA sequences are concatenated.
2. **Type separation** — `all_rrna.fasta` is split into `sep_5s.fasta`, `sep_5_8s.fasta`, `sep_18s.fasta`, `sep_28s.fasta`.
3. **Parameter sweep** — The rRNA config (`-r`) defines multiple probe design parameter combinations (e.g., `5S_40_0.8_40`, `28S_59_0.8_59`). Each combination is processed independently in its own subdirectory.
4. **Probe generation & clustering** — Same sliding-window, CD-HIT-EST, and internal-removal steps as single-species.
5. **Coverage evaluation** — BLASTN alignment against the Silva reference, coverage calculated per target sequence.
6. **Phylogenetic grouping** — Probes are assigned to taxonomic levels based on BLAST hit distribution. Core probes target conserved (high-hit-ratio) regions; plus probes target variable (low-hit-ratio) regions.
7. **Config scoring** — Each config is scored as `(1 / probe_num) × q3`. The best configuration (highest score) per rRNA type is automatically selected.
8. **Decision tree** — A hierarchical report shows probe counts and coverage across all taxonomic ranks.

### C. Probe Evaluation (`myco_rprobe_evaluate.sh`)

```
RNA-seq (R1/R2) ──► fastp (QC/trim) ──► Clean reads
    │
    ├── Bowtie2 → rRNA reference (per type: 5S/5.8S/18S/28S) → depth stats
    │       └── rRNA mapping rate → rrna_state.tsv
    │
    ├── Bowtie2 → Whole genome → sorted BAM
    │       └── Genome mapping rate → genome_state.tsv
    │
    └── Coverage plots (R: stat_coverage.R)
```

1. **Quality control** — fastp trims adapters and filters low-quality reads.
2. **rRNA mapping** — Clean reads are mapped to each rRNA reference (Bowtie2). Per-base depth and alignment rates are recorded.
3. **Genome mapping** — Reads are mapped to the whole genome for background evaluation.
4. **Coverage statistics** — Aggregated rRNA/genome mapping rates are written to `rrna_state.tsv` and `genome_state.tsv`.
5. **Visualization** — R script generates per-strain coverage plots.

### D. circRNA Identification (`myco_rprobe_circiden.sh`)

```
Leftover reads (non-rRNA) ──► BWA index + HISAT2 index ──► CIRIquant ──► circRNA GTF / BAM
```

### E. circRNA Sequence Extraction (`myco_rprobe_circseq.sh`)

```
circRNA GTF + Genome FASTA ──► bedtools getfasta ──► circRNA FASTA
```

---

## Requirements

### Software
- [Singularity](https://docs.sylabs.io/) (v3.x) for containerized tool execution
- Required container:
  - `singularity/pan_fungi.sif` - Consolidated image with all bioinformatics tools (barrnap, seqkit, cd-hit-est, BLAST, Bowtie2, Samtools, fastp, HISAT2, BWA, Python environment with Biopython/primer3-py/pandas/matplotlib/numpy, R with tidyverse/purrr)

### Python Dependencies (in container)
- Biopython
- primer3-py
- pandas
- matplotlib
- numpy

### R Dependencies (in container)
- tidyverse
- purrr

---

## Installation

```bash
git clone https://github.com/your-org/Myco_rProbe.git
cd Myco_rProbe
```

Ensure your Singularity image files are available at the specified paths.

---

## Usage

### 1. Single-species Probe Design

```bash
bash myco_rprobe_single.sh \
    -s Aspergillus_fumigatus \
    -g example/Aspergillus_fumigatus.ASM265v1.dna.toplevel.fa.gz \
    -t 20 \
    -o single_af \
    -m singularity/pan_fungi.sif
```

| Parameter | Description |
|-----------|-------------|
| `-s` | Species name (underscore-separated, e.g., `Aspergillus_fumigatus`) |
| `-g` | Genome FASTA file path (supports `.gz`) |
| `-t` | Number of CPU threads |
| `-o` | Output directory |
| `-m` | Singularity image path |
| `-l` | Probe length (default: 40 nt) |
| `-w` | Sliding window step (default: 40 nt) |
| `-c` | CD-HIT-EST identity threshold (default: 0.8) |

### 2. Multi-species Pan-fungal Probe Design

```bash
bash myco_rprobe_multi.sh \
    -g example/fungi_genome.list \
    -t 10 \
    -o pathgen_fungi \
    -m singularity/pan_fungi.sif \
    -r data/rRNA_list.tsv
```

**Genome list format** (TSV, no header):
```
<genome_path>  <genus>  <species_identifier>
```

**rRNA config format** (TSV):
```
<rna_type>  <sliding_step>  <consistency>  <width>
```

Example:
```
5S    40    0.8    40
5S    50    0.85   50
18S   59    0.8    59
28S   59    0.8    59
```

### 3. Probe Evaluation with RNA-seq Data

```bash
bash myco_rprobe_evaluate.sh \
    -g example/batch_lib_config.tsv \
    -t 10 \
    -o batch_mapping \
    -m singularity/pan_fungi.sif \
    -r pathgen_fungi
```

**Library config format** (TSV):
```
<strain_name>  <prefix>  <R1_path>  <R2_path>  <genome_path>  <gtf_path>
```

### 4. circRNA Identification

```bash
bash myco_rprobe_circiden.sh \
    -g example/batch_lib_config.tsv \
    -t 10 \
    -o batch_mapping \
    -m singularity/pan_fungi.sif
```

### 5. circRNA Sequence Extraction

```bash
bash myco_rprobe_circseq.sh \
    -g example/circ_seq.tsv \
    -t 5 \
    -o output_dir \
    -m singularity/pan_fungi.sif
```

---

## Output Structure

```
output_dir/
├── sep_5s.fasta              # Separated 5S rRNA sequences
├── sep_5_8s.fasta            # Separated 5.8S rRNA sequences
├── sep_18s.fasta             # Separated 18S rRNA sequences
├── sep_28s.fasta             # Separated 28S rRNA sequences
├── all_rrna.fasta            # All concatenated rRNA sequences
├── coverage_summary_results.txt  # Coverage summary
│
├── {rna_type}_{step}_{identity}_{width}/
│   ├── sliding.fasta         # Sliding-window probes
│   ├── cdhit_merge.fasta     # CD-HIT-EST merged probes
│   ├── cdhit_final.fasta     # Filtered final probes
│   ├── cdhit_final.alig.out  # BLAST alignment results
│   ├── rrna_probes.fasta     # Probe sequences (RC, RNA-to-DNA)
│   ├── rRNA_coverage_summary.txt  # Coverage per target
│   ├── probe_groups.tsv      # Probe phylogenetic grouping
│   ├── probe_filtered.tsv    # Filtered highest-level grouping
│   ├── fungal_tree_report.txt    # Taxonomic decision tree
│   └── {rna_type}_{level}_grouped.tsv  # Coverage per taxonomic level
│
├── coverage/                 # Per-strain coverage data
├── index/                    # Bowtie2/BLAST indices
│   ├── rrna/                 # rRNA reference indices
│   └── genome/               # Genome indices
├── rrna_state.tsv            # rRNA mapping statistics
└── genome_state.tsv          # Genome mapping statistics
```

---

## Key Scripts

| Script | Description |
|--------|-------------|
| `reassemble.sh` | Rebuild `singularity/pan_fungi.sif` from split parts |
| `myco_rprobe_single.sh` | Single-species probe design |
| `myco_rprobe_multi.sh` | Multi-species pan-fungal probe design |
| `myco_rprobe_evaluate.sh` | In silico probe evaluation with RNA-seq |
| `myco_rprobe_circiden.sh` | Circular RNA identification |
| `myco_rprobe_circseq.sh` | circRNA sequence extraction |
| `src/circ_blast.sh` | BLAST array job script |
| `src/plus_design/s01_run_mpile.sh` | Pileup generation for plus design |
| `src/plus_design/s02_mismatch_pileup.py` | Mismatch pileup analysis |

### Python Analysis Scripts (src/)

| Script | Description |
|--------|-------------|
| `remove_dimer_para.py` | Remove dimer-prone probes using primer3 (multithreaded) |
| `remove_inner.py` | Remove overlapping internal probes |
| `parse_coverage.py` | Parse BLAST results and calculate coverage statistics |
| `group_probes.py` | Group probes by phylogenetic levels |
| `group_probes_plus.py` | Advanced probe grouping with NCBI taxonomy |
| `grouped_coverage.py` | Calculate coverage per phylogenetic group |
| `grouped_coverage_plus.py` | Advanced coverage calculation with NCBI taxonomy |
| `decision_tree.py` | Generate hierarchical taxonomic decision tree report |
| `mismatch_pileup.py` | Analyze pileup for consensus and probe design |

---

---

## Citation

If you use Myco_rProbe in your research, please cite:

> (Citation information to be added)

---

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.

---

## Contact

For questions or support, please open an issue on GitHub.
