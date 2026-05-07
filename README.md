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

```
                         ┌─────────────────────┐
                         │   Genome FASTA(s)    │
                         └──────────┬──────────┘
                                    │
                         ┌──────────▼──────────┐
                         │  Barrnap (rRNA      │
                         │  identification)    │
                         └──────────┬──────────┘
                                    │
                         ┌──────────▼──────────┐
                         │  Silva DB Integration│
                         └──────────┬──────────┘
                                    │
              ┌─────────────────────┼─────────────────────┐
              │                     │                     │
   ┌──────────▼──────────┐ ┌───────▼────────┐ ┌─────────▼──────────┐
   │ Single-species      │ │ Multi-species  │ │ circRNA            │
   │ Probe Design        │ │ Probe Design   │ │ Identification     │
   │ (single.sh)         │ │ (multi.sh)     │ │ (circiden.sh)      │
   └──────────┬──────────┘ └───────┬────────┘ └─────────┬──────────┘
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ CD-HIT-EST      │            │
              │           │ Clustering      │            │
              │           └────────┬────────┘            │
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ Internal Probe  │            │
              │           │ Removal         │            │
              │           └────────┬────────┘            │
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ BLAST Coverage  │            │
              │           │ Evaluation      │            │
              │           └────────┬────────┘            │
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ Probe Grouping  │            │
              │           │ (Core vs Plus)  │            │
              │           └────────┬────────┘            │
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ Evaluation      │            │
              │           │ Pipeline        │◄───────────┤
              │           │ (evaluate.sh)   │            │
              │           └────────┬────────┘            │
              │                    │                     │
              │           ┌────────▼────────┐            │
              │           │ Decision Tree   │            │
              │           │ Report          │            │
              │           └─────────────────┘            │
```

---

## Requirements

### Software
- [Singularity](https://docs.sylabs.io/) (v3.x) for containerized tool execution
- Required containers (provided as `.sif` files):
  - `pan_fungi3.sif` / `pan_fungi5.sif` - Python environment with bioinformatics tools
  - `ciriquant-v1.1.3.sif` - CIRIquant for circRNA analysis
  - `blast-2.13.0.sif` - BLAST for sequence alignment
  - `seqkit-2.2.0.sif` - FASTA/Q sequence manipulation
  - `samtools-1.14.sif` - SAM/BAM processing

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
    -m /path/to/pan_fungi3.sif
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
    -m /path/to/pan_fungi3.sif \
    -r data/rRNA_list2.tsv
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
    -m /path/to/pan_fungi5.sif \
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
    -m /path/to/ciriquant.sif
```

### 5. circRNA Sequence Extraction

```bash
bash myco_rprobe_circseq.sh \
    -g example/circ_seq.tsv \
    -t 5 \
    -o output_dir \
    -m /path/to/pan_fungi5.sif
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
| `myco_rprobe_single.sh` | Single-species probe design |
| `myco_rprobe_multi.sh` | Multi-species pan-fungal probe design |
| `myco_rprobe_evaluate.sh` | In silico probe evaluation with RNA-seq |
| `myco_rprobe_circiden.sh` | Circular RNA identification |
| `myco_rprobe_circseq.sh` | circRNA sequence extraction |
| `merge_circ_byid.sh` | Merge circRNA results by group (SLURM) |
| `slurm_circ.sh` | SLURM array job for circRNA |
| `slurm_evaluate.sh` | SLURM array job for evaluation |
| `transfer_toT640.sh` | Data transfer to remote server |
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

## SLURM Cluster Usage

For large-scale jobs on SLURM clusters:

```bash
# Array job for circRNA identification (adjust --array in slurm_circ.sh)
sbatch slurm_circ.sh

# Array job for probe evaluation
sbatch slurm_evaluate.sh

# circRNA homology search
sbatch slurm_circ_homo.sh
```

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
