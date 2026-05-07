# Myco_rProbe Parameter Reference

This document provides detailed descriptions of all parameters, input file formats, and output files for the Myco_rProbe pipeline.

---

## 1. myco_rprobe_single.sh — Single-species Probe Design

### Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `-s` | String | Yes | — | Species Latin name with underscores (e.g., `Aspergillus_fumigatus`). Used for FASTA header annotation and Silva database search. |
| `-g` | Path | Yes | — | Genome FASTA file path. Supports gzipped (`.gz`) and plain FASTA formats. |
| `-t` | Integer | Yes | — | Number of CPU threads for parallel execution. Must be a positive integer. |
| `-o` | Path | Yes | — | Output directory for intermediate and final results. Will be created if it does not exist. |
| `-m` | Path | Yes | — | Path to Singularity `.sif` image file containing bioinformatics tools (barrnap, seqkit, cd-hit-est, python). |
| `-l` | Integer | No | 40 | Probe length in nucleotides. |
| `-w` | Integer | No | 40 | Sliding window step size in nucleotides. Recommended to match probe length. |
| `-c` | Float | No | 0.8 | CD-HIT-EST sequence identity threshold for clustering (range: 0.0–1.0, recommended: 0.8–0.9). |
| `-h` | Flag | No | — | Display help message and exit. |

### Input Files

| File | Description |
|------|-------------|
| Genome FASTA | Fungal genome assembly in FASTA format (plain or `.gz` compressed). |
| `data/fungi_sliva_18S.fasta` | SILVA 18S rRNA reference database for fungi (shipped with pipeline). |
| `data/fungi_sliva_28S.fasta` | SILVA 28S rRNA reference database for fungi (shipped with pipeline). |

### Output Files

| File | Description |
|------|-------------|
| `{species}.barrnap.rrna.fa` | Raw rRNA sequences predicted by Barrnap (intermediate, removed after processing). |
| `{species}.barrnap.gff` | Barrnap rRNA annotation in GFF format (intermediate, removed). |
| `{species}.rrna.fa` | Deduplicated rRNA sequences with Silva database integration, with species-name-prefixed headers. |
| `{species}_sliding.fa` | Sliding-window probe fragments generated from merged rRNA sequences. |
| `{species}_cdhit_merge.fa` | CD-HIT-EST clustered non-redundant probes. |
| `{species}_rrna_probes.fa` | **Final output**: Dimer-filtered, ready-to-use reverse-complement probe sequences. |

---

## 2. myco_rprobe_multi.sh — Multi-species Probe Design

### Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `-g` | Path | Yes | — | Genome list file (TSV: genome_path, genus, species_id). Each line represents one strain. |
| `-t` | Integer | Yes | — | Number of CPU threads. |
| `-o` | Path | Yes | — | Output directory for multi-species probe design results. |
| `-m` | Path | Yes | — | Singularity image path. |
| `-r` | Path | No | `data/rRNA_list.tsv` | rRNA configuration file (TSV: rna_type, sliding_step, consistency, width). Defines all probe parameter combinations to test. |
| `-h` | Flag | No | — | Display help message. |

### Genome List Format (`-g`)

Tab-separated, 3 columns, no header:
```
/path/to/genome1.fa.gz  GenusA  GenusA_species1
/path/to/genome2.fa.gz  GenusB  GenusB_species2
```

| Column | Description |
|--------|-------------|
| 1 | Full path to genome FASTA (supports `.gz`) |
| 2 | Genus name (e.g., `Aspergillus`) |
| 3 | Species identifier used for file naming (e.g., `Aspergillus_fumigatus`) |

### rRNA Config Format (`-r`)

Tab-separated, 4 columns, no header:
```
5S   40   0.8   40
5S   50   0.85  50
18S  59   0.8   59
28S  59   0.8   59
```

| Column | Description |
|--------|-------------|
| 1 | rRNA type: `5S`, `5_8S`, `18S`, or `28S` |
| 2 | Sliding window step size (nt). Determines probe density. |
| 3 | CD-HIT-EST identity threshold. Higher values yield more specific probe clusters. |
| 4 | Probe width (nt). Typically equals sliding step for non-overlapping probes. |

Note: The output directory name encodes these parameters as `{type}_{step}_{identity}_{width}` (e.g., `5S_40_0.8_40`).

### Output Files

| Output | Description |
|--------|-------------|
| `all_rrna.fasta` | Concatenated rRNA sequences from all processed strains (removed after separation). |
| `sep_5s.fasta` | 5S rRNA sequences separated by seqkit grep. |
| `sep_5_8s.fasta` | 5.8S rRNA sequences separated by seqkit grep. |
| `sep_18s.fasta` | 18S rRNA sequences separated by seqkit grep. |
| `sep_28s.fasta` | 28S rRNA sequences separated by seqkit grep. |
| `coverage_summary_results.txt` | Coverage statistics and best-configuration summary for all probe sets. |

**Per-configuration subdirectory** (`{type}_{step}_{identity}_{width}/`):

| File | Description |
|------|-------------|
| `sliding.fasta` | Sliding-window fragments before clustering. |
| `cdhit_merge.fasta` | After CD-HIT-EST clustering. |
| `cdhit_final.fasta` | After internal overlapping probe removal. |
| `cdhit_final.alig.out` | BLAST alignment results (probes vs. rRNA reference). |
| `rRNA_coverage_summary.txt` | Coverage summary per target sequence. |
| `rrna_probes.fasta` | Final reverse-complement probe sequences (RNA-to-DNA converted). |
| `probe_groups.tsv` | Probe-to-phylogenetic-level mapping. |
| `probe_filtered.tsv` | Filtered probe groups with highest taxonomic level per probe. |
| `fungal_tree_report.txt` | Hierarchical taxonomic decision tree with probe counts and coverage. |
| `{type}_{level}_grouped.tsv` | Coverage statistics at a specific taxonomic level (e.g., `28S_class_grouped.tsv`). |

### Coverage Summary Format

`coverage_summary_results.txt` columns:

| Column | Description |
|--------|-------------|
| `full_class_name` | Full configuration name (e.g., `5S_40_0.8_40`). |
| `high_class_name` | High-level class (rRNA type extracted from config name, e.g., `5S`). |
| `probe_num` | Number of probes in the final probe set. |
| `q1` | First quartile of per-target coverage ratios. |
| `q2` | Median (second quartile) of coverage ratios. |
| `q3` | Third quartile of coverage ratios. |
| `score` | Composite score = `(1 / probe_num) * q3`. Higher is better. |

---

## 3. myco_rprobe_evaluate.sh — Probe Evaluation

### Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `-g` | Path | Yes | — | Library config TSV file (strain, prefix, R1, R2, genome, gtf). |
| `-r` | Path | Yes | — | Path to rRNA directory containing `sep_*.fasta` files (output from `myco_rprobe_multi.sh`). |
| `-t` | Integer | Yes | — | Number of CPU threads. |
| `-o` | Path | Yes | — | Output directory. |
| `-m` | Path | Yes | — | Singularity image path. |
| `-h` | Flag | No | — | Display help message. |

### Library Config Format

Tab-separated, 6 columns, no header:
```
Candida_albicans  P1  /path/to/R1.fq.gz  /path/to/R2.fq.gz  /path/to/genome.fa  /path/to/annotation.gtf
```

| Column | Description |
|--------|-------------|
| 1 | Strain/species name (used for file naming) |
| 2 | Sample prefix (used for directory naming) |
| 3 | Path to R1 FastQ (raw reads) |
| 4 | Path to R2 FastQ (raw reads) |
| 5 | Path to genome FASTA |
| 6 | Path to genome annotation GTF |

### Output Structure

| Directory/File | Description |
|----------------|-------------|
| `{prefix}/0_clean/` | Quality-controlled FastQ files (fastp output). |
| `{prefix}/1_rrna/` | rRNA mapping results (Bowtie2 SAM/BAM/depth per rRNA type). |
| `{prefix}/2_aln/` | Whole-genome mapping results (Bowtie2 sorted BAM). |
| `index/rrna/` | Bowtie2 indices for each strain's rRNA references. |
| `index/genome/` | Bowtie2 indices for each strain's genome. |
| `coverage/` | Coverage summary files per rRNA type + per-strain depth data. |
| `coverage/*.probe.tsv` | Probe position information (start/end per target). |
| `coverage/*.length` | Reference sequence lengths. |
| `rrna_state.tsv` | Aggregated rRNA mapping rates across all strains. |
| `genome_state.tsv` | Aggregated genome mapping rates across all strains. |
| `{prefix}/1_rrna/{strain}.{rna_type}.depth` | Per-base sequencing depth for each rRNA target. |
| `{prefix}/1_rrna/{strain}.{rna_type}.log` | Bowtie2 mapping log with overall alignment rate. |

---

## 4. myco_rprobe_circiden.sh — circular RNA Identification

### Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `-g` | Path | Yes | — | Library config TSV (same format as evaluation). |
| `-t` | Integer | Yes | — | Number of CPU threads. |
| `-o` | Path | Yes | — | Output directory. |
| `-m` | Path | Yes | — | Singularity image (CIRIquant container). |
| `-h` | Flag | No | — | Display help message. |

### Output Structure

| Directory/File | Description |
|----------------|-------------|
| `{prefix}/3_circ/` | CIRIquant output directory with circRNA detection results. |
| `{prefix}/{strain}_config.yml` | CIRIquant configuration file (auto-generated). |
| `index/genome/{prefix}.*` | BWA and HISAT2 genome indices. |

---

## 5. myco_rprobe_circseq.sh — circRNA Sequence Extraction

### Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `-g` | Path | Yes | — | circRNA config TSV (genome, gtf). |
| `-t` | Integer | Yes | — | Number of CPU threads. |
| `-o` | Path | Yes | — | Output directory. |
| `-m` | Path | Yes | — | Singularity image path. |
| `-h` | Flag | No | — | Display help message. |

### Input Format

Tab-separated, 2 columns:
```
/path/to/genome.fa  /path/to/circ_annotation.gtf
```

### Output Files

| File | Description |
|------|-------------|
| `{prefix}.circfl.fa` | Extracted circRNA sequences (via bedtools getfasta). |

---

## 6. Python Script Parameters

### remove_dimer_para.py

```
Usage: python remove_dimer_para.py <input_fasta> <output_fasta> <min_dimer_dG> <max_workers> <max_n_proportion>
```

| Argument | Type | Description |
|----------|------|-------------|
| `input_fasta` | Path | Input FASTA with probe sequences. |
| `output_fasta` | Path | Output FASTA filtered for dimers. |
| `min_dimer_dG` | Float | Minimum dimerization free energy threshold (kcal/mol). Probes with ΔG below this are removed. |
| `max_workers` | Integer | Number of parallel worker threads. |
| `max_n_proportion` | Float | Maximum proportion of ambiguous bases ('N') allowed (0.0–1.0). |

### remove_inner.py

```
Usage: python remove_inner.py <input_file> <output_file>
```

| Argument | Type | Description |
|----------|------|-------------|
| `input_file` | Path | Input FASTA with sliding-window probes. |
| `output_file` | Path | Output FASTA with overlapping probes removed. |

### parse_coverage.py

```
Usage: python parse_coverage.py --blast <blast_file> --fasta <fasta_file>
                                --species_18S <18S_db> --species_28S <28S_db>
                                --outputdir <output_dir>
```

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `--blast` | Path | Yes | BLAST alignment results (outfmt 6). |
| `--fasta` | Path | Yes | FASTA file for reference sequence lengths. |
| `--species_18S` | Path | No | SILVA 18S database (for species name lookup). |
| `--species_28S` | Path | No | SILVA 28S database (for species name lookup). |
| `--outputdir` | Path | Yes | Output directory for coverage files and plots. |

### group_probes.py

```
Usage: python group_probes.py <alignment_file> <fasta_file> <taxonomy_file>
                              <output_file> <filter_output_file>
                              [--ratio_threshold FLOAT] [--total_threshold INT]
                              [--num_workers INT]
```

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `alignment_file` | Path | Yes | — | BLAST alignment results (TSV: probe_id, target_id, ...). |
| `fasta_file` | Path | Yes | — | FASTA with semicolon-separated phylogenetic lineage in headers. |
| `taxonomy_file` | Path | Yes | — | Taxonomy info file (TSV: tax_id, rank, name). |
| `output_file` | Path | Yes | — | Output TSV for all probe-to-taxon mappings. |
| `filter_output_file` | Path | Yes | — | Filtered output with highest taxonomic level per probe. |
| `--ratio_threshold` | Float | No | 0.7 | Minimum hit ratio for a probe to be assigned to a taxon. |
| `--total_threshold` | Int | No | 0 | Minimum total species at a level for filtering. |
| `--num_workers` | Int | No | CPU count | Number of parallel workers. |

### group_probes_plus.py

```
Usage: python group_probes_plus.py <alignment_file> <fasta_file>
                                   <taxonomy_file> <taxonomy_nodes_file>
                                   <output_file> <filter_output_file>
                                   [--ratio_threshold FLOAT]
                                   [--total_threshold INT]
                                   [--num_workers INT]
```

Additional argument compared to `group_probes.py`:

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `taxonomy_nodes_file` | Path | Yes | NCBI taxonomy nodes file (TSV: tax_id \|\| parent_id \|\| rank \|\| ...). Enables lineage resolution. |

### grouped_coverage.py

```
Usage: python grouped_coverage.py <filter_output_file> <hit_threshold>
                                  <classification_level> <fasta_file>
                                  <blast_file> <output_file>
                                  <probe_fasta_file> <output_probe_prefix>
```

| Argument | Type | Description |
|----------|------|-------------|
| `filter_output_file` | Path | Filtered probe groups from group_probes.py. |
| `hit_threshold` | Float | HitRatio filter threshold. |
| `classification_level` | String | Taxonomic rank to classify probes by (e.g., `phylum`, `class`). |
| `fasta_file` | Path | Background reference FASTA. |
| `blast_file` | Path | BLAST results file. |
| `output_file` | Path | Output TSV for coverage summary. |
| `probe_fasta_file` | Path | FASTA file with all probe sequences. |
| `output_probe_prefix` | String | Prefix for per-group probe FASTA output files. |

### grouped_coverage_plus.py

```
Usage: python grouped_coverage_plus.py <filter_output_file> <hit_threshold>
                                       <classification_level_rank> <fasta_file>
                                       <blast_file> <names_dmp_file>
                                       <nodes_dmp_file> <probe_fasta_file>
                                       <output_tsv_file> <output_probe_prefix>
```

| Argument | Type | Description |
|----------|------|-------------|
| `names_dmp_file` | Path | NCBI names.dmp file for taxonomy name resolution. |
| `nodes_dmp_file` | Path | NCBI nodes.dmp file for taxonomy hierarchy. |

### decision_tree.py

```
Usage: python decision_tree.py --nodes <nodes_file> --info <info_file>
                               --grouped-dir <grouped_dir>
                               --probes-fasta <probes_fasta>
                               [--species-map-fasta <map_fasta>]
                               [--fungi-root-name <name>]
                               [--output-file <output_file>]
```

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `--nodes` | Path | Yes | — | Taxonomy nodes file (TSV). |
| `--info` | Path | Yes | — | Taxonomy info/names file (TSV). |
| `--grouped-dir` | Path | Yes | — | Directory containing `*_grouped.tsv` files. |
| `--probes-fasta` | Path | Yes | — | FASTA file containing all probe sequences. |
| `--species-map-fasta` | Path | No | — | Optional FASTA mapping probe IDs to species lineages. |
| `--fungi-root-name` | String | No | `Fungi` | Root taxon name to start tree from. |
| `--output-file` | Path | No | stdout | Output file path for the decision tree report. |

### s02_mismatch_pileup.py (plus_design)

```
Usage: python s02_mismatch_pileup.py <pileup_file> <depth_file> <probe_length>
```

| Argument | Type | Description |
|----------|------|-------------|
| `pileup_file` | Path | Samtools mpileup output. |
| `depth_file` | Path | Per-base depth file. |
| `probe_length` | Integer | Desired probe length for output probes. |

---

## 7. SLURM Script Parameters

### slurm_circ.sh

| Variable | Description |
|----------|-------------|
| `outdir` | Base output directory (set inside script). |
| `saminfo` | Path to library config TSV. |
| `--array` | SLURM array range (e.g., `20` for 1–20 tasks). |

### slurm_evaluate.sh

| Variable | Description |
|----------|-------------|
| `outdir` | Base output directory (set inside script). |
| `saminfo` | Path to library config TSV. |
| `probe_dir` | Path to rRNA probe directory from multi-species design. |
| `--array` | SLURM array range (e.g., `1-16%4` for 16 tasks, 4 concurrent). |

### slurm_circ_homo.sh

| Variable | Description |
|----------|-------------|
| `outdir` | Base output directory (set inside script). |
| `saminfo` | Path to library config TSV. |
| `--array` | SLURM array range. |

---

## 8. Data Files Reference

| File | Format | Description |
|------|--------|-------------|
| `data/fungi_sliva_18S.fasta` | FASTA | SILVA fungal 18S rRNA reference sequences with phylogenetic lineage in headers. |
| `data/fungi_sliva_28S.fasta` | FASTA | SILVA fungal 28S rRNA reference sequences with phylogenetic lineage in headers. |
| `data/fungi_taxonomy_info.tsv` | TSV | Fungal taxonomy information (tax_id, rank, name). |
| `data/fungi_taxonomy_nodes.tsv` | TSV | Fungal taxonomy node relationships (tax_id \|\| parent_id \|\| rank \|\| ...). |
| `data/rRNA_list.tsv` | TSV | Default rRNA probe configuration (length 40–59 nt). |
| `data/rRNA_list2.tsv` | TSV | Extended rRNA probe configuration (length 55–65 nt, finer granularity). |
| `data/rRNA_list.sh` | Shell | Helper script for rRNA list generation. |
