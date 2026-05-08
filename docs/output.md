# Output files

## Multi-species design output structure

```
output_dir/
├── sep_5s.fasta                  # Separated 5S rRNA sequences
├── sep_5_8s.fasta                # Separated 5.8S rRNA sequences
├── sep_18s.fasta                 # Separated 18S rRNA sequences
├── sep_28s.fasta                 # Separated 28S rRNA sequences
├── all_rrna.fasta                # All concatenated rRNA sequences
├── coverage_summary_results.txt  # Coverage summary across all configs
│
├── {rna_type}_{step}_{identity}_{width}/
│   ├── sliding.fasta             # Sliding-window probe fragments
│   ├── cdhit_merge.fasta         # CD-HIT-EST clustered probes
│   ├── cdhit_final.fasta         # Final probes after internal removal
│   ├── cdhit_final.alig.out      # BLAST alignment results
│   ├── rrna_probes.fasta         # Reverse-complement probes (RNA-to-DNA)
│   ├── rRNA_coverage_summary.txt # Coverage per target sequence
│   ├── probe_groups.tsv          # Probe-to-phylogenetic-level mapping
│   ├── probe_filtered.tsv        # Filtered highest-level probe assignments
│   ├── fungal_tree_report.txt    # Taxonomic decision tree
│   └── {type}_{level}_grouped.tsv   # Coverage stats per taxonomic level
│
├── coverage/                     # Per-strain coverage data
├── index/                        # Bowtie2/BLAST indices
│   ├── rrna/                     # rRNA reference indices
│   └── genome/                   # Genome indices
├── rrna_state.tsv                # rRNA mapping statistics
└── genome_state.tsv              # Genome mapping statistics
```

## Single-species design output

| File | Description |
|------|-------------|
| `{species}.rrna.fa` | Deduplicated rRNA sequences with Silva annotations |
| `{species}_sliding.fa` | Sliding-window probe fragments |
| `{species}_cdhit_merge.fa` | CD-HIT-EST clustered non-redundant probes |
| `{species}_rrna_probes.fa` | **Final probes**: dimer-filtered, reverse-complemented |

## Probe evaluation output

| Directory/File | Description |
|----------------|-------------|
| `{prefix}/0_clean/` | Quality-controlled FastQ files (fastp output) |
| `{prefix}/1_rrna/` | rRNA mapping results (Bowtie2 SAM/BAM/depth per type) |
| `{prefix}/2_aln/` | Whole-genome mapping results (Bowtie2 sorted BAM) |
| `coverage/*.probe.tsv` | Probe position information (start/end per target) |

## Coverage summary format

`coverage_summary_results.txt` columns:

| Column | Description |
|--------|-------------|
| `full_class_name` | Config identifier (e.g., `5S_40_0.8_40`) |
| `high_class_name` | rRNA type (e.g., `5S`) |
| `probe_num` | Number of final probes |
| `q1` | First quartile of per-target coverage |
| `q2` | Median coverage |
| `q3` | Third quartile of coverage |
| `score` | Composite score = `(1 / probe_num) × q3` |
