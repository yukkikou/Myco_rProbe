# Usage

## Command overview

Myco_rProbe provides five main entry-point scripts, each corresponding to a distinct workflow. All scripts accept a Singularity image path via `-m` and output directory via `-o`.

## Common workflow

### 1. Single-species probe design

Design probes for a single fungal species from its genome assembly.

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
| `-s` | Species name (underscore-separated) |
| `-g` | Genome FASTA (`.gz` supported) |
| `-t` | CPU threads |
| `-o` | Output directory |
| `-m` | Singularity image path |
| `-l` | Probe length (default: 40) |
| `-w` | Sliding window step (default: 40) |
| `-c` | CD-HIT-EST identity (default: 0.8) |

### 2. Multi-species pan-fungal probe design

Design probes targeting conserved rRNA regions across multiple species.

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

### 3. Probe evaluation with RNA-seq data

Evaluate probe performance using RNA-seq alignment statistics.

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

### 4. circRNA identification

Identify circular RNAs from RNA-seq data.

```bash
bash myco_rprobe_circiden.sh \
    -g example/batch_lib_config.tsv \
    -t 10 \
    -o batch_mapping \
    -m singularity/pan_fungi.sif
```

### 5. circRNA sequence extraction

Extract circRNA sequences from annotation.

```bash
bash myco_rprobe_circseq.sh \
    -g example/circ_seq.tsv \
    -t 5 \
    -o output_dir \
    -m singularity/pan_fungi.sif
```

## Utility scripts

| Script | Description |
|--------|-------------|
| `reassemble.sh` | Rebuild `pan_fungi.sif` from split parts |
| `src/circ_blast.sh` | BLAST array job for circRNA homology search |
| `src/plus_design/s01_run_mpile.sh` | Pileup generation for plus probe design |
| `src/plus_design/s02_mismatch_pileup.py` | Mismatch pileup analysis for consensus calling |
