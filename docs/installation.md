# Installation

## Availability

Myco_rProbe is available on GitHub at [https://github.com/yukkikou/Myco_rProbe](https://github.com/yukkikou/Myco_rProbe).

## Dependencies

Install Singularity (v3.x) for containerized tool execution. All bioinformatics tools are bundled in the provided Singularity image.

## Install from source

```bash
git clone https://github.com/yukkikou/Myco_rProbe.git
cd Myco_rProbe
```

## Singularity image

The pipeline requires a consolidated Singularity image (`pan_fungi.sif`) containing all required tools. The image is provided as split archives in the repository:

```bash
# Reassemble the Singularity image from split parts
bash reassemble.sh
```

This creates `singularity/pan_fungi.sif` (1.6 GB), which includes:
- barrnap, seqkit, cd-hit-est
- BLAST, Bowtie2, Samtools, fastp, HISAT2, BWA
- Python 3 with Biopython, primer3-py, pandas, matplotlib, numpy
- R with tidyverse, purrr

## Required reference files

1. **Silva rRNA databases** — Provided in `data/`:
   - `data/fungi_sliva_18S.fasta`
   - `data/fungi_sliva_28S.fasta`

2. **Taxonomy files** — Provided in `data/`:
   - `data/fungi_taxonomy_info.tsv`
   - `data/fungi_taxonomy_nodes.tsv.gz` (decompress before use)
