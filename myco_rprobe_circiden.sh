#!/bin/bash
#
# Copyright 2024 Myco_rProbe
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# =============================================================================
#
# Myco_rProbe - Circular RNA identification pipeline
#
# This script identifies circular RNAs (circRNAs) from RNA-seq data using
# CIRIquant. It handles index building (BWA and HISAT2) and circRNA
# quantification for each sample.
#
# Usage:
#   bash myco_rprobe_circiden.sh -g <genome_mapping_list> -t <threads> \
#                                 -o <output_dir> -m <mirror_path>
#
# Author: Myco_rProbe Team
# Date:   2024
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Display help information
#
# Globals:
#   $0 - script name
# Returns:
#   None (exits with code 1)
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: bash $0 -g <genome_mapping_list> -t <threads> -o <output_dir> -m <mirror_path>

Required arguments:
  -g <genome_mapping_list>    Genome list file (TSV: strain<TAB>prefix<TAB>R1<TAB>R2<TAB>genome<TAB>gtf).
                              Example: 'example/batch_lib_config.tsv'
  -t <threads>                Number of CPU threads (e.g., 4).
  -o <output_dir>             Output directory path (e.g., 'mapping').
  -m <mirror_path>            Path to the Singularity image file.

Optional arguments:
  -h                          Show this help message and exit.

Examples:
  bash myco_rprobe_circiden.sh -g example/batch_lib_config.tsv -t 10 \\
      -o batch_mapping -m /path/to/ciriquant.sif
EOF
    exit 1
}

# -----------------------------------------------------------------------------
# Initialize and parse arguments
# -----------------------------------------------------------------------------
GENOME_LIST=""
THREADS=""
OUTPUT_DIR=""
MIRROR_PATH=""

while getopts ":g:t:o:m:h" opt; do
    case $opt in
        g)  GENOME_LIST="$OPTARG" ;;
        t)  THREADS="$OPTARG" ;;
        o)  OUTPUT_DIR="$OPTARG" ;;
        m)  MIRROR_PATH="$OPTARG" ;;
        h)  usage ;;
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Validate arguments
if [[ -z "$GENOME_LIST" || -z "$THREADS" || -z "$OUTPUT_DIR" || -z "$MIRROR_PATH" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

if [[ ! -f "$GENOME_LIST" ]]; then
    echo "Error: Genome list '$GENOME_LIST' does not exist." >&2
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
    echo "Error: Threads must be a positive integer." >&2
    exit 1
fi
if [[ ! -f "$MIRROR_PATH" ]]; then
    echo "Error: Singularity file '$MIRROR_PATH' does not exist." >&2
    exit 1
fi

echo "=== Parameters ==="
echo "Genome Mapping   : $GENOME_LIST"
echo "Threads          : $THREADS"
echo "Output Directory : $OUTPUT_DIR"
echo "Mirror PATH      : $MIRROR_PATH"
echo "==================="

readonly SIF_CMD="singularity -q exec -B /media --cleanenv $MIRROR_PATH"
mkdir -p "$OUTPUT_DIR"

echo "[$(date '+%F %T')]: Starting circRNA identification..."

# -----------------------------------------------------------------------------
# Build BWA and HISAT2 indices for a genome
#
# Args:
#   $1 - strain_prefix: sample prefix for index naming
#   $2 - genome_fasta: path to genome FASTA
#   $3 - index_dir: output directory for indices
#   $4 - thread: number of threads
# Returns:
#   0 on success, 1 on index build failure
# -----------------------------------------------------------------------------
build_indices() {
    local strain_prefix="$1"
    local genome_fasta="$2"
    local index_dir="$3"
    local thread="$4"
    local index_prefix="$index_dir/$strain_prefix"

    echo "[$(date '+%F %T')]: Checking indices for $strain_prefix..."
    mkdir -p "$index_dir"

    if [[ -f "${index_prefix}.bwt" ]]; then
        echo "[$(date '+%F %T')]: BWA index already exists for $strain_prefix."
    else
        echo "[$(date '+%F %T')]: Building BWA index for $strain_prefix..."
        $SIF_CMD bwa index -p "$index_prefix" "$genome_fasta" || {
            echo "[$(date '+%F %T')]: ERROR: BWA index build failed." >&2
            return 1
        }
    fi

    if [[ -f "${index_prefix}.1.ht2" ]]; then
        echo "[$(date '+%F %T')]: HISAT2 index already exists for $strain_prefix."
    else
        echo "[$(date '+%F %T')]: Building HISAT2 index for $strain_prefix..."
        $SIF_CMD hisat2-build -p "$thread" "$genome_fasta" "$index_prefix" || {
            echo "[$(date '+%F %T')]: ERROR: HISAT2 index build failed." >&2
            return 1
        }
    fi
    return 0
}

# -----------------------------------------------------------------------------
# Run CIRIquant for circRNA identification
#
# Args:
#   $1 - strain_name: strain name
#   $2 - strain_prefix: sample prefix
#   $3 - fq1: path to R1 FastQ
#   $4 - fq2: path to R2 FastQ
#   $5 - yml: path to YAML config
#   $6 - thread: number of threads
# Returns:
#   0 on success, 1 on CIRIquant failure
# -----------------------------------------------------------------------------
circrna_identification() {
    local strain_name="$1"
    local strain_prefix="$2"
    local fq1="$3"
    local fq2="$4"
    local yml="$5"
    local thread="$6"
    local circ_dir="$OUTPUT_DIR/$strain_prefix/3_circ"

    mkdir -p "$circ_dir"

    $SIF_CMD CIRIquant -t "$thread" -v \
        -1 "$fq1" \
        -2 "$fq2" \
        --config "$yml" \
        -p "${strain_name}_circRNA" \
        -o "$circ_dir" -l 2

    if [[ $? -eq 0 ]]; then
        echo "[$(date '+%F %T')]: circRNA identification completed for $strain_name."
        return 0
    else
        echo "[$(date '+%F %T')]: ERROR: circRNA identification failed for $strain_name." >&2
        return 1
    fi
}

# -----------------------------------------------------------------------------
# Main processing loop
# -----------------------------------------------------------------------------
while IFS=$'\t' read -r strain prefix raw1 raw2 genome gtf; do
    clean1="$OUTPUT_DIR/${prefix}/0_clean/${strain}_cleaned_R1.fastq.gz"
    clean2="$OUTPUT_DIR/${prefix}/0_clean/${strain}_cleaned_R2.fastq.gz"
    yml_file="$OUTPUT_DIR/${prefix}/${strain}_config.yml"

    # Check required files
    missing_files=false
    [[ ! -f "$genome" ]] && { echo "[$(date '+%F %T')]: ERROR: Genome not found: $genome" >&2; missing_files=true; }
    [[ ! -f "$gtf" ]] && { echo "[$(date '+%F %T')]: ERROR: GTF not found: $gtf" >&2; missing_files=true; }
    [[ ! -f "$clean1" ]] && { echo "[$(date '+%F %T')]: ERROR: R1 not found: $clean1" >&2; missing_files=true; }
    [[ ! -f "$clean2" ]] && { echo "[$(date '+%F %T')]: ERROR: R2 not found: $clean2" >&2; missing_files=true; }
    $missing_files && continue

    # Build genome indices
    index_base_dir="$OUTPUT_DIR/index/genome"
    index_prefix_path="$index_base_dir/$prefix"

    if ! build_indices "$prefix" "$genome" "$index_base_dir" "$THREADS"; then
        echo "[$(date '+%F %T')]: ERROR: Index building failed for $strain. Skipping." >&2
        continue
    fi

    # Generate CIRIquant YAML config
    echo "[$(date '+%F %T')]: Generating YAML config: $yml_file"
    cat > "$yml_file" << EOF
name: $prefix
tools:
  bwa: /venv/bin/bwa
  hisat2: /venv/bin/hisat2
  stringtie: /venv/bin/stringtie
  samtools: /venv/bin/samtools
reference:
  fasta: $genome
  gtf: $gtf
  bwa_index: $index_prefix_path
  hisat_index: $index_prefix_path
EOF

    # Run circRNA identification
    if circrna_identification "$strain" "$prefix" "$clean1" "$clean2" "$yml_file" "$THREADS"; then
        echo "[$(date '+%F %T')]: Strain $strain processed successfully."
    else
        echo "[$(date '+%F %T')]: ERROR: Failed to process $strain." >&2
    fi

done < <(sed 's/\r$//' "$GENOME_LIST")

echo "[$(date '+%F %T')]: CircRNA identification completed."
