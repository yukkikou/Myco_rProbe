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
# Myco_rProbe - Single-species rRNA probe design pipeline
#
# This script performs rRNA probe design for a single fungal species.
# It identifies rRNA sequences from a genome using Barrnap, integrates
# Silva database references, generates sliding-window probes, clusters
# them via CD-HIT-EST, and removes potential dimer artifacts.
#
# Usage:
#   bash myco_rprobe_single.sh -s <species_name> -g <genome_file> \
#                               -t <threads> -o <output_dir> -m <mirror_path> [options]
#
# Author: Myco_rProbe Team
# Date:   2024
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Default settings
# -----------------------------------------------------------------------------
readonly DEFAULT_PROBE_LENGTH=40    # Default probe length (nt)
readonly DEFAULT_SLIDING_WINDOW=40  # Default sliding window step (nt)
readonly DEFAULT_CONSISTENCY=0.8    # Default CD-HIT-EST sequence identity threshold

# -----------------------------------------------------------------------------
# Display help information
#
# Prints the usage guide detailing all required and optional arguments.
#
# Globals:
#   $0 - script name
# Returns:
#   None (exits with code 1)
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: bash $0 -s <species_name> -g <genome_file> -t <threads> -o <output_dir> -m <mirror_path> [options]

Required arguments:
  -s <species_name>      Species Latin name without spaces (e.g., 'Aspergillus_fumigatus').
  -g <genome_file>       Path to the genome FASTA file (e.g., 'example/Aspergillus_fumigatus.ASM265v1.dna.toplevel.fa.gz').
  -t <threads>           Number of CPU threads (e.g., 4).
  -o <output_dir>        Output directory path (e.g., 'single_af').
  -m <mirror_path>       Path to the Singularity image file (binds /media by default).

Optional arguments:
  -l <probe_length>      Probe length in nucleotides. Default: ${DEFAULT_PROBE_LENGTH}.
  -w <sliding_window>    Sliding window step size. Default: ${DEFAULT_SLIDING_WINDOW}.
  -c <consistency>       CD-HIT-EST sequence identity threshold (0.8-1.0). Default: ${DEFAULT_CONSISTENCY}.
  -h                     Show this help message and exit.

Examples:
  bash myco_rprobe_single.sh -s Aspergillus_fumigatus -g genome.fa.gz -t 20 -o single_af -m /path/to/image.sif
  bash myco_rprobe_single.sh -s Candida_albicans -g genome.fa.gz -t 10 -o single_ca -m /path/to/image.sif -l 50 -w 50 -c 0.85
EOF
    exit 1
}

# -----------------------------------------------------------------------------
# Initialize variables
# -----------------------------------------------------------------------------
SPECIES_NAME=""
GENOME_FILE=""
THREADS=""
OUTPUT_DIR=""
MIRROR_PATH=""
PROBE_LENGTH=${DEFAULT_PROBE_LENGTH}
SLIDING_WINDOW=${DEFAULT_SLIDING_WINDOW}
CONSISTENCY=${DEFAULT_CONSISTENCY}

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
while getopts ":s:g:t:o:m:l:w:c:h" opt; do
    case $opt in
        s)  SPECIES_NAME="$OPTARG" ;;
        g)  GENOME_FILE="$OPTARG" ;;
        t)  THREADS="$OPTARG" ;;
        o)  OUTPUT_DIR="$OPTARG" ;;
        m)  MIRROR_PATH="$OPTARG" ;;
        l)  PROBE_LENGTH="$OPTARG" ;;
        w)  SLIDING_WINDOW="$OPTARG" ;;
        c)  CONSISTENCY="$OPTARG" ;;
        h)  usage ;;
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate required arguments
# -----------------------------------------------------------------------------
if [[ -z "$SPECIES_NAME" || -z "$GENOME_FILE" || -z "$THREADS" || -z "$OUTPUT_DIR" || -z "$MIRROR_PATH" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

# -----------------------------------------------------------------------------
# Validate parameter values
# -----------------------------------------------------------------------------
if [[ ! -f "$GENOME_FILE" ]]; then
    echo "Error: Genome file '$GENOME_FILE' does not exist." >&2
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
    echo "Error: Threads must be a positive integer." >&2
    exit 1
fi
if [[ "$PROBE_LENGTH" -le 0 ]]; then
    echo "Error: Probe length must be a positive integer." >&2
    exit 1
fi
if [[ "$SLIDING_WINDOW" -le 0 ]]; then
    echo "Error: Sliding window must be a positive integer." >&2
    exit 1
fi
if [[ $(echo "$CONSISTENCY < 0" | bc -l) -eq 1 || $(echo "$CONSISTENCY > 1" | bc -l) -eq 1 ]]; then
    echo "Error: Consistency must be a number between 0 and 1." >&2
    exit 1
fi
if [[ ! -f "$MIRROR_PATH" ]]; then
    echo "Error: Singularity file '$MIRROR_PATH' does not exist." >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Print configuration parameters
# -----------------------------------------------------------------------------
echo "=== Parameters ==="
echo "Species Name       : $SPECIES_NAME"
echo "Genome File        : $GENOME_FILE"
echo "Threads            : $THREADS"
echo "Output Directory   : $OUTPUT_DIR"
echo "Mirror PATH        : $MIRROR_PATH"
echo "Probe Length       : $PROBE_LENGTH"
echo "Sliding Window Size: $SLIDING_WINDOW"
echo "Consistency        : $CONSISTENCY"
echo "==================="

# =============================================================================
# STEP 1: Identify rRNA sequences from the genome
# =============================================================================
echo "[$(date '+%F %T')]: Running probe design tool with the provided configuration..."

# Singularity execution wrapper
readonly SIF_CMD="singularity -q exec -B /media --cleanenv $MIRROR_PATH"

# File paths for intermediate results
readonly RRNA_FA="${OUTPUT_DIR}/${SPECIES_NAME}.barrnap.rrna.fa"
readonly RRNA_GFF="${OUTPUT_DIR}/${SPECIES_NAME}.barrnap.gff"
readonly RRNA_FILE="${OUTPUT_DIR}/${SPECIES_NAME}.rrna.fa"
readonly TMP_FILE="${OUTPUT_DIR}/${SPECIES_NAME}.tmp"
readonly SLIVA_18S="data/fungi_sliva_18S.fasta"
readonly SLIVA_28S="data/fungi_sliva_28S.fasta"

echo "[$(date '+%F %T')]: Identifying rRNA sequences with Barrnap..."

# Decompress genome if it is gzipped
if [[ "$GENOME_FILE" =~ \.gz$ ]]; then
    echo "[$(date '+%F %T')]: Uncompressing genome FASTA..."
    readonly UNCOMPRESSED_FILE="${OUTPUT_DIR}/${SPECIES_NAME}.fa"
    if [[ ! -f "$UNCOMPRESSED_FILE" ]]; then
        gunzip -c "$GENOME_FILE" > "$UNCOMPRESSED_FILE" || {
            echo "[$(date '+%F %T')]: Failed to uncompress genome FASTA file!" >&2
            exit 1
        }
    fi
    GENOME_FILE="$UNCOMPRESSED_FILE"
fi

# Run Barrnap for rRNA annotation
if ! $SIF_CMD barrnap --quiet --kingdom euk --threads "$THREADS" \
    --outseq "$RRNA_FA" --lencutoff 0.8 --reject 0.6 \
    "$GENOME_FILE" > "$RRNA_GFF"; then
    echo "[$(date '+%F %T')]: ERROR: Barrnap failed!" >&2
    rm -f "$RRNA_FA" "$RRNA_GFF"
    exit 2
fi

[[ -s "$RRNA_FA" ]] || { echo "[$(date '+%F %T')]: ERROR: Barrnap generated empty file!" >&2; exit 3; }

# Extract 5S and 5.8S rRNA sequences
$SIF_CMD seqkit grep -r -n -i -p "5S|5_8S" "$RRNA_FA" > "$TMP_FILE"

# Remove duplicate sequences
echo "[$(date '+%F %T')]: Removing duplicate sequences..."
$SIF_CMD seqkit rmdup -s -i -o "$RRNA_FILE" "$TMP_FILE" || {
    echo "[$(date '+%F %T')]: ERROR: Failed to remove duplicates!" >&2
    exit 1
}

# Prepend species name to FASTA headers
sed -i "s/>/>${SPECIES_NAME}_/g" "$RRNA_FILE" || {
    echo "[$(date '+%F %T')]: ERROR: Failed to modify FASTA headers!" >&2
    exit 1
}

# Search Silva database for matching 18S/28S sequences
SPECIES_NAME_MODIFIED="${SPECIES_NAME//_/ }"
echo "[$(date '+%F %T')]: Searching Silva database for '${SPECIES_NAME_MODIFIED}'..."

$SIF_CMD seqkit grep -r -n -i -p "${SPECIES_NAME_MODIFIED}" "$SLIVA_18S" >> "$RRNA_FILE" || {
    echo "[$(date '+%F %T')]: WARNING: Failed to find 18S sequences in Silva database." >&2
}
$SIF_CMD seqkit grep -r -n -i -p "${SPECIES_NAME_MODIFIED}" "$SLIVA_28S" >> "$RRNA_FILE" || {
    echo "[$(date '+%F %T')]: WARNING: Failed to find 28S sequences in Silva database." >&2
}

# Clean up intermediate files
rm -f "$RRNA_FA" "$RRNA_GFF" "$TMP_FILE"
if [[ -f "$UNCOMPRESSED_FILE" ]]; then
    rm -f "$UNCOMPRESSED_FILE" "${UNCOMPRESSED_FILE}.fai"
fi

# =============================================================================
# STEP 2: Generate rRNA probes
# =============================================================================
echo "[$(date '+%F %T')]: Generating rRNA probes..."

# Cluster rRNA sequences with CD-HIT-EST
$SIF_CMD cd-hit-est -i "${RRNA_FILE}" -c "$CONSISTENCY" -T "$THREADS" -ap 1 -M 0 -d 0 -o "${RRNA_FILE}.merge"

# Generate sliding-window probes
readonly SLIDING_FILE="${OUTPUT_DIR}/${SPECIES_NAME}_sliding.fa"
$SIF_CMD seqkit sliding -g -s "$SLIDING_WINDOW" -W "$PROBE_LENGTH" -j "$THREADS" -o "$SLIDING_FILE" "${RRNA_FILE}.merge"
rm -f "${RRNA_FILE}.merge.clstr"

# =============================================================================
# STEP 3: Cluster sliding probes with CD-HIT-EST
# =============================================================================
echo "[$(date '+%F %T')]: Merging rRNA probes..."
readonly MERGE_FILE="${OUTPUT_DIR}/${SPECIES_NAME}_cdhit_merge.fa"
$SIF_CMD cd-hit-est -i "$SLIDING_FILE" -c "$CONSISTENCY" -T "$THREADS" -ap 1 -M 0 -d 0 -n 4 -o "$MERGE_FILE"

# =============================================================================
# STEP 4: Remove potential dimer probes
# =============================================================================
echo "[$(date '+%F %T')]: Removing potential dimer probes..."
readonly REV_DIMER="src/remove_dimer_para.py"
readonly FINAL_FILE="${OUTPUT_DIR}/${SPECIES_NAME}_rrna_probes.fa"

$SIF_CMD python "$REV_DIMER" "$MERGE_FILE" "$FINAL_FILE" -15 "$THREADS" 0.1

# Clean up temporary files
rm -f "$SLIDING_FILE" "${MERGE_FILE}.clstr"

echo "[$(date '+%F %T')]: Completed! Final probes written to: $FINAL_FILE"
