#!/bin/bash
#
# Copyright 2025 The Myco_rProbe Authors.
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
# Description:
#   Main workflow driver for circular RNA (circRNA) probe design.
#   Reads a genome-to-GTF mapping file, validates inputs, and extracts
#   circRNA sequences using bedtools getfasta inside a Singularity
#   container.  Each input line produces one FASTA output file.
#
# Usage:
#   bash myco_rprobe_circseq.sh -g <genome_mapping_list> -t <threads> \
#                                -o <output_dir> -m <mirror_path>

set -euo pipefail

#######################################
# Display usage information and exit.
# Globals:
#   $0 - script name
# Arguments:
#   None
# Returns:
#   Exits with status 1.
#######################################
usage() {
    echo "Usage: bash $0 -g <genome_mapping_list> -r <rrna_dir> -t <threads> -o <output_dir> -m <mirror_path> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -g <genome_mapping_list>    Genome list file (e.g., 'example/genome_mapping.tsv')."
    echo "  -t <threads>                Number of threads (e.g., 4)."
    echo "  -o <output_dir>             Path to the output directory (e.g., 'mapping')."
    echo "  -m <mirror_path>            Singularity file path provided (bind mirror)."
    echo ""
    echo "Optional arguments:"
    echo "  -h                      Show this help message and exit."
    echo ""
    exit 1
}

# Initialize variables.
GENOME_LIST=""
THREADS=""
OUTPUT_DIR=""
MIRROR_PATH=""

# Parse command-line arguments.
while getopts ":r:g:t:o:m:h" opt; do
    case $opt in
        g)  # Genome file path.
            GENOME_LIST="$OPTARG"
            ;;
        t)  # Number of threads.
            THREADS="$OPTARG"
            ;;
        o)  # Output directory path.
            OUTPUT_DIR="$OPTARG"
            ;;
        m)  # Singularity mirror image path.
            MIRROR_PATH="$OPTARG"
            ;;
        h)  # Display help and exit.
            usage
            ;;
        \?) # Invalid option.
            echo "Error: Invalid option -$OPTARG" >&2
            usage
            ;;
        :)  # Missing required argument.
            echo "Error: Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

#######################################
# Validate that all required arguments are provided.
#######################################
if [[ -z "$GENOME_LIST" || -z "$THREADS" || -z "$OUTPUT_DIR" || -z "$MIRROR_PATH" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

# Validate that the genome mapping list exists.
if [[ ! -f "$GENOME_LIST" ]]; then
    echo "Error: Genome file '$GENOME_LIST' does not exist." >&2
    exit 1
fi

# Validate that threads is a positive integer.
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
    echo "Error: Threads must be a positive integer." >&2
    exit 1
fi

# Validate that the Singularity image exists.
if [[ ! -f "$MIRROR_PATH" ]]; then
    echo "Error: Singularity file '$MIRROR_PATH' does not exist." >&2
    exit 1
fi

# Mark validated parameters as read-only.
readonly GENOME_LIST THREADS OUTPUT_DIR MIRROR_PATH

#######################################
# Print the configuration parameters for debugging.
#######################################
echo "=== Parameters ==="
echo "Genome Mapping     : $GENOME_LIST"
echo "Threads            : $THREADS"
echo "Output Directory   : $OUTPUT_DIR"
echo "Mirror PATH        : $MIRROR_PATH"
echo "==================="

#######################################
# Core workflow: extract circRNA sequences.
#######################################
echo "[$(date '+%F %T')]: Running probe design tool with the provided configuration..."
echo "Command: ./myco_rprobe_circseq.sh \
    --genome_mapping_list '$GENOME_LIST' --threads '$THREADS' \
    --output '$OUTPUT_DIR' --mirror '$MIRROR_PATH'"

mkdir -p "$OUTPUT_DIR"

echo "[$(date '+%F %T')]: Extract targeted circRNA sequences ..."

# Build the Singularity command prefix as an array to avoid word-splitting issues.
SIF_CMD=(singularity -q exec -B /media --cleanenv "$MIRROR_PATH")

# Read the genome mapping list line by line.
# Expected format: tab-separated columns (genome_path, gtf_path).
while IFS=$'\t' read -r genome gtf; do
    # Skip empty lines and comment lines.
    [[ -z "$genome" || "$genome" =~ ^[[:space:]]*# ]] && continue

    if [[ ! -f "$genome" || ! -f "$gtf" ]]; then
        echo "[$(date '+%F %T')]: ERROR: Missing required files. Skipping."
        [[ ! -f "$genome" ]] && echo "  - Genome not found at: $genome"
        [[ ! -f "$gtf" ]] && echo "  - GTF not found at: $gtf"
        continue
    fi

    gtf_basename=$(basename "$gtf")
    prefix="${gtf_basename%.*}"
    circ_seq="$OUTPUT_DIR/${prefix}.circfl.fa"

    "${SIF_CMD[@]}" bedtools getfasta -name -fi "$genome" -bed "$gtf" -fo "$circ_seq"

done < <(sed 's/\r$//' "$GENOME_LIST")

echo "[$(date '+%F %T')]: Extract Finished"
