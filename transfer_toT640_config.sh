#!/bin/bash
#
# Description:
#   Rewrites paths in the library configuration TSV for use on the
#   remote host (T640).  For each specified line, it extracts the
#   prefix, strips the trailing strain number (e.g., -1, -2, -3),
#   and rebuilds the genome and GTF paths under a remote reference
#   base directory.  The result is written to a new output file.
#
# Usage:
#   bash transfer_toT640_config.sh
#
# Configuration (edit before running):
#   INPUT_FILE       - Path to the source library configuration TSV.
#   OUTPUT_FILE      - Path where the rewritten TSV is saved.
#   TARGET_LINES     - Array of 1-based line numbers to convert.
#   REMOTE_REF_BASE  - Remote base path for reference files.

# --- Configuration ---
INPUT_FILE="example/cand_lib_config.tsv"
OUTPUT_FILE="example/cand_lib_config_toT640.tsv"
# Line numbers to process.
TARGET_LINES=(2 8 13 14 15)

# Remote base path for reference files.
REMOTE_REF_BASE="/media/data4/hxy/PanFungi/0_reference"

echo "Extracting and converting paths..."

# Truncate or create the output file.
> "$OUTPUT_FILE"

for LINE_NUM in "${TARGET_LINES[@]}"; do
    # Fetch the line content from the input file.
    LINE_CONTENT=$(sed -n "${LINE_NUM}p" "$INPUT_FILE")

    if [[ -z "$LINE_CONTENT" ]]; then
        echo "Warning: Line $LINE_NUM does not exist, skipping."
        continue
    fi

    # Use awk to rewrite genome and GTF paths.
    # Field layout: $1=species, $2=prefix, $5=genome, $6=gtf.
    # Logic: extract strain from prefix, derive basename, and
    #        prepend the remote reference base path.
    echo "$LINE_CONTENT" | awk -v ref_base="$REMOTE_REF_BASE" 'BEGIN {OFS="\t"} {
        # 1. Parse strain (strip trailing -1, -2, or -3 from the prefix).
        strain = $2;
        sub(/-[123]$/, "", strain);

        # 2. Get the filename (basename) from each path.
        n = split($5, a, "/"); genome_file = a[n];
        m = split($6, b, "/"); gtf_file = b[m];

        # 3. Construct new paths under the remote reference base.
        $5 = ref_base "/" strain "/" genome_file;
        $6 = ref_base "/" strain "/" gtf_file;

        print $0;
    }' >> "$OUTPUT_FILE"
done

echo "Conversion complete! New file saved to: $OUTPUT_FILE"
echo "------------------------------------------"
cat "$OUTPUT_FILE"
