#!/bin/bash
#
# Description:
#   Transfers reference genome, GTF annotation, and clean data
#   files for specified samples to a remote host (T640).  Reads a
#   configuration file and processes only the line numbers listed
#   in TARGET_LINES.  Deduplication is applied so that the same
#   strain's reference files are transferred only once.
#
# Usage:
#   bash transfer_toT640.sh
#
# Configuration (edit before running):
#   CONFIG_FILE  - Path to the library configuration TSV.
#   TARGET_HOST  - Remote host (SSH alias or hostname).
#   TARGET_LINES - Array of 1-based line numbers to process.

# --- Configuration ---
CONFIG_FILE="example/cand_lib_config.tsv"
TARGET_HOST="toT640"
# Line numbers to process, space-separated.
TARGET_LINES=(2 8 13 14 15)

# Associative array tracking synced strains (deduplication).
declare -A SYNCED_STRAINS

echo "Processing specified lines: ${TARGET_LINES[*]}"

for LINE_NUM in "${TARGET_LINES[@]}"; do
    # Fetch content of the specified line from the config file.
    LINE_CONTENT=$(sed -n "${LINE_NUM}p" "$CONFIG_FILE")

    # Skip if the line is empty.
    if [[ -z "$LINE_CONTENT" ]]; then
        echo "Warning: Line $LINE_NUM does not exist, skipping."
        continue
    fi

    # Parse columns (1:species, 2:prefix, 5:genome, 6:gtf).
    SPECIES=$(echo "$LINE_CONTENT" | awk '{print $1}')
    PREFIX=$(echo "$LINE_CONTENT" | awk '{print $2}')
    GENOME=$(echo "$LINE_CONTENT" | awk '{print $5}')
    GTF=$(echo "$LINE_CONTENT" | awk '{print $6}')

    # Parse strain: strip trailing -1, -2, or -3 from the prefix.
    STRAIN=$(echo "$PREFIX" | sed 's/-[123]$//')

    echo "------------------------------------------"
    echo "Processing line $LINE_NUM: Prefix=$PREFIX -> Strain=$STRAIN"

    # --- Step 1: Transfer genome and GTF (deduplicated) ---
    if [[ -z "${SYNCED_STRAINS[$STRAIN]}" ]]; then
        REMOTE_REF_DIR="/media/data4/hxy/PanFungi/0_reference/${STRAIN}"

        echo "Syncing reference genome to: $REMOTE_REF_DIR"
        ssh "$TARGET_HOST" "mkdir -p $REMOTE_REF_DIR"
        scp "$GENOME" "$GTF" "${TARGET_HOST}:${REMOTE_REF_DIR}/"

        SYNCED_STRAINS[$STRAIN]=1
    else
        echo "Reference files for strain $STRAIN already synced, skipping."
    fi

    # --- Step 2: Transfer clean data directory (corrected path structure) ---
    LOCAL_CLEAN_DIR="/media/share/data5/1610305236/Panfungi/1_data/Candida/${PREFIX}/0_clean"
    # Remote target is set to the prefix level.
    REMOTE_PREFIX_DIR="/media/data4/hxy/PanFungi/1_data/Candida/${PREFIX}"

    if [ -d "$LOCAL_CLEAN_DIR" ]; then
        echo "Syncing clean data to: $REMOTE_PREFIX_DIR/0_clean"
        # Create the prefix directory on the remote host first.
        ssh "$TARGET_HOST" "mkdir -p $REMOTE_PREFIX_DIR"
        # Copy the entire local 0_clean directory to the remote prefix directory.
        scp -r "$LOCAL_CLEAN_DIR" "${TARGET_HOST}:${REMOTE_PREFIX_DIR}/"
    else
        echo "Error: Local directory does not exist $LOCAL_CLEAN_DIR"
    fi

done

echo "------------------------------------------"
echo "All tasks completed."
