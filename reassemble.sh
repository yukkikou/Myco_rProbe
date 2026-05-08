#!/bin/bash
# Reassemble pan_fungi.sif from split parts
# Usage: bash reassemble.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PARTS_DIR="$SCRIPT_DIR/singularity"
OUTPUT="$PARTS_DIR/pan_fungi.sif"

if [ -f "$OUTPUT" ]; then
    echo "Error: $OUTPUT already exists. Remove it first if you want to reassemble." >&2
    exit 1
fi

cat "$PARTS_DIR"/pan_fungi.sif.part* > "$OUTPUT"

if [ $? -eq 0 ]; then
    echo "Successfully reassembled: $OUTPUT"
    ls -lh "$OUTPUT"
else
    echo "Error: Reassembly failed." >&2
    exit 1
fi
