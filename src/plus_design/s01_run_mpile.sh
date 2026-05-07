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
#   Runs samtools mpileup on a 28S BAM file to generate a pileup
#   for probe design.  Expects a prefix argument that is used to
#   derive the BAM, reference FASTA, and output file names.
#
# Usage:
#   bash s01_run_mpile.sh <prefix>
#   Example: bash s01_run_mpile.sh Mucor

prefix=$1

# --- Input parameters ---
BAM_FILE="${prefix}.28S.bam"
REF_FASTA="${prefix}_28S_sliva_28S.fa"
# Extract the chromosome name from the reference FASTA header.
CHROMOSOME=$(grep "^>" "$REF_FASTA" | head -1 | cut -d " " -f 1 | sed 's/>//g')
MIN_BASE_QUALITY="20"
OUTPUT_PILEUP="${prefix}_28s.pileup"

# Get sequence length from the reference FASTA.
#SEQUENCE_LENGTH=3639

SEQUENCE_LENGTH=$(singularity exec -B /media:/media /media/share/data4/container/seqkit/seqkit-2.2.0.sif seqkit stat -T Mucor_28S_sliva_28S.fa  | awk 'NR>1 {print $8}')

# --- Pre-flight checks ---

if [ ! -f "$BAM_FILE" ]; then
  echo "Error: BAM file '$BAM_FILE' not found."
  exit 1
fi

if [ ! -f "$BAM_FILE".fai ]; then
  echo "Building index"
  singularity exec -B /media:/media /media/share/data4/container/samtools/samtools-1.14.sif samtools index "$BAM_FILE"
fi

if [ ! -f "$REF_FASTA" ]; then
  echo "Error: Reference FASTA '$REF_FASTA' not found."
  exit 1
fi

if [ -z "$CHROMOSOME" ]; then
  echo "Error: Could not determine chromosome name from FASTA file."
  exit 1
fi

if [ -z "$SEQUENCE_LENGTH" ]; then
  echo "Error: Could not determine sequence length from FASTA file."
  exit 1
fi

#######################################
# Build and execute the samtools mpileup command.
#######################################
MPILEUP_COMMAND="singularity exec -B /media:/media /media/share/data4/container/samtools/samtools-1.14.sif samtools mpileup -aa -q $MIN_BASE_QUALITY -f $REF_FASTA -r $CHROMOSOME:1-$SEQUENCE_LENGTH -o $OUTPUT_PILEUP $BAM_FILE"

echo "Running: $MPILEUP_COMMAND"

# Execute the mpileup command.
eval "$MPILEUP_COMMAND"

# Check if the command was successful.
if [ $? -eq 0 ]; then
  echo "samtools mpileup completed successfully. Output is in '$OUTPUT_PILEUP'."
else
  echo "samtools mpileup failed. Check the error messages above."
  exit 1
fi

exit 0
