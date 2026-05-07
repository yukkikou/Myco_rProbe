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
#   SLURM array job script for running BLAST on circRNA FASTA files.
#   Each array task reads its FASTA file path from FASTA_LIST
#   (one path per line, indexed by SLURM_ARRAY_TASK_ID), runs
#   blastn inside a Singularity container, and writes results
#   to RESULTS_DIR.
#
# Usage:
#   Submitted as a SLURM array job.  The following environment
#   variables must be set (typically exported before sbatch):
#     SLURM_ARRAY_TASK_ID  - Array task index (1-based).
#     SLURM_CPUS_PER_TASK  - CPUs allocated to this task.
#     FASTA_LIST           - File listing FASTA paths, one per line.
#     RESULTS_DIR          - Output directory for BLAST results.
#     BLAST_DB             - Path to the BLAST nucleotide database.
#     BLAST_SIF            - Path to the BLAST Singularity image.

# Get the FASTA file for this array task.
FASTA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FASTA_LIST")

if [ -z "$FASTA_FILE" ] || [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: Invalid or missing fasta file for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract prefix from the filename (strip .fasta extension).
BASENAME=$(basename "$FASTA_FILE" .fasta)
OUTPUT_FILE="$RESULTS_DIR/${BASENAME}.blast"

echo "[$(date '+%F %T')]: Starting BLAST for $BASENAME"
echo "Input: $FASTA_FILE"
echo "Output: $OUTPUT_FILE"
echo "Database: $BLAST_DB"

#######################################
# Run BLAST against the nucleotide database.
#######################################
singularity exec -B /media:/media "$BLAST_SIF" blastn \
    -query "$FASTA_FILE" \
    -db "$BLAST_DB" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-5 \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -max_target_seqs 100 \
    -out "$OUTPUT_FILE"

if [ $? -eq 0 ]; then
    echo "[$(date '+%F %T')]: BLAST completed successfully for $BASENAME"
    echo "Results saved to: $OUTPUT_FILE"
else
    echo "[$(date '+%F %T')]: ERROR: BLAST failed for $BASENAME"
    exit 1
fi
