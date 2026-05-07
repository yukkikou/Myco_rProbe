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
# Myco_rProbe - Multi-species rRNA probe design pipeline
#
# This script handles multi-species rRNA probe design. It processes a list
# of fungal genomes, identifies rRNA sequences via Barrnap, generates
# sliding-window probes across multiple rRNA types (5S, 5.8S, 18S, 28S),
# clusters probes with CD-HIT-EST, evaluates coverage via BLAST, and groups
# probes by phylogenetic levels (core vs. plus).
#
# Usage:
#   bash myco_rprobe_multi.sh -g <genome_list> -t <threads> \
#                              -o <output_dir> -m <mirror_path> [options]
#
# Author: Myco_rProbe Team
# Date:   2024
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Constants and default values
# -----------------------------------------------------------------------------
readonly DEFAULT_RRNA_CONF="data/rRNA_list.tsv"

# -----------------------------------------------------------------------------
# Display help information
#
# Prints usage guide with required and optional arguments.
#
# Globals:
#   $0 - script name
# Returns:
#   None (exits with code 1)
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: bash $0 -g <genome_list> -t <threads> -o <output_dir> -m <mirror_path> [options]

Required arguments:
  -g <genome_list>       Path to the genome list file (TSV: genome_path<TAB>genus<TAB>species).
                         Example: 'example/fungi_genome.list'
  -t <threads>           Number of CPU threads (e.g., 4).
  -o <output_dir>        Output directory path (e.g., 'multi').
  -m <mirror_path>       Path to the Singularity image file (binds /media).

Optional arguments:
  -r <rRNA_config>       rRNA configuration file (TSV: rna_type sliding_step consistency width).
                         Default: 'data/rRNA_list.tsv'.
  -h                     Show this help message and exit.

Examples:
  bash myco_rprobe_multi.sh -g example/fungi_genome.list -t 10 -o pathgen_fungi -m /path/to/image.sif
  bash myco_rprobe_multi.sh -g example/fungi_genome.list -t 20 -o output -m /path/to/image.sif -r data/rRNA_list2.tsv
EOF
    exit 1
}

# -----------------------------------------------------------------------------
# Initialize variables
# -----------------------------------------------------------------------------
RRNA_CONF="${DEFAULT_RRNA_CONF}"
GENOME_LIST=""
THREADS=""
OUTPUT_DIR=""
MIRROR_PATH=""

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
while getopts ":r:g:t:o:m:h" opt; do
    case $opt in
        r)  RRNA_CONF="$OPTARG" ;;
        g)  GENOME_LIST="$OPTARG" ;;
        t)  THREADS="$OPTARG" ;;
        o)  OUTPUT_DIR="$OPTARG" ;;
        m)  MIRROR_PATH="$OPTARG" ;;
        h)  usage ;;
        \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate arguments
# -----------------------------------------------------------------------------
if [[ -z "$GENOME_LIST" || -z "$THREADS" || -z "$OUTPUT_DIR" || -z "$MIRROR_PATH" ]]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

if [[ ! -f "$RRNA_CONF" ]]; then
    echo "Error: rRNA config file '$RRNA_CONF' does not exist." >&2
    exit 1
fi
if [[ ! -f "$GENOME_LIST" ]]; then
    echo "Error: Genome list file '$GENOME_LIST' does not exist." >&2
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

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Print configuration parameters
# -----------------------------------------------------------------------------
echo "=== Parameters ==="
echo "rRNA config        : $RRNA_CONF"
echo "Genome list        : $GENOME_LIST"
echo "Threads            : $THREADS"
echo "Output Directory   : $OUTPUT_DIR"
echo "Mirror PATH        : $MIRROR_PATH"
echo "==================="

# =============================================================================
# File path definitions
# =============================================================================
readonly SIF_CMD="singularity -q exec -B /media --cleanenv $MIRROR_PATH"
readonly SLIVA_18S="data/fungi_sliva_18S.fasta"
readonly SLIVA_28S="data/fungi_sliva_28S.fasta"
readonly ALL_RRNA="${OUTPUT_DIR}/all_rrna.fasta"
readonly SEP_5S="${OUTPUT_DIR}/sep_5s.fasta"
readonly SEP_5_8S="${OUTPUT_DIR}/sep_5_8s.fasta"
readonly SEP_18S="${OUTPUT_DIR}/sep_18s.fasta"
readonly SEP_28S="${OUTPUT_DIR}/sep_28s.fasta"
readonly REV_PROBE="src/remove_inner.py"
readonly COV_CAL="src/parse_coverage.py"
readonly GROUP_PROBES="src/group_probes.py"
readonly GROUP_PROBES_PLUS="src/group_probes_plus.py"
readonly GROUPED_CAL="src/grouped_coverage.py"
readonly GROUPED_CAL_PLUS="src/grouped_coverage_plus.py"
readonly DEC_TREE="src/decision_tree.py"
readonly FUNGI_TAX="data/fungi_taxonomy_info.tsv"
readonly FUNGI_NODE="data/fungi_taxonomy_nodes.tsv"

echo "[$(date '+%F %T')]: Starting probe design pipeline..."

# =============================================================================
# STEP 1: Process individual strains
# =============================================================================

# -----------------------------------------------------------------------------
# Process a single strain: run Barrnap for rRNA identification
#
# Args:
#   $1 - strain_name: species identifier
#   $2 - genome_file: path to genome FASTA (may be .gz)
#   $3 - outdir: output directory
# Returns:
#   0 on success, 1 if gunzip fails, 2 if Barrnap fails
# -----------------------------------------------------------------------------
process_strain() {
    local strain_name="$1"
    local genome_file="$2"
    local outdir="$3"
    local rrna_fa="$outdir/${strain_name}.barrnap.rrna.fa"
    local rrna_gff="$outdir/${strain_name}.barrnap.gff"
    local uncompressed_file="$outdir/${strain_name}.fasta"

    if [[ "$genome_file" == *.gz ]]; then
        echo "[$(date '+%F %T')]: Unzipping genome for strain $strain_name..."
        if [[ ! -f "$uncompressed_file" ]]; then
            gunzip -c "$genome_file" > "$uncompressed_file" || {
                echo "[$(date '+%F %T')]: ERROR: Failed to unzip genome for $strain_name." >&2
                return 1
            }
        fi
        genome_file="$uncompressed_file"
    fi

    echo "[$(date '+%F %T')]: Running Barrnap for strain $strain_name..."

    if ! $SIF_CMD barrnap --kingdom euk --threads "$THREADS" \
            --outseq "$rrna_fa" --lencutoff 0.8 --reject 0.6 "$genome_file" > "$rrna_gff"; then
        echo "[$(date '+%F %T')]: ERROR: Barrnap failed for strain $strain_name." >&2
        rm -f "$rrna_fa" "$rrna_gff"
        return 2
    fi

    if [[ -f "$uncompressed_file" ]]; then
        rm -f "$uncompressed_file" "${uncompressed_file}.fai"
    fi

    echo "[$(date '+%F %T')]: Barrnap completed for strain $strain_name." && rm -f "$rrna_gff"
    return 0
}

# -----------------------------------------------------------------------------
# Deduplicate rRNA sequences for a strain
#
# Removes sequence duplicates, extracts 5S/5.8S rRNA, and appends
# matching Silva 18S/28S sequences from the reference database.
#
# Args:
#   $1 - strain_name: species identifier
#   $2 - rrna_seq: path to Barrnap output FASTA
#   $3 - outdir: output directory
# Returns:
#   0 on success, 1 if seqkit dedup fails, 2 if seqkit grep fails,
#   3 if output file is empty
# -----------------------------------------------------------------------------
rrna_uniq() {
    local strain_name="$1"
    local rrna_seq="$2"
    local outdir="$3"
    local tmp_file="$outdir/${strain_name}.rrna.tmp"
    local output_file="$outdir/${strain_name}.rrna.rmdup.fa"

    echo "[$(date '+%F %T')]: Deduplicating rRNA for strain $strain_name..."

    if [[ ! -f "$output_file" ]]; then
        if ! $SIF_CMD seqkit rmdup -s -i -o "$tmp_file" "$rrna_seq"; then
            echo "[$(date '+%F %T')]: ERROR: seqkit rmdup failed for $strain_name." >&2
            return 1
        fi

        sed -i "s/>/>${strain_name}_/g" "$tmp_file"

        if ! $SIF_CMD seqkit grep -r -n -i -p "5S|5_8S" "$tmp_file" > "$output_file"; then
            echo "[$(date '+%F %T')]: ERROR: seqkit grep failed for $strain_name." >&2
            return 2
        fi

        # Search Silva database for matching sequences
        local first_part
        first_part=$(echo "$strain_name" | cut -d'_' -f1)
        local second_part
        second_part=$(echo "$strain_name" | cut -d'_' -f2 | cut -d'.' -f1)
        local species_name="${first_part} ${second_part}"

        echo "[$(date '+%F %T')]: Searching Silva database for $species_name..."
        $SIF_CMD seqkit grep -r -n -i -p "${species_name}" "$SLIVA_18S" >> "$output_file" 2>/dev/null || true
        $SIF_CMD seqkit grep -r -n -i -p "${species_name}" "$SLIVA_28S" >> "$output_file" 2>/dev/null || true
        rm -f "$rrna_seq" "$tmp_file"
    fi

    [[ -s "$output_file" ]] || {
        echo "[$(date '+%F %T')]: ERROR: Deduplicated rRNA file is empty for $strain_name." >&2
        return 3
    }

    echo "[$(date '+%F %T')]: Deduplication completed for strain $strain_name."
    return 0
}

# Process all strains from the genome list
while IFS=$'\t' read -r genome genus species; do
    if [[ -z "$genome" || -z "$genus" || -z "$species" ]]; then
        echo "[$(date '+%F %T')]: ERROR: Insufficient fields in genome list entry" >&2
        continue
    fi
    genome="$(realpath -e "$genome")" || {
        echo "[$(date '+%F %T')]: ERROR: Genome file not found: $genome." >&2
        continue
    }
    echo "[$(date '+%F %T')]: Processing strain: Genome=$genome Genus=$genus Species=$species..."
    if process_strain "$species" "$genome" "$OUTPUT_DIR"; then
        if rRNA_uniq "$species" "$OUTPUT_DIR/${species}.barrnap.rrna.fa" "$OUTPUT_DIR"; then
            echo "[$(date '+%F %T')]: Pipeline completed for strain $species."
        else
            echo "[$(date '+%F %T')]: ERROR: rRNA deduplication failed for $species." >&2
        fi
    else
        echo "[$(date '+%F %T')]: ERROR: Failed to process strain $species." >&2
    fi
done < "$GENOME_LIST"
echo "[$(date '+%F %T')]: All strains processing completed."

# =============================================================================
# STEP 2: Generate and evaluate probes for each rRNA type
# =============================================================================
echo "[$(date '+%F %T')] #############################################################"
echo "[$(date '+%F %T')]: Generating rRNA probes..."

# Concatenate all strain rRNA sequences
cat "$OUTPUT_DIR"/*.rrna.rmdup.fa > "$ALL_RRNA"
rm -f "$OUTPUT_DIR"/*.rrna.rmdup.fa

echo "[$(date '+%F %T')]: Separating rRNA by type..."
declare -A RNA_FILES=(
    ["5S"]="${SEP_5S}"
    ["5_8S"]="${SEP_5_8S}"
    ["18S"]="${SEP_18S}"
    ["28S"]="${SEP_28S}"
)

# Separate by rRNA type and cluster with CD-HIT-EST
for rna_type in "${!RNA_FILES[@]}"; do
    rna_file="${RNA_FILES[$rna_type]}"

    echo "[$(date '+%F %T')]: Separating ${rna_type} -> ${rna_file}..."
    $SIF_CMD seqkit grep -r -n -i -p "$rna_type" "${ALL_RRNA}" > "${rna_file}"
    $SIF_CMD cd-hit-est -i "${rna_file}" -c 0.8 -T "$THREADS" -ap 1 -M 0 -d 0 -n 4 -o "${rna_file}.merge"

    rm -f "${rna_file}.merge.clstr"
    original_count=$(grep "^>" "${rna_file}" | wc -l)
    merged_count=$(grep "^>" "${rna_file}.merge" | wc -l)
    echo "[$(date '+%F %T')]: ${rna_type}: ${original_count} seqs merged to ${merged_count}"
done

# Generate sliding probes and evaluate for each configuration
while IFS=$'\t' read -r rna_type sliding_step constance width; do
    # Map rRNA type to its FASTA file
    declare -A RNA_IDENTIFY=(
        ["5S"]="${SEP_5S}"
        ["5_8S"]="${SEP_5_8S}"
        ["18S"]="${SEP_18S}"
        ["28S"]="${SEP_28S}"
    )

    rna_file="${RNA_IDENTIFY[$rna_type]}"
    if [[ -z "$rna_file" ]]; then
        echo "Error: Unsupported rRNA type: $rna_type" >&2
        exit 1
    fi

    echo "[$(date '+%F %T')]: Processing $rna_type (step=$sliding_step, identity=$constance, width=$width)"

    out_prefix="${OUTPUT_DIR}/${rna_type}_${sliding_step}_${constance}_${width}"
    mkdir -p "${out_prefix}"

    sliding_file="${out_prefix}/sliding.fasta"
    merge_file="${out_prefix}/cdhit_merge.fasta"
    final_file="${out_prefix}/cdhit_final.fasta"
    align_res="${out_prefix}/cdhit_final.alig.out"

    # Generate sliding-window sequences
    echo "[$(date '+%F %T')]: Generating sliding probes..."
    $SIF_CMD seqkit sliding -g -s "$sliding_step" -W "$width" -j "$THREADS" -o "${sliding_file}" "${rna_file}.merge"

    # Cluster sliding probes
    echo "[$(date '+%F %T')]: Merging probes with CD-HIT-EST..."
    $SIF_CMD cd-hit-est -i "$sliding_file" -c "$constance" -T "$THREADS" -ap 1 -M 0 -d 0 -n 4 -o "${merge_file}"

    # Remove internal overlapping probes
    echo "[$(date '+%F %T')]: Removing internal probes..."
    $SIF_CMD python "$REV_PROBE" "${merge_file}" "${final_file}"

    # Build BLAST database if needed
    db_files=(
        "${OUTPUT_DIR}/db/${rna_type}/${rna_type}.nhr"
        "${OUTPUT_DIR}/db/${rna_type}/${rna_type}.nin"
        "${OUTPUT_DIR}/db/${rna_type}/${rna_type}.nsq"
    )
    db_exists=true
    for db_file in "${db_files[@]}"; do
        if [[ ! -f "$db_file" ]]; then
            db_exists=false
            break
        fi
    done

    if ! $db_exists; then
        echo "[$(date '+%F %T')]: Building BLAST database for ${rna_type}..."
        mkdir -p "${OUTPUT_DIR}/db/${rna_type}"
        $SIF_CMD makeblastdb -in "$rna_file" -dbtype nucl -out "${OUTPUT_DIR}/db/${rna_type}/${rna_type}"
    else
        echo "[$(date '+%F %T')]: BLAST database for ${rna_type} already exists. Skipping."
    fi

    # Evaluate probe coverage via BLAST
    echo "[$(date '+%F %T')]: BLASTing probes against rRNA references..."
    $SIF_CMD blastn -query "${final_file}" \
        -db "${OUTPUT_DIR}/db/${rna_type}/${rna_type}" \
        -evalue 1e-6 \
        -outfmt 6 \
        -out "${align_res}" \
        -num_threads "$THREADS"

    echo "[$(date '+%F %T')]: Calculating coverage..."
    $SIF_CMD python "${COV_CAL}" --blast "${align_res}" \
        --species_18S "${SLIVA_18S}" --species_28S "${SLIVA_28S}" \
        --fasta "${rna_file}" \
        --outputdir "${out_prefix}"

    # Clean intermediate files
    rm -f "${sliding_file}" "${merge_file}" "${merge_file}.clstr"

done < "$RRNA_CONF"

rm -f "$ALL_RRNA"

# =============================================================================
# STEP 3: Summarize probe results and select best configurations
# =============================================================================
echo "[$(date '+%F %T')]: Summarizing probe results..."
readonly COV_STAT="${OUTPUT_DIR}/coverage_summary_results.txt"
echo -e "full_class_name\thigh_class_name\tprobe_num\tq1\tq2\tq3\tscore" > "$COV_STAT"

declare -A MAX_SCORES
declare -A BEST_DIRS

# -----------------------------------------------------------------------------
# Extract high-level class name from directory name
#
# Args:
#   $1 - dir_basename: directory name (e.g., "5S_40_0.8_40")
# Sets globals:
#   full_class_name, high_class_name
# -----------------------------------------------------------------------------
extract_classes() {
    local dir_basename="$1"
    full_class_name="$dir_basename"
    high_class_name=$(echo "$full_class_name" | cut -d'_' -f1)
}

# -----------------------------------------------------------------------------
# Calculate statistics for a probe set
#
# Computes probe count and coverage quartiles (Q1, Q2, Q3) from the
# coverage summary file, then derives a composite score.
#
# Args:
#   $1 - fasta_file: path to cdhit_final.fasta
#   $2 - coverage_file: path to rRNA_coverage_summary.txt
# Sets globals:
#   probe_num, q1, q2, q3, score
# -----------------------------------------------------------------------------
calculate_stats() {
    local fasta_file="$1"
    local coverage_file="$2"

    probe_num="NA"; q1="NA"; q2="NA"; q3="NA"; score="NA"

    if [[ -f "$fasta_file" ]]; then
        local line_count
        line_count=$(wc -l < "$fasta_file")
        probe_num=$((line_count / 2))
    fi

    if [[ -f "$coverage_file" && $(awk '{print NF}' "$coverage_file" | head -1) -ge 3 ]]; then
        local ratio
        ratio=$(awk 'NR > 1 {print $3/$2*100}' "$coverage_file")
        if [[ -n $ratio ]]; then
            local csvtk_res
            csvtk_res=$(echo "$ratio" | $SIF_CMD csvtk summary -t -f 1:q1,1:q2,1:q3 | tail -1)
            q1=$(echo "$csvtk_res" | cut -f1)
            q2=$(echo "$csvtk_res" | cut -f2)
            q3=$(echo "$csvtk_res" | cut -f3)
        fi
    fi

    if [[ "$probe_num" != "NA" && "$q3" != "NA" ]]; then
        score=$(awk -v pn="$probe_num" -v q3="$q3" 'BEGIN {print (1 / pn) * q3}')
    fi
}

# -----------------------------------------------------------------------------
# Update max scores tracking
#
# Args:
#   $1 - class: high-level class name
#   $2 - score: computed score
#   $3 - dir: directory path
#   $4 - full_class: full class name
# -----------------------------------------------------------------------------
update_max_scores() {
    local class="$1"
    local score="$2"
    local dir="$3"
    local full_class="$4"

    if [[ "$score" != "NA" ]]; then
        if [[ -z ${MAX_SCORES["$class"]} || $(awk -v s1="${MAX_SCORES[$class]}" -v s2="$score" \
            'BEGIN {print (s2 > s1)}') == 1 ]]; then
            MAX_SCORES["$class"]="$score"
            BEST_DIRS["$class"]="$dir#$full_class"
        elif [[ "${MAX_SCORES[$class]}" == "$score" ]]; then
            local existing_probe_num
            existing_probe_num=$(wc -l < "${BEST_DIRS[$class]%%#*}/cdhit_final.fasta" | awk '{print $1 / 2}')
            if [[ "$probe_num" -lt "$existing_probe_num" ]]; then
                BEST_DIRS["$class"]="$dir#$full_class"
            fi
        fi
    fi
}

# Evaluate all configuration directories
for dir in "${OUTPUT_DIR}"/*/; do
    dir=${dir%/}
    dir_basename=$(basename "$dir")

    fasta_file="${dir}/cdhit_final.fasta"
    coverage_file="${dir}/rRNA_coverage_summary.txt"

    if [[ ! -f "$coverage_file" ]]; then
        echo "Skipping directory (no coverage file): $dir_basename"
        continue
    fi

    full_class_name=""; high_class_name=""
    extract_classes "$dir_basename"
    calculate_stats "$fasta_file" "$coverage_file"

    echo -e "$full_class_name\t$high_class_name\t$probe_num\t$q1\t$q2\t$q3\t$score" >> "$COV_STAT"
    update_max_scores "$high_class_name" "$score" "$dir" "$full_class_name"
done

# Remove non-optimal configurations
for dir in "${OUTPUT_DIR}"/*/; do
    dir=${dir%/}
    dir_basename=$(basename "$dir")
    full_class_name="$dir_basename"
    extract_classes "$full_class_name"

    if [[ "${BEST_DIRS[$high_class_name]%%#*}" != "$dir" ]]; then
        echo "Removing suboptimal directory: $dir"
        rm -rf "$dir"
    fi
done

# Write summary statistics
{
    echo -e "\n*** Summary of Max Scores ***"
    echo -e "High_Class\tMax_Score\tBest_Full_Class\tBest_Dir"
    for high_class in "${!MAX_SCORES[@]}"; do
        echo -e "$high_class\t${MAX_SCORES[$high_class]}\t${BEST_DIRS[$high_class]#*#}\t${BEST_DIRS[$high_class]%%#*}"
    done
} >> "$COV_STAT"

echo "[$(date '+%F %T')]: Summary written to $COV_STAT"

# =============================================================================
# STEP 4: Group probes into core and plus sets
# =============================================================================
for dir in "${OUTPUT_DIR}"/*/; do
    dir=${dir%/}
    rrna_type=$(echo "$(basename "$dir")" | cut -d'_' -f1)
    valid_types=("5S" "5" "18S" "28S")

    if [[ ! " ${valid_types[@]} " =~ " ${rrna_type} " ]]; then
        continue
    fi

    declare -A RNA_IDENTIFY=(
        ["5S"]="${SEP_5S}"
        ["5"]="${SEP_5_8S}"
        ["18S"]="${SEP_18S}"
        ["28S"]="${SEP_28S}"
    )

    rna_file="${RNA_IDENTIFY[$rrna_type]}"
    if [[ -z "$rna_file" ]]; then
        echo "Error: No RNA file for type '$rrna_type'. Skipping..." >&2
        continue
    fi

    final_file="${dir}/cdhit_final.fasta"
    align_res="${dir}/cdhit_final.alig.out"
    group_file="${dir}/probe_groups.tsv"
    group_filter="${dir}/probe_filtered.tsv"
    probe_file="${dir}/rrna_probes.fasta"
    tree_report="${dir}/fungal_tree_report.txt"
    ratio_threshold=0.5

    echo "[$(date '+%F %T')]: Grouping probes for $rrna_type..."

    if [[ "$rrna_type" == "18S" || "$rrna_type" == "28S" ]]; then
        # Use standard grouping for 18S/28S
        $SIF_CMD python "${GROUP_PROBES}" "${align_res}" "${rna_file}" \
            "${FUNGI_TAX}" "${group_file}" \
            "${group_filter}" \
            --ratio_threshold "$ratio_threshold" \
            --num_workers "$THREADS"

        # Generate reverse-complement (RNA-to-DNA) probes
        $SIF_CMD seqkit seq -r -p --rna2dna "${final_file}" > "${probe_file}"

        # Evaluate coverage at each taxonomic level
        echo "[$(date '+%F %T')]: Calculating level coverage..."
        declare -a LEVELS=("subkingdom" "phylum" "subphylum" "class" "subclass" "order" "suborder" "family" "subfamily")

        for level in "${LEVELS[@]}"; do
            echo "[$(date '+%F %T')]: Coverage for $level..."
            grp_res="${dir}/${rrna_type}_${level}"

            $SIF_CMD python "${GROUPED_CAL}" \
                "${group_filter}" "$ratio_threshold" "$level" \
                "${rna_file}" "${align_res}" \
                "${grp_res}_grouped.tsv" \
                "${probe_file}" \
                "${grp_res}"
        done

    elif [[ "$rrna_type" == "5S" || "$rrna_type" == "5" ]]; then
        # Use advanced grouping for 5S/5.8S (with NCBI taxonomy nodes)
        $SIF_CMD python "${GROUP_PROBES_PLUS}" "${align_res}" "${rna_file}" \
            "${FUNGI_TAX}" "${FUNGI_NODE}" \
            "${group_file}" \
            "${group_filter}" \
            --ratio_threshold "$ratio_threshold" \
            --num_workers "$THREADS"

        $SIF_CMD seqkit seq -r -p --rna2dna "${final_file}" > "${probe_file}"

        echo "[$(date '+%F %T')]: Calculating level coverage..."
        declare -a LEVELS=("subkingdom" "phylum" "subphylum" "class" "subclass" "order" "suborder" "family" "subfamily")

        for level in "${LEVELS[@]}"; do
            echo "[$(date '+%F %T')]: Coverage for $level..."
            grp_res="${dir}/${rrna_type}_${level}"

            $SIF_CMD python "${GROUPED_CAL_PLUS}" \
                "${group_filter}" "$ratio_threshold" "$level" \
                "${rna_file}" "${align_res}" \
                "${FUNGI_TAX}" "${FUNGI_NODE}" \
                "${probe_file}" \
                "${grp_res}_grouped.tsv" \
                "${grp_res}"
        done
    else
        echo "Error: Unsupported RNA type: $rrna_type" >&2
        exit 1
    fi

    # Generate taxonomic decision tree report
    $SIF_CMD python "${DEC_TREE}" \
        --nodes "${FUNGI_NODE}" \
        --info "${FUNGI_TAX}" \
        --grouped-dir "${dir}" \
        --probes-fasta "${probe_file}" \
        --species-map-fasta "${rna_file}" \
        --fungi-root-name "Fungi" \
        --output-file "${tree_report}"

done

echo "[$(date '+%F %T')]: Pipeline completed successfully!"
