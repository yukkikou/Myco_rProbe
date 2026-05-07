#!/bin/bash
#SBATCH -A ylab
#SBATCH -J cMerge
#SBATCH -p mix,fatnode
#SBATCH -D /media/share/workdir/1610305236
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -o /home/1610305236/log/%x-%A_%a.out
#SBATCH -e /home/1610305236/log/%x-%A_%a.err
#SBATCH --mem 20G

#
# Description:
#   SLURM batch job that merges circRNA data by sample group.
#   For each group (identified by the prefix before the dash in
#   sample directory names), it builds a sample list and runs
#   prep_CIRIquant inside a Singularity container.
#
# Usage:
#   Submitted via sbatch.  Expects CIRIquant Singularity image
#   and sample directories under resdir.

# --- Path definitions ---
resdir="/media/share/data5/1610305236/Panfungi/1_data/Candida"
outdir="$resdir/merge"

# Singularity image and command prefix.
sifdir="/media/share/resource/container"
cquanSif="singularity exec --bind /media:/media $sifdir/CIRIquant/ciriquant-v1.1.3.sif"

# --- Main execution ---
echo "***************************************************"
echo "*** Job Name:    $SLURM_JOB_NAME"
echo "*** Job ID:      $SLURM_JOB_ID"
echo "*** Node:        $SLURMD_NODENAME"
echo "*** Start Time:  $(date)"
echo "***************************************************"

mkdir -p "$outdir"

# Discover sample groups from directory names matching the pattern.
# Groups are identified by the prefix before the dash (e.g., C001, S002).
groups=$(ls "$resdir" | grep -E '^[CSPA][0-9]+-[0-9]+' | cut -d'-' -f1 | sort -u)

for group in $groups; do
    echo "Processing Group: $group"

    group_lst="$outdir/${group}_samples.lst"
    group_workdir="$outdir/${group}"
    group_prefix="$group_workdir/${group}"

    mkdir -p "$group_workdir"

    # Build the sample list for this group by scanning matching directories.
    > "$group_lst"
    for sample_dir in "$resdir/${group}-"*; do
        if [ -d "$sample_dir" ]; then
            sid=$(basename "$sample_dir")
            gtf=$(ls "$sample_dir/3_circ/"*_circRNA.gtf 2>/dev/null | head -n 1)

            if [ -n "$gtf" ] && [ -f "$gtf" ]; then
                echo -e "$sid\t$gtf\t$group" >> "$group_lst"
            else
                echo "Warning: No GTF found in $sample_dir/3_circ/"
            fi
        fi
    done

    # Run prep_CIRIquant if the sample list is non-empty.
    if [ -s "$group_lst" ]; then
        echo "Running prep_CIRIquant for $group..."
        $cquanSif prep_CIRIquant \
            -i "$group_lst" \
            --lib "${group_prefix}_library_info.csv" \
            --circ "${group_prefix}_circRNA_info.csv" \
            --bsj "${group_prefix}_circRNA_bsj.csv" \
            --ratio "${group_prefix}_circRNA_ratio.csv" \
            &> "${group_prefix}_CIRIquant.log"
    else
        echo "Error: Empty list for $group, skipping..."
    fi

    echo "Group $group finished."
    echo "---------------------------"
done
