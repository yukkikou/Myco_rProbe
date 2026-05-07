#!/bin/bash
#SBATCH -A ylab
#SBATCH -J cBlast
#SBATCH -p mix,fatnode
#SBATCH -D /media/share/data6/1610305236
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o /home/1610305236/log/%x-%A_%a.out
#SBATCH -e /home/1610305236/log/%x-%A_%a.err
#SBATCH --mem 12G

# ==================== Configuration ====================
outdir=/media/share/data5/1610305236/Panfungi/script/Myco_rProbe/batch_mapping
saminfo=$outdir/../example/batch_lib_config.tsv

sifdir="/media/share/data4/container"
blastSif="$sifdir/blast/blast-2.13.0.sif"

srcdir="/media/share/data5/1610305236/Panfungi/script/Myco_rProbe/src"
blast_script="$srcdir/circ_blast.sh"

thread=$SLURM_CPUS_PER_TASK
task_id=$SLURM_ARRAY_TASK_ID

# ==================== Main Script ====================
echo "***************************************************"
echo "*** Job Name:    $SLURM_JOB_NAME"
echo "*** Job ID:      $SLURM_JOB_ID"
echo "*** Array Task:  $SLURM_ARRAY_TASK_ID"
echo "*** Node:        $SLURMD_NODENAME"
echo "*** Start Time:  $(date)"
echo "***************************************************"

# Check input file existence
if [ ! -f "$saminfo" ]; then
    echo "ERROR: Input file $saminfo not found."
    exit 1
fi

# Check blast script existence
if [ ! -f "$blast_script" ]; then
    echo "ERROR: BLAST script $blast_script not found."
    exit 1
fi

# Get specific line for this array task
line=$(sed -n "${task_id}p" "$saminfo" | tr -d '\r')

if [ -z "$line" ]; then
    echo "ERROR: Task ID $task_id is out of bounds or file $saminfo is empty."
    exit 1
fi

echo "=== Task Information ==="
echo "Task ID  : $task_id"
echo "Threads  : $thread"
echo "Output   : $outdir"
echo "Input Line: $line"
echo "========================"

# ==================== Process circRNA sequences ====================
echo "[$(date '+%F %T')]: Starting circRNA sequence extraction and formatting"

# Count valid entries and collect fasta files
valid_fasta_files=()
strain_count=0

while IFS=$'\t' read -r strain prefix raw1 raw2 genome gtf; do
    # Skip header line if exists
    if [[ "$strain" == "strain" ]] || [[ "$strain" =~ ^#.* ]]; then
        continue
    fi
    
    # Skip empty lines
    if [ -z "$strain" ] || [ -z "$prefix" ]; then
        continue
    fi
    
    ((strain_count++))
    
    echo "[$(date '+%F %T')]: Processing strain $strain (prefix: $prefix)"
    
    # Define file paths
    circ_fa="$outdir/$prefix/3_circ/circ/${prefix}_circRNA_index.fa"
    formated_fa="$outdir/$prefix/3_circ/circ/${prefix}_circRNA_index.fasta"
    
    # Check if source file exists
    if [ ! -f "$circ_fa" ]; then
        echo "[$(date '+%F %T')]: WARNING: circRNA file not found for $prefix at: $circ_fa"
        continue
    fi
    
    # Format fasta headers - add prefix to avoid duplicate sequence names
    echo "[$(date '+%F %T')]: Formatting $circ_fa"
    sed "s/^>/>${prefix}_/" "$circ_fa" > "$formated_fa"
    
    if [ -f "$formated_fa" ]; then
        valid_fasta_files+=("$formated_fa")
        echo "[$(date '+%F %T')]: Successfully formatted: $formated_fa"
    else
        echo "[$(date '+%F %T')]: ERROR: Failed to create formatted file: $formated_fa"
    fi
    
done < <(sed 's/\r$//' "$saminfo")

echo "[$(date '+%F %T')]: Processed $strain_count strains, found ${#valid_fasta_files[@]} valid fasta files"

# Check if we have any valid files to process
if [ ${#valid_fasta_files[@]} -eq 0 ]; then
    echo "ERROR: No valid circRNA fasta files found. Exiting."
    exit 1
fi

# ==================== Merge all fasta files ====================
echo "[$(date '+%F %T')]: Merging all fasta files"

circ_homo_dir="$outdir/circ_homo"
mkdir -p "$circ_homo_dir"
cd "$circ_homo_dir" || exit 1

# Clear previous merged file if exists
[ -f "circ_all.fasta" ] && rm -f "circ_all.fasta"

# Merge all formatted fasta files
cat "${valid_fasta_files[@]}" > circ_all.fasta

# Check if merge was successful
if [ ! -s "circ_all.fasta" ]; then
    echo "ERROR: Failed to create merged fasta file or file is empty"
    exit 1
fi

total_seqs=$(grep -c "^>" circ_all.fasta)
echo "[$(date '+%F %T')]: Merged fasta contains $total_seqs sequences"

# ==================== Create BLAST database ====================
echo "[$(date '+%F %T')]: Creating BLAST database"

# Remove old database files if exist
rm -f circ_homo.n*

# Create BLAST database using singularity
singularity exec -B /media:/media $blastSif makeblastdb \
    -in circ_all.fasta \
    -dbtype nucl \
    -parse_seqids \
    -hash_index \
    -out circ_homo

# Check if database was created successfully
if [ ! -f "circ_homo.nhr" ] || [ ! -f "circ_homo.nin" ] || [ ! -f "circ_homo.nsq" ]; then
    echo "ERROR: Failed to create BLAST database"
    exit 1
fi

echo "[$(date '+%F %T')]: BLAST database created successfully"

# ==================== Prepare BLAST jobs ====================
echo "[$(date '+%F %T')]: Preparing BLAST jobs"

# Create directory for BLAST results
blast_results_dir="$circ_homo_dir/blast_results"
mkdir -p "$blast_results_dir"

# Calculate number of chunks based on number of files
num_chunks=${#valid_fasta_files[@]}
echo "[$(date '+%F %T')]: Will submit $num_chunks BLAST jobs (one per strain)"

# Create a file list for the array job
fasta_list_file="$circ_homo_dir/fasta_list.txt"
printf "%s\n" "${valid_fasta_files[@]}" > "$fasta_list_file"

# ==================== Submit BLAST array jobs ====================
echo "[$(date '+%F %T')]: Submitting BLAST array jobs"

SUB_BLAST_JID=$(sbatch \
    --job-name=cblast \
    -p mix,fatnode \
    --output=/home/1610305236/log/circ_blast_%A_%a.out \
    --error=/home/1610305236/log/circ_blast_%A_%a.err \
    --nodes=1 \
    --cpus-per-task=4 \
    --mem=12G \
    --array=1-"$num_chunks"%10 \
    --export=ALL,BLAST_DB="$circ_homo_dir/circ_homo",FASTA_LIST="$fasta_list_file",RESULTS_DIR="$blast_results_dir",BLAST_SIF="$blastSif" \
    "$blast_script" \
    | awk '{print $4}')

if [ -z "$SUB_BLAST_JID" ]; then
    echo "ERROR: Failed to submit BLAST subtasks. Exiting."
    exit 1
fi

echo "[$(date '+%F %T')]: Submitted BLAST subtasks with Job ID: $SUB_BLAST_JID"

# ==================== Wait for BLAST jobs to complete ====================
echo "[$(date '+%F %T')]: Waiting for BLAST subtasks to complete..."

# Submit a wait job that depends on BLAST completion
WAIT_JOB=$(sbatch \
    --dependency=afterok:${SUB_BLAST_JID} \
    --job-name=wait_blast \
    --wrap="echo 'BLAST jobs completed at $(date)'" \
    -p mix,fatnode \
    --output=/home/1610305236/log/wait_blast_%j.out \
    --error=/home/1610305236/log/wait_blast_%j.err \
    | awk '{print $4}')

if [ -z "$WAIT_JOB" ]; then
    echo "ERROR: Failed to submit wait job"
    exit 1
fi

echo "[$(date '+%F %T')]: Wait job ID: $WAIT_JOB"

# Monitor wait job status
while true; do
    if ! squeue -j "$WAIT_JOB" 2>/dev/null | grep -q "$WAIT_JOB"; then
        break
    fi
    echo -n "."
    sleep 60
done
echo

# Check final status of BLAST jobs
echo "[$(date '+%F %T')]: Checking BLAST job status"

WAIT_JOB_STATE=$(sacct -j "$WAIT_JOB" --format=State --noheader | tail -1 | xargs)

if [ "$WAIT_JOB_STATE" != "COMPLETED" ]; then
    echo "ERROR: BLAST subtasks may have failed!"
    echo "Wait job state: $WAIT_JOB_STATE"
    
    echo "BLAST job details:"
    sacct -j "$SUB_BLAST_JID" --format=JobID,State,ExitCode,Elapsed,NodeList
    
    # Check for any completed results
    completed_count=$(ls -1 "$blast_results_dir"/*.blast 2>/dev/null | wc -l)
    echo "Completed BLAST results: $completed_count / $num_chunks"
    
    exit 1
fi

echo "[$(date '+%F %T')]: All BLAST subtasks completed successfully"

# ==================== Merge BLAST results ====================
echo "[$(date '+%F %T')]: Merging BLAST results"

merged_blast="$circ_homo_dir/all_vs_all_blast_results.txt"
filter_blast="$circ_homo_dir/all_vs_all_blast_results.txt"

cat "$blast_results_dir"/*.blast > "$merged_blast"

if [ ! -s "$merged_blast" ]; then
    echo "WARNING: Merged BLAST results file is empty"
else
    result_lines=$(wc -l < "$merged_blast")
    echo "[$(date '+%F %T')]: Merged BLAST results contain $result_lines lines"
fi

awk '{
    split($1, a, ":");
    split($2, b, ":");
    split(a[1], species1, "_");
    split(b[1], species2, "_");
    
    pident = $3;
    align_len = $4;
    qlen = $13;
    slen = $14;
    
    query_cov = (align_len / qlen) * 100;
    target_cov = (align_len / slen) * 100;
    
    if (species1[1] != species2[1] && pident > 40 && query_cov > 40 && target_cov > 40) {
        print $0
    }
}' $merged_blast > $filter_blast

# ==================== Cleanup and Summary ====================
echo "***************************************************"
echo "*** Job Completed Successfully"
echo "*** End Time: $(date)"
echo "*** Results Directory: $circ_homo_dir"
echo "*** Merged BLAST Results: $merged_blast"
echo "***************************************************"

