#!/bin/bash
#SBATCH -A ylab 
#SBATCH -J cQuan
#SBATCH -p mix,fatnode
#SBATCH -D /media/share/workdir/1610305236
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH -o /home/1610305236/log/%x-%A_%a.out
#SBATCH -e /home/1610305236/log/%x-%A_%a.err
#SBATCH --mem 240G
#SBATCH --array 20

outdir=/media/share/data5/1610305236/Panfungi/1_data/Cryptococcus
saminfo=/media/share/data5/1610305236/Panfungi/script/Myco_rProbe/example/cryp_lib_config.tsv
#saminfo=/media/share/data5/1610305236/Panfungi/script/Myco_rProbe/example/cand_lib_config.tsv

sifdir="/media/share/resource/container"
cquanSif="$sifdir/CIRIquant/ciriquant-v1.1.3.sif"

srcdir="/media/share/data5/1610305236/Panfungi/script/Myco_rProbe"
circ_src="$srcdir/circ_tmp2.sh"
#circ_src="$srcdir/myco_rprobe_circiden.sh"

thread=$SLURM_CPUS_PER_TASK
task_id=$SLURM_ARRAY_TASK_ID

# main
echo "***************************************************"
echo "*** Job Name:    $SLURM_JOB_NAME"
echo "*** Job ID:      $SLURM_JOB_ID"
echo "*** Array Task:  $SLURM_ARRAY_TASK_ID"
echo "*** Node:        $SLURMD_NODENAME"
echo "*** Start Time:  $(date)"
echo "***************************************************"

line=$(sed -n "${task_id}p" "$saminfo" | tr -d '\r')

if [ -z "$line" ]; then
    echo "ERROR: Task ID $task_id is out of bounds or file $saminfo is empty/not found."
    exit 1
fi

echo "=== Task Information ==="
echo "Task ID  : $task_id"
echo "Threads  : $thread"
echo "Output   : $outdir"
echo "Input Line: $line"
echo "========================"

mkdir -p $outdir && cd $outdir

TASK_INPUT_FILE="${outdir}/task_${task_id}_input.info"
echo "$line" > "$TASK_INPUT_FILE"

if [ ! -s "$TASK_INPUT_FILE" ]; then
    echo "ERROR: Failed to create or write to temporary input file: $TASK_INPUT_FILE"
    exit 1
fi

echo "Temporary input file for task $task_id created at: $TASK_INPUT_FILE"

bash "$circ_src" -g "$TASK_INPUT_FILE"  -t "$thread" -o "$outdir" -m "$cquanSif"

EXIT_CODE=$?

rm "$TASK_INPUT_FILE"

echo "[$(date '+%F %T')] Task $task_id finished with exit code $EXIT_CODE."
