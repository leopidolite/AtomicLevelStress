#!/bin/bash
#SBATCH -J equil_chain
#SBATCH -p serial
#SBATCH -t 24:00:00
#SBATCH -N 1 -n 24
#SBATCH -c 1
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -q normal
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=leopold_li@brown.edu

# Equilibration Chain --> equilibrate until total steps completed, NVT --> run Production chain

TEMP=${TEMP:-850}
EQ_TOTAL=3000000
PROD_TOTAL=3000000
cd $SLURM_SUBMIT_DIR
SCRIPT=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K\S+')

latest() { ls $1.* 2>/dev/null | grep -oE '[0-9]+$' | sort -n | tail -1; }
next_stage() { sbatch -J prod_T${TEMP} --export=ALL,TEMP=$TEMP,PROD_TOTAL=$PROD_TOTAL,EFILE=restart_equil_T${TEMP}.$1 prod_chain.slurm; }

DONE=$(latest restart_equil_T${TEMP}); DONE=${DONE:-0}
[ "$DONE" -ge "$EQ_TOTAL" ] && { next_stage $DONE; exit 0; }

[ "$DONE" -eq 0 ] && RF=none || RF=restart_equil_T${TEMP}.$DONE
srun lmp -in in.equil -log none -v temp $TEMP -v total $EQ_TOTAL -v rfile $RF || exit 1

NEW=$(latest restart_equil_T${TEMP}); NEW=${NEW:-0}
[ "$NEW" -gt "$DONE" ] || exit 1

if [ "$NEW" -ge "$EQ_TOTAL" ]; then
    next_stage $NEW
else
    sbatch -J equil_T${TEMP} --dependency=afterok:$SLURM_JOB_ID --export=ALL,TEMP=$TEMP $SCRIPT
fi
