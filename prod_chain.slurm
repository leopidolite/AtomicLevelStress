#!/bin/bash
#SBATCH -J prod_chain
#SBATCH -p serial
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -N 1 -n 24
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -q normal
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=leopold_li@brown.edu

module load /home/leopold/mylammps/build/lmp
module load openmpi
module load gcc

LAMMPS_EXE="/home/leopold/mylammps/build/lmp"


# Production chain: run production in NVT until total completed --> Submit born matrix calculation

TEMP=${TEMP:?}
PROD_TOTAL=${PROD_TOTAL:?}
EFILE=${EFILE:?}
cd $SLURM_SUBMIT_DIR

latest() { ls $1.* 2>/dev/null | grep -oE '[0-9]+$' | sort -n | tail -1; }
next_stage() { sbatch -J born_T${TEMP} --array=0 --export=ALL,TEMPS=$TEMP born_species.slurm; }

LAST=$(latest restart_prod_T${TEMP})
DONE=${LAST:-0}
[ -n "$LAST" ] && [ "$DONE" -ge "$PROD_TOTAL" ] && { next_stage; exit 0; }

[ -z "$LAST" ] && RF=none || RF=restart_prod_T${TEMP}.$DONE
srun "${LAMMPS_EXE}" -in in.prod -log log_prod_T${TEMP}.${DONE}.lammps \
     -v temp $TEMP -v total $PROD_TOTAL -v rfile $RF -v efile $EFILE || exit 1

NEW=$(latest restart_prod_T${TEMP})
[ -n "$NEW" ] && { [ -z "$LAST" ] || [ "$NEW" -gt "$DONE" ]; } || exit 1

if [ "$NEW" -ge "$PROD_TOTAL" ]; then
    next_stage
else
    sbatch -J prod_T${TEMP} --dependency=afterok:$SLURM_JOB_ID --export=ALL,TEMP=$TEMP,EFILE=$EFILE prod_chain.slurm
fi
