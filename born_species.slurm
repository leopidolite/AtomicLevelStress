#!/bin/bash
#SBATCH -J bornpa
#SBATCH -p serial
#SBATCH -q normal
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o born.out
#SBATCH -e born.err
#SBATCH --array=0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leopold_li@brown.edu

module purge
module load gcc
module load openmpi

LAMMPS_EXE=/home/leopold/mylammps/build/lmp
SCRATCH=/scratch/hpcl-mat269/lls/PRODUCTION_DATA

temps=(1500)
# temps) 
T=${temps[$SLURM_ARRAY_TASK_ID]}

srun "${LAMMPS_EXE}" -in in.born_species \
     -var temp ${T} \
     -var scratch ${SCRATCH} \
     -var delta 1.0e-3 \
     -var every 1000 \
     -var nframe 1000
