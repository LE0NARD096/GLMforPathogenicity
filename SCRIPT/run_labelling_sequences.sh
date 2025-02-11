#!/bin/bash
#SBATCH -J filtering_variants
#SBATCH -o /NFSHOME/lmasci/DNABERT_2/DATA/log/log.out
#SBATCH -n 4
#SBATCH -p normal
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
srun python /NFSHOME/lmasci/DNABERT_2/SCRIPT/preprocess.py
