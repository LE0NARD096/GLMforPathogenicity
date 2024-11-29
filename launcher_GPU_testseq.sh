#!/bin/bash -l
#SBATCH -J dnabert
#SBATCH -o logs/%j.out
#SBATCH -n 1
#SBATCH -p cuda
#SBATCH -c 8
#SBATCH --gres=gpu:large

# Set environment variables for threading
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

#srun python test_finetunedmodel.py
srun python test_finetunedmodel_single.py
