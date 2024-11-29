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

# Variablesokmijnuhb987654321
data_path="/NFSHOME/lmasci/DNABERT_2/sample_BRCATOT_500bp"  
lr=0.00001

#            --eval_steps 200 \
#            --save_steps 200 \
#           loadbestmodel added

srun python /NFSHOME/lmasci/DNABERT_2/finetune/train_running.py \
            --model_name_or_path zhihan1996/DNABERT-2-117M \
            --data_path ${data_path} \
            --kmer -1 \
            --model_max_length 512 \
            --per_device_train_batch_size 8 \
            --per_device_eval_batch_size 16 \
            --gradient_accumulation_steps 1 \
            --learning_rate ${lr} \
            --num_train_epochs 5 \
            --output_dir /NFSHOME/lmasci/DNABERT_2/finetune/output/BRCATOT_500bp_lr1 \
            --evaluation_strategy steps \
            --warmup_steps 50 \
            --eval_steps 200 \
            --save_steps 200 \
            --logging_steps 100000 \
            --overwrite_output_dir True \
            --log_level info \
            --find_unused_parameters False \
            --save_model True 


