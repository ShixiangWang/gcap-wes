#!/bin/bash

#SBATCH -p cn
#SBATCH -J cv
#SBATCH --array=1-3
#SBATCH -o output-cv-%J.o
#SBATCH -n 80
#SBATCH --mem=102400

# echo SLURM_JOB_ID $SLURM_JOB_ID                # 单个子任务 job id
# echo SLURM_ARRAY_JOB_ID $SLURM_ARRAY_JOB_ID    # 主任务 job id
# echo SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID  # 设定的 array 元素值

task_id=$SLURM_ARRAY_TASK_ID

Rscript ~/proj/ecDNA/03-modeling-renew/02-hyper-train.R $task_id
