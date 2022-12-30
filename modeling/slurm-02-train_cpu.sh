#!/bin/bash
#SBATCH -N 1
#SBATCH -o output-%J.o
#SBATCH -n 40
#SBATCH --mem=200gb
#SBATCH -p cpuPartition  #fat-4820-Partition
#SBATCH -J wsx
# sbatch slurm-02-train_others.sh 1
# sbatch slurm-02-train_others.sh 1 rev
Rscript 02-hyper-train.R $1 $2

