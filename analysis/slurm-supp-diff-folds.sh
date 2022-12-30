#!/bin/bash
#SBATCH -N 1
#SBATCH -o output-%J.o
#SBATCH -n 40
#SBATCH --mem=200gb
#SBATCH -p cpuPartition  #fat-4820-Partition
#SBATCH -J wsx

Rscript supp-cv-with-different-folds.R $1

