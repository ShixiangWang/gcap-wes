#!/bin/bash
#SBATCH -N 1
#SBATCH -o output-%J.o
#SBATCH -n 80
#SBATCH --mem=300gb
#SBATCH -p fat-8260-Partition
#SBATCH -J wsx
Rscript 02-hyper-train.R $1 $2
