# This script serves a template

library(readr)
library(gcap)

setwd("~/proj/ecDNA/")
df <- read_csv("data/nonEC_wes_pair_info.csv")

tfile <- file.path("~/data/gdc", df$file_id_tumor, df$file_name_tumor)
nfile <- file.path("~/data/gdc", df$file_id_normal, df$file_name_normal)

tn <- df$tumor
nn <- df$normal
id <- df$pair_id_uniq

gcap.runASCAT(
  tumourseqfile = tfile, normalseqfile = nfile,
  tumourname = tn, normalname = nn, jobname = id,
  outdir = "data/nonEC_result",
  skip_finished_ASCAT = TRUE
)


# HPC template:
# #!/bin/bash
# #SBATCH -p cn
# #SBATCH -J nonEC
# #SBATCH -o nonEC.o
# #SBATCH -n 22
# #SBATCH --mem=102400
# 
# Rscript xx/call-gcap.R


