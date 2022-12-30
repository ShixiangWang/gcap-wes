# Ref: https://github.com/ShixiangWang/sigminer/blob/master/data-raw/transcript_region.R
setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(data.table)
library(stringr)
library(magrittr)

extract_col <- function(x, name) {
  str_extract(x, paste0(name, " ([^;]+);")) %>%
    str_remove(paste0(name, " ")) %>%
    str_remove_all("\"") %>%
    str_remove(";")
}

# The raw gtf data is downloaded from GENCODE
gtf_hg38 <- data.table::fread("~/data/refs/gencode.v38.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

exon_hg38 <- gtf_hg38[V3 == "exon"]
gtf_hg38 <- gtf_hg38[V3 == "gene"]

gtf_hg38[, `:=`(gene_name = extract_col(V9, "gene_name"), gene_id = extract_col(V9, "gene_id"), gene_type = extract_col(V9, "gene_type"))]

gene_hg38 <- gtf_hg38[, .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_hg38)[1:4] <- c("chrom", "start", "end", "strand")

length(unique(gene_hg38$gene_name))
length(unique(gene_hg38$gene_id))



# https://bedtools.readthedocs.io/en/latest/content/general-usage.html Use gene_id for the name
gene_bed <- gene_hg38[, .(chrom, start, end, gene_id)]
gene_bed[, `:=`(start, start - 1)]

data.table::fwrite(gene_bed, "data/hg38_genes.bed", col.names = FALSE, sep = "\t")
saveRDS(gene_hg38, "data/hg38_gene_info.rds")

exon_hg38[, `:=`(gene_name = extract_col(V9, "gene_name"), exon_id = extract_col(V9, "exon_id"), exon_number = extract_col(V9, "exon_number"))]

exon_hg38 <- exon_hg38[, .(V1, V4, V5, gene_name, exon_id, exon_number)]
colnames(exon_hg38)[1:3] <- c("chrom", "start", "end")
exon_hg38$exon_number <- as.integer(exon_hg38$exon_number)

length(unique(exon_hg38$exon_id))
nrow(exon_hg38)

which(duplicated(exon_hg38$exon_id[1:1000]))
exon_hg38[exon_id == exon_id[47]]

# Remove duplication
exon_hg38_2 <- unique(exon_hg38[, 1:5])

which(duplicated(exon_hg38_2$exon_id))
exon_hg38_2[exon_id == exon_id[762416]]

# Remove the duplication in Y
exon_hg38_3 <- exon_hg38_2[!duplicated(exon_hg38_2$exon_id)]
nrow(exon_hg38_3)
length(unique(exon_hg38$exon_id))

#saveRDS(exon_hg38_3, file = "data/hg38_exon_info.rds")

# hg19 --------------------------------------------------------------------
# https://www.gencodegenes.org/human/release_38lift37.html

gtf_hg19 <- data.table::fread("~/data/refs/gencode.v38lift37.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)
gtf_hg19 <- gtf_hg19[V3 == "gene"]

gtf_hg19[, `:=`(gene_name = extract_col(V9, "gene_name"), gene_id = extract_col(V9, "gene_id"), gene_type = extract_col(V9, "gene_type"))]

gene_hg19 <- gtf_hg19[, .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_hg19)[1:4] <- c("chrom", "start", "end", "strand")

length(unique(gene_hg19$gene_name))
length(unique(gene_hg19$gene_id))

# https://bedtools.readthedocs.io/en/latest/content/general-usage.html Use gene_id for the name
gene_bed <- gene_hg19[, .(chrom, start, end, gene_id)]
gene_bed[, `:=`(start, start - 1)]

data.table::fwrite(gene_bed, "data/hg19_genes.bed", col.names = FALSE, sep = "\t")
saveRDS(gene_hg19, "data/hg19_gene_info.rds")
