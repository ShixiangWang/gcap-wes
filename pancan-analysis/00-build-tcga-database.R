library(UCSCXenaShiny)
library(data.table)
library(IDConverter)

# Use ENSEMBL id is better than symbol for join different data types.

tcga = list()
tcga$survival = tcga_surv

tcga_purity$cancer_type = NULL
tcga$purity = tcga_purity
tcga$instability = tcga_genome_instability

# gene non-silent mutation ----------------------
mutation = data.table::fread("/data3/wsx_data/Xena/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz")
mutation = mutation[!is.na(sample), lapply(.SD, function(x) max(x, na.rm = TRUE)), by = .(sample)] %>% 
  tibble::column_to_rownames("sample") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample")

data.table::setDT(mutation)
# Remove genes without mutation
mutation2 = mutation[, lapply(.SD, function(x) {
  if (is.numeric(x)) {
    if (sum(x) > 0) {
      x
    } else {
      NULL
    }
  } else {
    x
  }
})]

geneMap = data.table::fread("/data3/wsx_data/Xena/probeMap%2Fgencode.v23.annotation.gene.probemap")
geneMap[, id := substr(id, 1, 15)]
new_ids = convert_custom(colnames(mutation2)[-1], "gene", "id", dt = geneMap)
sum(is.na(new_ids))
mutation2 = mutation2[, c(TRUE, !is.na(new_ids)), with = FALSE]
colnames(mutation2)[-1] = convert_custom(colnames(mutation2)[-1], "gene", "id", dt = geneMap)

tcga$mutation = mutation2

rm(mutation, mutation2); gc()

# Gene expression ---------------------------------------------------------
expr = data.table::fread("/data3/wsx_data/Xena/tcga_RSEM_gene_tpm.gz")
expr[, sample := substr(sample, 1, 15)]

sum(duplicated(expr$sample))

expr2 = melt(expr, id.vars = "sample")[
  !is.na(sample), .(value = mean(value, na.rm = TRUE)), by = .(sample, variable)] %>% 
  dcast(variable ~ sample)
colnames(expr2)[1] = "sample"

tcga$expression = expr2

rm(expr2, expr); gc()

# Copy number variation ----------------------------------------------------

cnv = data.table::fread("/data3/wsx_data/Xena/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz")
cnv_thresholded = data.table::fread("/data3/wsx_data/Xena/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")

cnv[1:10, 1:10]
cnv_thresholded[1:10, 1:10]

sum(duplicated(cnv$Sample))
sum(duplicated(cnv_thresholded$Sample))

cnv[, Sample := convert_custom(Sample, "gene", "id", dt = geneMap)]
cnv = cnv[!is.na(Sample)]

cnv_thresholded[, Sample := convert_custom(Sample, "gene", "id", dt = geneMap)]
cnv_thresholded = cnv_thresholded[!is.na(Sample)]

cnv2 = melt(cnv, id.vars = "Sample")[!is.na(Sample)] %>% 
  dcast(variable ~ Sample)
colnames(cnv2)[1] = "sample"

cnv_thresholded2 = melt(cnv_thresholded, id.vars = "Sample")[!is.na(Sample)] %>% 
  dcast(variable ~ Sample)
colnames(cnv_thresholded2)[1] = "sample"

cnv2[1:10, 1:10]
cnv_thresholded2[1:10, 1:10]

tcga$cnv = cnv2
tcga$cnv_thresholded = cnv_thresholded2

rm(cnv, cnv2, cnv_thresholded, cnv_thresholded2); gc()

# Protein expression ------------------------------------------------------

prot = data.table::fread("/data3/wsx_data/Xena/TCGA-RPPA-pancan-clean.xena.gz")
prot[1:10, 1:10]

sum(duplicated(prot$SampleID))

prot2 = melt(prot, id.vars = "SampleID")[!is.na(SampleID)] %>% 
  dcast(variable ~ SampleID, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(prot2)[1] = "sample"

prot2[1:10, 1:10]

tcga$protein = prot2

rm(prot, prot2); gc()

# Immune ------------------------------------------------------------------

immune_subtype = data.table::fread("/data3/wsx_data/Xena/Subtype_Immune_Model_Based.txt.gz")
immune_signature = data.table::fread("/data3/wsx_data/Xena/TCGA_pancancer_10852whitelistsamples_68ImmuneSigs.xena.gz")
  
immune_signature[1:10, 1:10]
immune_signature$V1

immune_signature2 = melt(immune_signature, id.vars = "V1")[!is.na(V1)] %>% 
  dcast(variable ~ V1, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(immune_signature2)[1] = "sample"

tcga$immune =  merge(immune_subtype, immune_signature2, by = "sample", all = TRUE)

rm(immune_signature, immune_signature2, immune_subtype); gc()

# Stemness ----------------------------------------------------------------

stem_rna = data.table::fread("/data3/wsx_data/Xena/StemnessScores_RNAexp_20170127.2.tsv.gz")
stem_rna[, 1:10]

stem_mey = data.table::fread("/data3/wsx_data/Xena/StemnessScores_DNAmeth_20170210.tsv.gz")
stem_mey[, 1:10]

stemness = merge(
  melt(stem_rna, id.vars = "sample")[!is.na(sample)] %>% 
    dcast(variable ~ sample, fun.aggregate = function(x) mean(x, na.rm = TRUE)),
  melt(stem_mey, id.vars = "sample")[!is.na(sample)] %>% 
    dcast(variable ~ sample, fun.aggregate = function(x) mean(x, na.rm = TRUE)),
  by = "variable", all = TRUE
)
colnames(stemness)[1] = "sample"

tcga$stemness = stemness

rm(stem_mey, stem_rna, stemness); gc()


# Genes/Pathways ----------------------------------------------------------------

gene_program = data.table::fread("/data3/wsx_data/Xena/Pancan12_GenePrograms_drugTargetCanon_in_Pancan33.tsv.gz")
gene_program[1:10, 1:10]

gene_program2 = melt(gene_program, id.vars = "sample")[!is.na(sample)] %>% 
  dcast(variable ~ sample, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(gene_program2)[1] = "sample"

gene_set_ssGSEA = data.table::fread("/data3/wsx_data/Xena/PanCan33_ssGSEA_1387GeneSets_NonZero_sample_level.txt.gz")
gene_set_ssGSEA[1:10, 1:10]

gene_set_ssGSEA2 = melt(gene_set_ssGSEA, id.vars = "V1")[!is.na(V1)] %>% 
  dcast(variable ~ V1, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(gene_set_ssGSEA2)[1] = "sample"

RABIT = data.table::fread("/data3/wsx_data/Xena/RABIT%2Fpancan%2FRABIT_pancan.HiSeq.V2.gz")
RABIT[1:10, 1:10]

tcga$RABIT = RABIT
tcga$gene_program = gene_program2
tcga$gene_set_ssGSEA = gene_set_ssGSEA2

rm(RABIT, gene_program, gene_set_ssGSEA, gene_set_ssGSEA2, gene_program2); gc()

# Others ------------------------------------------------------------------

HRD = data.table::fread("/data3/wsx_data/Xena/TCGA.HRD_withSampleID.txt.gz")
HRD[1:4, 1:10]

HRD2 = melt(HRD, id.vars = "sampleID")[!is.na(sampleID)] %>% 
  dcast(variable ~ sampleID, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(HRD2)[1] = "sample"
HRD2[1:10]

iCluster = data.table::fread("/data3/wsx_data/Xena/TCGA_PanCan33_iCluster_k28_tumor.gz")
iCluster[1:10]

Subtype =  data.table::fread("/data3/wsx_data/Xena/TCGASubtype.20170308.tsv.gz")
Subtype = merge(Subtype, iCluster, by.x = "sampleID", by.y = "sample", all = TRUE)

tcga$HRD = HRD2
tcga$Subtype = Subtype

rm(HRD, HRD2, Subtype); gc()

tcga$clinical = as.data.table(tcga_clinical)

# Gene methylation --------------------------------------------------------

methy = data.table::fread("/data3/wsx_data/Xena/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena", showProgress = TRUE)
loc_map = data.table::fread("/data3/wsx_data/Xena/probeMap%2FilluminaMethyl450_hg19_GPL16304_TCGAlegacy")
loc_map = loc_map[!is.na(gene) & gene != "."]
colnames(loc_map)[1] = "sample"
loc_map2 = loc_map[, .(sample, gene)] %>% 
  tidyr::separate_rows(gene, sep = ",") %>% 
  dplyr::mutate(gene = convert_custom(gene, "gene", "id", dt = geneMap)) %>% 
  dplyr::filter(!is.na(gene))

methy = methy[sample %in% loc_map2$sample]
methy = methy[, !duplicated(colnames(methy)), with = FALSE]

gc()
methy = merge(methy, as.data.table(loc_map2), by = "sample", all.x = TRUE)
gc()

methy[1:10, 1:10]
methy_mean = methy[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(gene), .SDcols = colnames(methy)[!colnames(methy) %in% c("gene", "sample")]]
#methy_median = methy[, lapply(.SD, function(x) median(x, na.rm = TRUE)), by = .(gene), .SDcols = colnames(methy)[!colnames(methy) %in% c("gene", "sample")]]

sum(duplicated(methy_mean$gene))
methy_mean2 = melt(methy_mean, id.vars = "gene")[!is.na(gene)] %>% 
  dcast(variable ~ gene, fun.aggregate = function(x) mean(x, na.rm = TRUE))
colnames(methy_mean2)[1] = "sample"

tcga$methylation = methy_mean2
rm(methy, methy_mean, methy_mean2); gc()

# Fusion downloaded from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753333/
fusion = readxl::read_excel("/data3/wsx_data/tcga_gene_fusion.xlsx", 1)
fusion$Sample = substr(fusion$Sample, 1, 15) # Keep gene name here, convert if when necessary

tcga$fusion = fusion

save(tcga, file = "/data3/wsx_data/tcga.RData")

# # Query method
# library(dplyr)
# library(dbplyr)
# 
# tcga = dbConnect(drv = RSQLite::SQLite(), "/data3/wsx_data/tcga.sqlite")
# tbl(tcga, "survival") %>% 
#   select(OS) %>% 
#   collect()
# 
# dbDisconnect(tcga)
