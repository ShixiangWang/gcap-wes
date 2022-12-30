setwd(file.path(PROJ_DIR, "modeling"))

library(UCSCXenaShiny)
library(dplyr)
library(IDConverter) # https://github.com/ShixiangWang/IDConverter

load("../preprocessing-and-EDA/data/pancan_amplicon_list_and_summary.RData")

# Obtain basic clinical information ---------------------------------------

samp_info <- data_summary_tcga %>%
  dplyr::select(sample_barcode) %>%
  # Use case id instead of sample id to fully match all samples
  dplyr::mutate(case_id = substr(sample_barcode, 1, 12)) %>%
  dplyr::left_join(tcga_clinical %>%
                     dplyr::select(patient, type, age_at_initial_pathologic_diagnosis, gender) %>%
                     unique(), by = c("case_id" = "patient")) %>% 
  dplyr::ungroup()

colnames(samp_info) <- c("sample", "case", "type", "age", "gender")

samp_info

table(samp_info$type)

saveRDS(samp_info, file = "data/model_tcga_amplicon_related_clinical_data.rds")

# Obtain ASCN related information -----------------------------------------

tcga_wes_cn <- readRDS("data/tcga_ascat_wes_tumor_cn.rds")

## There are some tumor-normal pairs for same tumor sample
## We use the CNV from tumor-normal pair with highest purity
keep_idx <- dplyr::tibble(
  sample = names(tcga_wes_cn$purity),
  purity = tcga_wes_cn$purity %>% as.numeric()
) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::filter(sample %in% data_summary_tcga$sample_barcode) %>% # Only keep samples with amplicon records
  dplyr::group_by(sample) %>%
  dplyr::slice_max(purity, with_ties = FALSE) %>% # Make sure only peak one pair for one sample
  dplyr::ungroup() %>%
  dplyr::pull(idx) %>% sort()

purity_ploidy <- dplyr::tibble(
  sample = names(tcga_wes_cn$purity)[keep_idx],
  purity = tcga_wes_cn$purity[keep_idx] %>% as.numeric(),
  ploidy = tcga_wes_cn$ploidy[keep_idx] %>% as.numeric(),
)

purity_ploidy
tcga_wes_cn2 = tcga_wes_cn$data
df_index = dplyr::tibble(
  source = unique(tcga_wes_cn2$source),
  idx = seq_along(source)
)

tcga_wes_cn2 = dplyr::left_join(tcga_wes_cn2, df_index, by = "source")

tcga_wes_cn2 <- tcga_wes_cn2 %>%
  dplyr::filter(idx %in% keep_idx)

all.equal(unique(tcga_wes_cn2$sample), purity_ploidy$sample)
tcga_wes_cn = tcga_wes_cn2

tcga_wes_cn$source <- tcga_wes_cn$idx <- NULL
tcga_wes_cn$chromosome <- paste0("chr", tcga_wes_cn$chromosome)
colnames(tcga_wes_cn)[4] <- "segVal"

library(sigminer)
score_aneuploidy <- sigminer::get_Aneuploidy_score(
  tcga_wes_cn,
  ploidy_df = purity_ploidy %>% dplyr::select(sample, ploidy),
  genome_build = "hg38",
  rm_black_arms = TRUE
)
score_aneuploidy <- score_aneuploidy[, .(sample, AScore)]

score_pLOH <- sigminer::get_pLOH_score(tcga_wes_cn, genome_build = "hg38")

cn <- read_copynumber(
  tcga_wes_cn,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  max_copynumber = 10000L,
  join_adj_seg = FALSE,
  genome_build = "hg38",
  add_loh = TRUE,
  loh_min_len = 1e3
)

cn_tally <- sig_tally(cn, method = "S")

saveRDS(tcga_wes_cn, file = "data/tcga_ascat_tidy_wes_tumor_cn.rds")

cn_act <- sig_fit(t(cn_tally$all_matrices$CN_48), sig_db = "CNS_TCGA", sig_index = "ALL", return_class = "data.table")

sapply(cn_act[, -1], summary)

cn@summary.per.sample[, .(sample, cna_burden)]

## merge measures
score_df <- purrr::reduce(list(
  purity_ploidy, score_aneuploidy,
  score_pLOH, cn@summary.per.sample[, .(sample, cna_burden)],
  cn_act
), dplyr::full_join, by = "sample")

sapply(score_df[, -1], function(x) sum(is.na(x)))

saveRDS(score_df, file = "data/model_tcga_train_samples_features.rds")

rm(list = setdiff(ls(), "PROJ_DIR"))
gc()

# Combine features to generate gene-based model input data

# "data/model_tcga_amplicon_related_clinical_data.rds"
# "data/model_tcga_train_samples_features.rds"
# "data/model_target_gene_refs.rds"
# "data/model_gene_amplicon_freq.rds"
# "data/model_Y_all.rds"
# "data/model_train_gene_ASCN.rds"
library(dplyr)
library(mltools)

# Obtain prior knowledge from WGS amplicon data ---------------------------
# Extract gene-level features
# We only focus on autosome

amplicon <- readRDS("../preprocessing-and-EDA/data/amplicon_hg19.rds")
## Obtain amplicon frequency from known amplicon list

# We already know total case number is *3212* in NG paper
# Extrachromosomal DNA is associated with oncogene amplification and poor outcome across multiple cancers
length(unique(amplicon$sample_barcode))

amplicon$amplicon_index <- NULL
amplicon <- amplicon[chromosome %in% paste0("chr", 1:22)]

# We only focus on protein coding genes/as they are the major target of WES
# The data source is hg19, so we use hg19 to calculate result on gene levels
ref_genes <- readRDS("../preprocessing-and-EDA/data/hg19_gene_info.rds")
ref_genes <- ref_genes[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
ref_genes$strand <- NULL
ref_genes[, gene_id := gsub("(\\..+)", "", gene_id)]
ref_genes$gene_type <- NULL

# Use gene_id is more robust
# > length(unique(ref_genes$gene_id))
# > length(unique(ref_genes$gene_name))
ref_genes$gene_name <- NULL

regions <- modules::use("../lib/regions.R")
amplicon_gene_df <- regions$overlaps(amplicon, ref_genes)

gene_freq <- amplicon_gene_df[intersect_ratio >= 1, .(freq = length(sample_barcode) / 3212L),
                              by = .(amplicon_classification, gene_id)]

# We time the freq by 1000 to represent frequency in 1000 cases.
gene_freq$freq <- gene_freq$freq * 1000L
colnames(gene_freq)[1] <- "type"
gene_freq
saveRDS(gene_freq, "data/model_gene_amplicon_freq.rds")

Y_all <- amplicon_gene_df[, .(y = max(intersect_ratio)), by = .(gene_id, sample_barcode, amplicon_classification)]
saveRDS(Y_all, file = "data/model_Y_all.rds")

# Obtain gene-level ASCN --------------------------------------------------
ref_genes_hg38 <- readRDS("../preprocessing-and-EDA/data/hg38_gene_info.rds")
ref_genes_hg38 <- ref_genes_hg38[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
ref_genes_hg38$strand <- NULL
ref_genes_hg38[, gene_id := gsub("(\\..+)", "", gene_id)]
ref_genes_hg38$gene_type <- NULL
ref_genes_hg38$gene_name <- NULL

# only keep gene id overlaps in hg19 and hg38
overlap_genes = intersect(ref_genes$gene_id, ref_genes_hg38$gene_id)
ref_genes_hg38 = ref_genes_hg38[gene_id %in% overlap_genes]

ascn_gene_df <- regions$overlaps(readRDS("data/tcga_ascat_tidy_wes_tumor_cn.rds"), ref_genes_hg38)

# nrow(unique(ascn_gene_df[, .(sample, gene_id, segVal, minor_cn)]))
# nrow(unique(ascn_gene_df[, .(sample, gene_id, segVal, minor_cn, intersect_ratio)]))
# which(duplicated(ascn_gene_df[, .(sample, gene_id, segVal, minor_cn)]))
# x <- ascn_gene_df[gene_id == "ENSG00000181143" & sample == "TCGA-02-2485-01"]
# x
# x[c(2, 1, 3)][, .(sample, gene_id, segVal, minor_cn, intersect_ratio)][
#   , .SD[order(intersect_ratio, decreasing = TRUE),][1], by = .(sample, gene_id)]

# Only keep the most overlap size segment for each gene
library(data.table)
ascn_gene_df2 <- ascn_gene_df[, .(sample, gene_id, segVal, minor_cn, intersect_ratio)][order(intersect_ratio, decreasing = TRUE)]
ascn_gene_df2 <- ascn_gene_df2[!duplicated(ascn_gene_df2[, .(sample, gene_id)])]

sum(duplicated(ascn_gene_df2[, .(sample, gene_id)]))
#saveRDS(ascn_gene_df2, file = "data/model_gene_ASCN.rds")

df_all_cli <- readRDS("data/model_tcga_amplicon_related_clinical_data.rds")
df_samp_fts <- readRDS("data/model_tcga_train_samples_features.rds")
df_amp_freq <- gene_freq
df_Y <- Y_all[sample_barcode %in% df_samp_fts$sample]
colnames(df_Y)[2:3] <- c("sample", "type")
df_gene_ascn <- ascn_gene_df2

fts_samp <- df_samp_fts %>%
  dplyr::left_join(
    df_all_cli %>%
      dplyr::select(sample, type, age, gender) %>%
      unique(),
    by = "sample"
  ) %>%
  dplyr::mutate(
    gender = ifelse(gender == "MALE", 1L, 0L)
  ) %>%
  data.table::as.data.table()

sum(is.na(fts_samp$gender))

fts_samp <- fts_samp %>%
  dplyr::select(-type) %>%
  dplyr::full_join(one_hot(fts_samp[, .(sample, type = factor(type))], sparsifyNAs = TRUE),
                   by = "sample"
  )

gene_hg19 = readRDS("../preprocessing-and-EDA/data/hg19_gene_info.rds")[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
gene_hg38 = readRDS("../preprocessing-and-EDA/data/hg38_gene_info.rds")[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
gene_hg19[, gene_id := gsub("(\\..+)", "", gene_id)]
gene_hg38[, gene_id := gsub("(\\..+)", "", gene_id)]
df_genes = data.table::data.table(gene_id = intersect(gene_hg19$gene_id, gene_hg38$gene_id))

fts_freq <- df_amp_freq %>%
  tidyr::pivot_wider(names_from = "type", values_from = "freq", values_fill = 0) %>%
  dplyr::right_join(df_genes %>% dplyr::select(gene_id), by = "gene_id") %>%
  dplyr::mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))
colnames(fts_freq)[-1] <- paste0("freq_", colnames(fts_freq)[-1])
colnames(fts_freq)[5] <- "freq_HR"

df_Y_wide <- df_Y %>%
  tidyr::pivot_wider(names_from = "type", values_from = "y", values_fill = 0)
colnames(df_Y_wide)[3:6] <- c("y_circle", "y_linear", "y_HR", "y_BFB")

data <- df_gene_ascn %>%
  dplyr::left_join(fts_samp, by = "sample") %>%
  dplyr::left_join(fts_freq, by = "gene_id") %>%
  dplyr::left_join(df_Y_wide, by = c("gene_id", "sample"))

data[, y_circle := ifelse(is.na(y_circle), 0, y_circle)]
data[, y_linear := ifelse(is.na(y_linear), 0, y_linear)]
data[, y_HR := ifelse(is.na(y_HR), 0, y_HR)]
data[, y_BFB := ifelse(is.na(y_BFB), 0, y_BFB)]

# Set up data for subsample model
data_bk = data.table::copy(data)
sum(data$y_circle)
sum(data$y_HR)
sum(data$y_linear)
# BFB may be circular amplicon, remove it from modeling
data2 = data[y_circle == 1 | (y_circle == 0 & (y_HR == 1 | y_linear == 1))]
# A focal amplicon is defined as amplification with CN above ploidy + 4
data2 = data2[segVal >= ploidy + 4]
# Randomly sample equal-sizes non-amplicon records
set.seed(2021)
data3 = dplyr::slice_sample(data[y_circle == 0 & y_HR == 0 & 
                                   y_BFB == 0 & y_linear == 0 & 
                                   (segVal < ploidy + 4)][order(sample, -segVal)], n = nrow(data2))
data2 = rbind(data2, data3)[order(sample)]

data2[, y := ifelse(y_circle == 1, 1, 0)]
data2[, c("y_circle", "y_HR", "y_BFB", "y_linear") := NULL]
colnames(data2)
table(data2$y)

data = data2
data_gene_load = data[, .(load = sum(y)), by = .(sample)][order(load, decreasing = TRUE), ]

# Order the data
data = data[order(match(sample, data_gene_load$sample))]
identical(unique(data$sample), data_gene_load$sample)

data_sample_record = rle(data$sample)
sum(data$y)
sum(data_gene_load$load)

saveRDS(data, file = "data/train_data_subsample_raw.rds")
# Exclude gene_id and sample id as feature
data$gene_id <- NULL
data$sample <- NULL
data$intersect_ratio <- NULL

data_mat <- as.matrix(data)
data_mat[1:5, ]

colnames(data_mat)[1] <- "total_cn"
saveRDS(data_mat, file = "data/train_data_subsample.rds")
save(data_sample_record, data_gene_load, file = "data/train_data_extra_info_subsample.RData")

# Version v4
data = data_bk
# BFB may be circular amplicon, remove it from modeling
data2 = data[y_circle == 1 | (y_circle == 0 & (y_HR == 1 | y_linear == 1))]
# A focal amplicon is defined as amplification with CN above ploidy + 4
data2 = data2[segVal >= ploidy + 4]
data <- rbind(data2,
              data[y_circle == 0 & y_HR == 0 & y_BFB == 0 & y_linear == 0 & (segVal < ploidy + 4)])[order(sample)]

data[, y := ifelse(y_circle == 1, 1, 0)]
data[, c("y_circle", "y_HR", "y_BFB", "y_linear") := NULL]

which(sapply(data, function(x) any(is.na(x))))
sum(is.na(fts_samp$age))

data_gene_load = data[, .(load = sum(y)), by = .(sample)][order(load, decreasing = TRUE), ]

# Order the data
data = data[order(match(sample, data_gene_load$sample))]
identical(unique(data$sample), data_gene_load$sample)

data_sample_record = rle(data$sample)
sum(data$y)
sum(data_gene_load$load)

saveRDS(data, file = "data/train_data_v4_raw.rds")
# Exclude gene_id and sample id as feature
data$gene_id <- NULL
data$sample <- NULL
data$intersect_ratio <- NULL

data_mat <- as.matrix(data)
data_mat[1:5, ]

colnames(data_mat)[1] <- "total_cn"
saveRDS(data_mat, file = "data/train_data_v4.rds")
save(data_sample_record, data_gene_load, file = "data/train_data_extra_info_v4.RData")

rm(list = setdiff(ls(), "PROJ_DIR"))
gc()
