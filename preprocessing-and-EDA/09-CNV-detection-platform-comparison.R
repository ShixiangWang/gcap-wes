setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(UCSCXenaShiny)
library(dplyr)
library(IDConverter) # https://github.com/ShixiangWang/IDConverter


# Comparison of different gene-level data from WGS/SNP/WES and the --------
# association with amplicon

load("data/pancan_amplicon_list_and_summary.RData")
target_samps <- data_summary_tcga$sample_barcode

## Obtain CNV from PCAWG TCGA samples (WGS)

pcawg_tcga_samps <- pcawg_purity %>%
  dplyr::mutate(
    submitter_id = convert_pcawg(
      icgc_specimen_id,
      from = "icgc_specimen_id",
      to = "submitted_specimen_id",
      db = "simple"
    ),
    submitter_id = substr(submitter_id, 1, 15)
  ) %>%
  dplyr::filter(startsWith(submitter_id, "TCGA")) %>%
  dplyr::select(icgc_specimen_id, submitter_id)

pcawg_tcga_map <- pcawg_tcga_samps$submitter_id
names(pcawg_tcga_map) <- pcawg_tcga_samps$icgc_specimen_id

# length(unique(pcawg_tcga_samps$submitter_id))
pcawg_cn <- data.table::fread("/data3/wsx_data/Xena/20170119_final_consensus_copynumber_sp")
pcawg_tcga_cn <- pcawg_cn[sampleID %in% pcawg_tcga_samps$icgc_specimen_id, .(sampleID, chr, start, end, total_cn, minor_cn)]
pcawg_tcga_cn$sampleID <- pcawg_tcga_map[pcawg_tcga_cn$sampleID]

length(unique(pcawg_tcga_cn$sampleID))
saveRDS(pcawg_tcga_cn, "/data3/wsx_data/Xena/pcawg_tcga_cn.rds")

# Filter samples in target cases
pcawg_tcga_cn <- pcawg_tcga_cn[sampleID %in% target_samps]
length(unique(pcawg_tcga_cn$sampleID))

## Obtain CNV from TCGA samples (SNP array)
tcga_ascat = readRDS("data/tcga_ascat.rds")
tcga_ascat = tcga_ascat[sample %in% target_samps]

## Overlaps region to genes
pcawg_tcga_cn <- pcawg_tcga_cn[, .(chr = paste0("chr", chr), start, end, total_cn, minor_cn, sampleID)]
tcga_snp6_cn <- tcga_ascat[, .(chr = paste0("chr", chr), start, end, total_cn, minor_cn, sampleID = sample)]

tcga_wes_cn <- readRDS("data/tcga_ascat_wes_tumor_cn.rds")$data
colnames(tcga_wes_cn) <- c("chr", "start", "end", "total_cn", "minor_cn", "repeatID", "sampleID")
tcga_wes_cn$chr <- paste0("chr", tcga_wes_cn$chr)
tcga_wes_cn$sampleID <- stringr::str_remove(tcga_wes_cn$sampleID, ".ASCAT.rds")

regions <- modules::use("../lib/regions.R")

overlap_dt_wgs <- regions$collapse_to_genes(pcawg_tcga_cn, ref_file = "data/hg19_genes.bed") # use hg19 version reference
overlap_dt_snp <- regions$collapse_to_genes(tcga_snp6_cn, ref_file = "data/hg38_genes.bed")
overlap_dt_wes <- regions$collapse_to_genes(tcga_wes_cn, ref_file = "data/hg38_genes.bed")

saveRDS(overlap_dt_wgs, file = "/data3/wsx_data/tcga_overlap_gene_wgs.rds")
saveRDS(overlap_dt_snp, file = "/data3/wsx_data/tcga_overlap_gene_snp.rds")
saveRDS(overlap_dt_wes, file = "/data3/wsx_data/tcga_overlap_gene_wes.rds")

compare_different_platform <- function(overlap_file1, overlap_file2, track_repeatID = FALSE, track_len = 15, only_cnv = TRUE) {
  ## When track_repeatID is TRUE, the file with repeatID column should be the second file
  ## track_len specifies how many characters the repeatID and sampleID in common
  library(magrittr)
  
  overlap_dt1 <- readRDS(overlap_file1)
  overlap_dt2 <- readRDS(overlap_file2)
  
  id_share <- if (track_repeatID) {
    intersect(unique(overlap_dt1$sampleID), unique(overlap_dt2$repeatID))
  } else {
    intersect(unique(overlap_dt1$sampleID), unique(overlap_dt2$sampleID))
  }
  
  message("Shared sample number: ", length(id_share))
  
  overlap_dt1 <- overlap_dt1[sampleID %in% id_share]
  overlap_dt2 <- if (track_repeatID) overlap_dt2[repeatID %in% id_share] else overlap_dt2[sampleID %in% id_share]
  
  overlap_dt1[, uid := gsub("(\\..+)", "", gene_id)]
  overlap_dt2[, uid := gsub("(\\..+)", "", gene_id)]
  
  uid_share <- intersect(unique(overlap_dt1$uid), unique(overlap_dt2$uid))
  
  # cosine similarity and rss
  if (track_repeatID) {
    uniq_id <- unique(overlap_dt2$sampleID)
    n_id <- length(uniq_id)
  } else {
    n_id <- length(id_share)
  }
  
  cmp_result <- data.frame(
    sample = vector("character", length = n_id),
    sim_total = vector("numeric", length = n_id),
    rmse_total = vector("numeric", length = n_id),
    sim_minor = vector("numeric", length = n_id),
    rmse_minor = vector("numeric", length = n_id),
    n = vector("numeric", length = n_id)
  )
  
  for (i in seq_len(n_id)) {
    track_id <- if (track_repeatID) uniq_id[i] else id_share[i]
    
    cli::cli_alert("comparing pair for id {track_id}")
    # Drop genes with low intersect ratio
    track_id_dt1 <- if (track_repeatID) substr(track_id, 1, track_len) else track_id
    dt1 <- overlap_dt1[sampleID %in% track_id_dt1 & intersect_ratio > 0.5 & uid %in% uid_share & chr != "chrY", .(uid, total_cn, minor_cn)] %>%
      na.omit() %>%
      unique() %>%
      tibble::column_to_rownames("uid")
    dt2 <- overlap_dt2[sampleID %in% track_id & intersect_ratio > 0.5 & uid %in% uid_share & chr != "chrY", .(uid, total_cn, minor_cn)] %>%
      na.omit() %>%
      unique() %>%
      tibble::column_to_rownames("uid")
    
    id_common <- intersect(rownames(dt1), rownames(dt2))
    dt1 <- dt1[id_common, ]
    dt2 <- dt2[id_common, ]
    
    if (only_cnv) {
      # 避免正常区域检测的影响
      # Focus genes with CNV in at least one data platform
      dt3 <- cbind(dt1, dt2)
      id_common <- rownames(dt3[!(dt3[[1]] == 2 & dt3[[2]] == 1 & dt3[[3]] == 2 & dt3[[4]] == 1), ])
      dt1 <- dt1[id_common, ]
      dt2 <- dt2[id_common, ]
    }
    
    rmse <- function(x, y) {
      n <- length(x)
      sqrt(sum((x - y)^2) / n)
    }
    
    cmp_result$sample[i] <- track_id
    cmp_result$n[i] <- nrow(dt1)
    cmp_result$sim_total[i] <- sigminer::cosine(dt1$total_cn, dt2$total_cn)
    cmp_result$rmse_total[i] <- rmse(dt1$total_cn, dt2$total_cn)
    cmp_result$sim_minor[i] <- sigminer::cosine(dt1$minor_cn, dt2$minor_cn)
    cmp_result$rmse_minor[i] <- rmse(dt1$minor_cn, dt2$minor_cn)
    
    cli::cli_alert(" ({cmp_result$n[i]}, {cmp_result$sim_total[i]}, {cmp_result$rmse_total[i]}, {cmp_result$sim_minor[i]}, {cmp_result$rmse_minor[i]})")
  }
  
  return(cmp_result)
}

overlap_file_wgs <- "/data3/wsx_data//tcga_overlap_gene_wgs.rds"
overlap_file_wes <- "/data3/wsx_data//tcga_overlap_gene_wes.rds"
overlap_file_snp <- "/data3/wsx_data//tcga_overlap_gene_snp.rds"

## WGS and SNP
cmp_wgs_snp <- compare_different_platform(overlap_file_wgs, overlap_file_snp)
# all.equal(cmp_wgs_snp, cmp_result_wgs_snp)
saveRDS(cmp_wgs_snp, file = "data/cnv_comparison_result_between_TCGA_WGS_and_SNParray.rds")
## WGS and WES
cmp_wgs_wes <- compare_different_platform(overlap_file_wgs, overlap_file_wes, track_repeatID = TRUE)
saveRDS(cmp_wgs_wes, file = "data/cnv_comparison_result_between_TCGA_WGS_and_WES.rds")
## SNP and WES
cmp_snp_wes <- compare_different_platform(overlap_file_snp, overlap_file_wes, track_repeatID = TRUE)
saveRDS(cmp_snp_wes, file = "data/cnv_comparison_result_between_TCGA_WES_and_SNParray.rds")

## Include all genes
## WGS and SNP
cmp_wgs_snp_all <- compare_different_platform(overlap_file_wgs, overlap_file_snp, only_cnv = FALSE)
saveRDS(cmp_wgs_snp_all, file = "data/cnv_comparison_result_between_TCGA_WGS_and_SNParray_allgenes.rds")

cmp_wgs_wes_all <- compare_different_platform(overlap_file_wgs, overlap_file_wes, track_repeatID = TRUE, only_cnv = FALSE)
saveRDS(cmp_wgs_wes_all, file = "data/cnv_comparison_result_between_TCGA_WGS_and_WES_allgenes.rds")

cmp_snp_wes_all <- compare_different_platform(overlap_file_snp, overlap_file_wes, track_repeatID = TRUE, only_cnv = FALSE)
saveRDS(cmp_snp_wes_all, file = "data/cnv_comparison_result_between_TCGA_WES_and_SNParray_allgenes.rds")

length(intersect(unique(data_summary_tcga$sample_barcode), unique(tcga_wes_cn$repeatID)))

## Visualization
rm(list = ls())
gc()

comp_df <- data.table::rbindlist(
  list(
    wgs_snp = readRDS("data/cnv_comparison_result_between_TCGA_WGS_and_SNParray.rds"),
    wgs_wes = readRDS("data/cnv_comparison_result_between_TCGA_WGS_and_WES.rds"),
    snp_wes = readRDS("data/cnv_comparison_result_between_TCGA_WES_and_SNParray.rds")
  ),
  idcol = "platform_pair"
)

comp_df_all <- data.table::rbindlist(
  list(
    wgs_snp = readRDS("data/cnv_comparison_result_between_TCGA_WGS_and_SNParray_allgenes.rds"),
    wgs_wes = readRDS("data/cnv_comparison_result_between_TCGA_WGS_and_WES_allgenes.rds"),
    snp_wes = readRDS("data/cnv_comparison_result_between_TCGA_WES_and_SNParray_allgenes.rds")
  ),
  idcol = "platform_pair"
)

plot_comparison <- function(comp_df) {
  count_df <- comp_df %>% dplyr::count(platform_pair)
  pair_levels <- c("wgs_snp", "wgs_wes", "snp_wes")
  names(pair_levels) <- c(
    paste0("WGS and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_snp") %>% dplyr::pull(n), ")")),
    paste0("WGS and WES", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_wes") %>% dplyr::pull(n), ")")),
    paste0("WES and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "snp_wes") %>% dplyr::pull(n), ")"))
  )
  
  
  df <- comp_df %>%
    tidyr::pivot_longer(cols = c("sim_total", "rmse_total", "sim_minor", "rmse_minor"), names_to = "measure", values_to = "score") %>%
    dplyr::mutate(
      platform_pair = forcats::fct_recode(platform_pair, !!!pair_levels),
      measure = forcats::fct_recode(measure,
                                    `Total CN similarity` = "sim_total",
                                    `Total CN RMSE` = "rmse_total",
                                    `Minor CN similarity` = "sim_minor",
                                    `Minor CN RMSE` = "rmse_minor"
      )
    )
  
  library(rstatix)
  stat.test <- df %>%
    group_by(measure) %>%
    wilcox_test(score ~ platform_pair, p.adjust.method = "bonferroni") %>%
    add_xy_position(scales = "free_y")
  stat.test$y.position[7:12] <- rep(c(1.01, 1.05, 1.1), 2)
  
  library(ggplot2)
  library(ggpubr)
  ggplot(df, aes(platform_pair, score)) +
    geom_violin(aes(fill = platform_pair)) + # cannot set fill in ggplot() for using stat_pvalue_manual
    facet_wrap(~measure, scales = "free", nrow = 2) +
    cowplot::theme_cowplot() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "none"
    ) +
    stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, tip.length = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # Add 10% spaces between the p-value labels and the plot border
}

plot_comparison(comp_df)
plot_comparison(comp_df_all)

ggsave("plots/platform_comparison_for_genes_with_cnv.pdf", plot = plot_comparison(comp_df), width = 6, height = 8)
ggsave("plots/platform_comparison_for_all_genes.pdf", plot = plot_comparison(comp_df_all), width = 6, height = 8)


plot_comparison2 <- function(comp_df) {
  count_df <- comp_df %>% dplyr::count(platform_pair)
  pair_levels <- c("wgs_snp", "wgs_wes", "snp_wes")
  names(pair_levels) <- c(
    paste0("WGS and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_snp") %>% dplyr::pull(n), ")")),
    paste0("WGS and WES", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_wes") %>% dplyr::pull(n), ")")),
    paste0("WES and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "snp_wes") %>% dplyr::pull(n), ")"))
  )
  
  
  df <- comp_df %>%
    tidyr::pivot_longer(cols = c("sim_total", "rmse_total"), names_to = "measure", values_to = "score") %>%
    dplyr::mutate(
      platform_pair = forcats::fct_recode(platform_pair, !!!pair_levels),
      measure = forcats::fct_recode(measure,
                                    `CN similarity` = "sim_total",
                                    `CN RMSE` = "rmse_total",
      )
    )
  
  library(ggplot2)
  library(ggpubr)
  ggplot(df, aes(platform_pair, score)) +
    geom_violin(aes(fill = platform_pair)) + # cannot set fill in ggplot() for using stat_pvalue_manual
    facet_wrap(~measure, scales = "free", nrow = 1) +
    cowplot::theme_cowplot() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "none"
    )
}

p1 = plot_comparison2(comp_df %>% dplyr::select(-dplyr::ends_with("minor")))
p2 = plot_comparison2(comp_df_all %>% dplyr::select(-dplyr::ends_with("minor")))

ggsave("plots/platform_comparison_for_genes_with_cnv_only_total.pdf", plot = p1, width = 6, height = 4)
ggsave("plots/platform_comparison_for_all_genes_only_total.pdf", plot = p2, width = 6, height = 4)

