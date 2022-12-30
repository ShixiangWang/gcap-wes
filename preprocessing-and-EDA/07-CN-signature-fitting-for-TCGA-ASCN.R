# Fitting copy number signatures with 19 reference signatures
# This is outdated, as Nature paper by Steel have illustrated this result,
# no need to modify/update further.

setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(sigminer)
library(dplyr)

tcga <- modules::use("https://biosisyphus.github.io/Rlib/tcga.R")

dt <- readRDS("data/tcga_ascat_snp6_array_result.rds")
dt$GDC_Aliquot <- NULL

# a sample may have multiple aliquots
# filter aliquot barcodes
keep_samps <- unique(dt$Sample)
keep_samps <- keep_samps[as.integer(substr(keep_samps, 14, 15)) < 10]
keep_samps <- tcga$filterReplicates(keep_samps)
# Check
which(duplicated(substr(keep_samps, 1, 15)))
keep_samps[which(duplicated(substr(keep_samps, 1, 15)))]

dt <- dt[dt$Sample %in% keep_samps, ]
# 5401 tumor samples

table(dt$Chromosome)

# Remove sex chrs
dt <- dt[!dt$Chromosome %in% c("chrX", "chrY"), ]

colnames(dt)[6] <- "minor_cn"
cn <- read_copynumber(
  dt,
  seg_cols = c("Chromosome", "Start", "End", "Copy_Number"),
  samp_col = "Sample",
  max_copynumber = 1000L,
  min_segnum = 22L,
  join_adj_seg = FALSE,
  genome_build = "hg38",
  add_loh = TRUE,
  loh_min_len = 1e3
)

cn@summary.per.sample

cn_tally <- sig_tally(cn, method = "S")

# Handling overlapping samples
load("data/pancan_amplicon_list_and_summary.RData")

tcga_samples <- rownames(cn_tally$all_matrices$CN_48)
names(tcga_samples) <- substr(tcga_samples, 1, 15)

overlap_samples <- intersect(names(tcga_samples), data_summary_tcga$sample_barcode)

overlap_mat <- t(cn_tally$all_matrices$CN_48[tcga_samples[overlap_samples], ])
colnames(overlap_mat) <- substr(colnames(overlap_mat), 1, 15)

fit_result <- sig_fit(overlap_mat, sig_db = "CNS_TCGA", sig_index = "ALL", return_class = "data.table")
data <- dplyr::inner_join(data_summary_tcga, fit_result, by = c("sample_barcode" = "sample"))

dir.create("data")
saveRDS(fit_result, file = "data/tcga_overlap_samples_cns_activity.rds")

data2 <- data
# Add some other extra measures
data <- data %>%
  dplyr::left_join(
    cn@summary.per.sample[, .(sample, n_loh, cna_burden)] %>% dplyr::mutate(sample = substr(sample, 1, 15)),
    by = c("sample_barcode" = "sample")
  ) %>%
  dplyr::left_join(
    get_cn_ploidy(cn) %>% dplyr::mutate(sample = substr(sample, 1, 15)),
    by = c("sample_barcode" = "sample")
  ) %>%
  dplyr::left_join(
    get_Aneuploidy_score(cn, rm_black_arms = TRUE) %>% dplyr::mutate(sample = substr(sample, 1, 15)) %>% dplyr::select(sample, AScore),
    by = c("sample_barcode" = "sample")
  )

plot(log(data$CN8 + 1), log(data$Circular + 1))
cor(data$CN8, data$Circular)

summary(lm(paste("Circular", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data))
summary(lm(paste("BFB", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data))
summary(lm(paste("Linear", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data))
summary(lm(paste("`Heavily-rearranged`", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data))

summary(lm(Circular ~ CN7 + CN8 + CN9 + CN13, data = data))

show_cor(data, x_vars = colnames(data)[2:5], y_vars = colnames(data)[6:24], cor_method = "pearson")

# Get model summary table
library(gtsummary)

tbl_regression(lm(paste("Circular", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data)) %>%
  bold_labels() %>%
  bold_p() %>%
  add_glance_table(
    # label = list(sigma ~ "\U03C3"),
    include = c(r.squared)
  )

tbl_regression(
  lm(paste("Circular", paste(c(paste0("CN", 1:19), c("n_loh", "cna_burden", "ploidy", "AScore")), collapse = "+"), sep = "~"), data = data)
) %>%
  bold_labels() %>%
  bold_p() %>%
  add_glance_table(
    # label = list(sigma ~ "\U03C3"),
    include = c(r.squared)
  )

tbl_regression(
  lm(paste("Circular", paste(c("n_loh", "cna_burden", "ploidy", "AScore"), collapse = "+"), sep = "~"), data = data)
) %>%
  bold_labels() %>%
  bold_p() %>%
  add_glance_table(
    # label = list(sigma ~ "\U03C3"),
    include = c(r.squared)
  )

# # Bootstrap fitting
# fit_result_bt <- sig_fit_bootstrap_batch(
#   overlap_mat,
#   sig_db = "CNS_TCGA",
#   sig_index = "ALL",
#   n = 100,
#   use_parallel = 4,
#   job_id = "tcga_overlap_samples_cns",
#   result_dir = "data/tcga_cns_activity"
# )
# 
# # Use median activity as robust estimation
# fit_result_md <- data.table::dcast(fit_result_bt$expo[type != "optimal"][, .(exposure = median(exposure)), by = .(sample, sig)], sample ~ sig, value.var = "exposure")
# 
# save(fit_result_bt, fit_result_md, file = "data/tcga_overlap_samples_cns_bt_activity.RData")
# 
# 
# data2 <- dplyr::inner_join(data_summary_tcga, fit_result_md, by = c("sample_barcode" = "sample"))
# 
# summary(lm(paste("Circular", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data2))
# summary(lm(paste("BFB", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data2))
# summary(lm(paste("Linear", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data2))
# summary(lm(paste("`Heavily-rearranged`", paste(paste0("CN", 1:19), collapse = "+"), sep = "~"), data = data2))
# 
# summary(lm(Circular ~ CN7 + CN8 + CN9 + CN13, data = data))
# 
# show_cor(data2, x_vars = colnames(data2)[2:5], y_vars = colnames(data2)[6:24], cor_method = "pearson")

# The results from optimal activity are very consistent with median bootstrapping activity
# So no need to run bootstrapping analysis for copy number signatures based on my experience.


# # Combine CNS
#
# fit <- lm(Circular ~ CN7 + CN8 + CN9 + CN13, data = data)
# coef(fit)
# coef_fit <- coef(fit)
#
# data$score <- coef_fit[1] + coef_fit[2] * data$CN7 + coef_fit[3] * data$CN8 +coef_fit[4] * data$CN9 + coef_fit[5] * data$CN13
#
# plot(Circular ~ score, data = data)
#
# cor(data$Circular, data$CN8)
# cor(data$Circular, data$score)
