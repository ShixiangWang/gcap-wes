library(data.table)

# pc3 = fread("/data3/wsx_data/cnvkit_result/tumor_vs_sampleB/pc3-1.call.cns2")
# snu16 = fread("/data3/wsx_data/cnvkit_result/tumor_vs_sampleB/SNU16-1.call.cns2")
# 
# dt = rbind(
#   # https://en.wikipedia.org/wiki/PC3
#   #  62-year-old Caucasian male
#   pc3[chromosome %in% paste0("chr", 1:22),
#       .(chromosome, start, end, total_cn = cn, minor_cn = NA, sample = "pc3_WES",
#         purity = 1)],
#   # https://www.atcc.org/products/crl-5974
#   # 33-year-old, female, Asian, stomach cancer patient
#   snu16[chromosome %in% paste0("chr", 1:22),
#         .(chromosome, start, end, total_cn = cn, minor_cn = NA, sample = "snu16_WES",
#           purity = 1)]
# )

fl = list.files("~/gcap-analysis/manuscript/data/segs/", full.names = TRUE)

dt = rbindlist(
  lapply(fl, function(x) {
    dt = data.table::fread(x)
    dt = dt[chromosome %in% paste0("chr", 1:22),
       .(chromosome, start, end, total_cn = cn, minor_cn = NA, 
         sample = sub("^([^.]+).*", "\\1", basename(x)),
         purity = 1)]
    dt[, sample := ifelse(sample %in% c("pc3-1", "SNU16-1"), paste0(sample, "_WES"), paste0(sample, "_WGS"))]
    dt
  })
)

#dt = rbind(dt, dt_fl)

library(gcap)
?gcap.ASCNworkflow()

rv = gcap.ASCNworkflow(
  dt,
  genome_build = "hg38",
  model = readRDS("../modeling/data/xgb_v3/XGB_NF11.rds"),
  tightness = 1,
  gap_cn = 4,
  overlap = 0.5,
  outdir = "/data3/wsx_data/raw_cell_line/gcap_result_cnvkit",
  result_file_prefix = "cellline"
)

saveRDS(rv, file = "data/cellline.rds")

rv$convertGeneID()
View(rv$data)
rv$sample_summary

# https://pubmed.ncbi.nlm.nih.gov/25026282/

data = readRDS("/data3/wsx_data/raw_cell_line/gcap_result_cnvkit/cellline_prediction_result.rds")
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
data[, gene_name := convert_hm_genes(gene_id)]
data2 = data[, list(sample, band, gene_id, gene_name, total_cn, purity, ploidy, AScore, cna_burden,
            freq_Linear, freq_BFB, freq_HR, freq_Circular, prob, background_cn, gene_class)]
data2
data2 = data2[order(gene_class, -total_cn, -prob)]
data2

fwrite(data2, file = "data/cell_line_gcap.tsv", sep = "\t")
