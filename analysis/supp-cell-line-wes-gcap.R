# library(dplyr)
# 
# info = dplyr::tribble(
#   ~id, ~type, ~path,
#   # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP292966&o=acc_s%3Aa
#   "SampleB_SRR13084972", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084972.bam",
#   "SampleB_SRR13084973", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084973.bam",
#   "SampleB_SRR13084974", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084974.bam",
#   "SampleB_SRR13084975", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084975.bam",
#   "SampleB_SRR13084976", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084976.bam",
#   "SampleB_SRR13084977", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084977.bam",
#   # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP078289&o=acc_s%3Aa&s=SRR3929836
#   "SNU16_SRR3929836", "tumor", "/data3/wsx_data/raw_cell_line/fq/bam/SRR3929836.bam",
#   # CCLE 女性 胃癌 33岁 亚裔
#   # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP186687&search=SRR8657217&o=acc_s%3Aa
#   "SNU16_SRR8657217", "tumor", "/data3/wsx_data/raw_cell_line/fq/bam/SRR8657217.bam",
#   # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP026538&o=acc_s%3Aa
#   "HCCBL_SRR925766", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR925766.bam",
#   "HCCBL_SRR925779", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR925779.bam",
#   # Own data
#   "pc3-1", "tumor", "/data3/wsx_data/raw_cell_line/own_data/bam/pc3-1.bam",
#   "SNU16-1", "tumor", "/data3/wsx_data/raw_cell_line/own_data/bam/SNU16-1.bam"
# )
# 
# all_info = expand.grid(
#   tumor = info %>% filter(type == "tumor") %>% pull(id),
#   normal = info %>% filter(type == "normal") %>% pull(id)
# )
# 
# data = left_join(all_info, info %>% filter(type == "tumor"), by = c("tumor"="id")) %>%
#   rename(tumor_path = path) %>%
#   left_join(info %>% filter(type == "normal"), by = c("normal"="id")) %>%
#   select(-starts_with("type")) %>%
#   mutate(id = paste(tumor, normal, sep = "__")) %>% 
#   rename(normal_path = path)
# 
# data.table::fwrite(data, file = "data/cell_line_pair_info.csv")

df = data.table::fread("data/cell_line_pair_info.csv")
df

library(gcap)
# cd
# for i in $(ls fq/bam/*.bam own_data/bam/*.bam); do samtools index -@ 20 $i; done

outdir = "/data3/wsx_data/raw_cell_line/gcap_result"


tfile <- df$tumor_path
nfile <- df$normal_path
tn <- df$tumor
nn <- df$normal
id <- df$id

gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile,
  tumourname = tn, normalname = nn, jobname = id,
  outdir = outdir,
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", 
  g1000allelesprefix = file.path(
    "~/data/1000G_loci_hg38/",
    "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
  ), 
  g1000lociprefix = file.path("~/data/1000G_loci_hg38/",
                              "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
  ),
  GCcontentfile = "~/data/GC_correction_hg38.txt",
  replictimingfile = "~/data/RT_correction_hg38.txt",
  skip_finished_ASCAT = TRUE,
  result_file_prefix = "cell_line",
  genome_build = "hg38",
  # Downsample model have better recall
  model = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample.rds"),
  tightness = 0,
  gap_cn = 3
)

gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile,
  tumourname = tn, normalname = nn, jobname = id,
  outdir = outdir,
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", 
  g1000allelesprefix = file.path(
    "~/data/1000G_loci_hg38/",
    "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
  ), 
  g1000lociprefix = file.path("~/data/1000G_loci_hg38/",
                              "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
  ),
  GCcontentfile = "~/data/GC_correction_hg38.txt",
  replictimingfile = "~/data/RT_correction_hg38.txt",
  skip_finished_ASCAT = TRUE,
  result_file_prefix = "cell_line_overlap50",
  genome_build = "hg38",
  # Downsample model have better recall
  model = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample.rds"),
  tightness = 0,
  gap_cn = 3,
  overlap = 0.5
)
