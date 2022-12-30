library(dplyr)

info = dplyr::tribble(
  ~id, ~type, ~path,
  # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP292966&o=acc_s%3Aa
  # https://europepmc.org/article/MED/33863366
  # had a typical depth of 235× with an average depth range from 151× to 402× in their targeted regions
  # 这两个参考好像有问题，是不是混了癌细胞？
  #"SampleB_SRR13084972", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084972.bam",
  #"SampleB_SRR13084973", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR13084973.bam",
  # https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP078289&o=acc_s%3Aa&s=SRR3929836
  "SNU16_SRR3929836", "tumor", "/data3/wsx_data/raw_cell_line/fq/bam/SRR3929836.bam",
  # 女性 胃癌 33岁 亚裔
  # Own data
  "pc3-1", "tumor", "/data3/wsx_data/raw_cell_line/own_data/bam/pc3-1.bam",
  "SNU16-1", "tumor", "/data3/wsx_data/raw_cell_line/own_data/bam/SNU16-1.bam",
  # Another normal
  "HCCBL_SRR925779", "normal", "/data3/wsx_data/raw_cell_line/fq/bam/SRR925779.bam"
)

all_info = expand.grid(
  tumor = info %>% filter(type == "tumor") %>% pull(id),
  normal = info %>% filter(type == "normal") %>% pull(id)
)

data = left_join(all_info, info %>% filter(type == "tumor"), by = c("tumor"="id")) %>%
  rename(tumor_path = path) %>%
  left_join(info %>% filter(type == "normal"), by = c("normal"="id")) %>%
  select(-starts_with("type")) %>%
  mutate(id = paste(tumor, normal, sep = "__")) %>%
  rename(normal_path = path)

data.table::fwrite(data, file = "data/cell_line_pair_info2.csv")

df = data.table::fread("data/cell_line_pair_info2.csv")
df

library(gcap)
# cd
# for i in $(ls fq/bam/*.bam own_data/bam/*.bam); do samtools index -@ 20 $i; done

outdir = "/data3/wsx_data/raw_cell_line/gcap_result2"

df = df %>% dplyr::arrange(tumor)
df

for (i in seq_len(nrow(df))) {
  tfile <- df$tumor_path[i]
  nfile <- df$normal_path[i]
  tn <- df$tumor[i]
  nn <- df$normal[i]
  id <- df$id[i]
  
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
    result_file_prefix = id,
    genome_build = "hg38",
    model = readRDS("../modeling/data/xgb_v3/XGB_NF11.rds"),
    tightness = 1,
    gap_cn = 4
  )
}

# x[POS > 127735434 & POS < 127742951]
#xx = readRDS("/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1__SampleB_SRR13084972.ASCAT.rds")
#xx = readRDS("/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1__SampleB_SRR13084973.ASCAT.rds")
xx = readRDS("/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1__HCCBL_SRR925779.ASCAT.rds")
xx = readRDS("/data3/wsx_data/raw_cell_line/gcap_result2/SNU16-1__HCCBL_SRR925779.ASCAT.rds")
xx1 = as.data.table(xx$segments)
xx1[chr == 8 & !endpos < 127700000 & startpos < 127800000]

# library(IDConverter)
# options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
# xx[, gene_id := convert_hm_genes(gene_id)]

library(ASCAT)
GCcontentfile = "~/data/GC_correction_hg38.txt"
replictimingfile = "~/data/RT_correction_hg38.txt"

ascat.bc <- ascat.loadData(
  "/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1_tumourLogR.txt",
  "/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1_tumourBAF.txt",
  "/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1_normalLogR.txt",
  "/data3/wsx_data/raw_cell_line/gcap_result2/pc3-1_normalBAF.txt"
)
ascat.bc$Tumor_LogR[startsWith(rownames(ascat.bc$Tumor_LogR), "8_12395"), , drop = FALSE]
ascat.bc$Tumor_LogR[startsWith(rownames(ascat.bc$Tumor_LogR), "8_127738"), , drop = FALSE] # MYC
ascat.bc$Tumor_BAF[startsWith(rownames(ascat.bc$Tumor_BAF), "8_12395"), , drop = FALSE]
ascat.bc$Tumor_BAF[startsWith(rownames(ascat.bc$Tumor_BAF), "8_127738"), , drop = FALSE] # MYC

#ascat.bc <- ascat.GCcorrect(ascat.bc, GCcontentfile, replictimingfile)
ascat.bc.a <- ascat.aspcf(ascat.bc, ascat.gg = NULL, penalty = 10)
ascat.output <- ascat.runAscat(ascat.bc.a, gamma = 1L, pdfPlot = FALSE)

ascat.bc.a$SNPpos["8_127740678", ]


ascat.output$nA["8_127740678", ]
ascat.output$nB["8_127740678", ]
as.data.table(ascat.output$segments_raw)[chr == 8 & startpos > 120000000 & endpos < 130000000]
as.data.table(ascat.output$segments)[chr == 8 & startpos > 120000000 & endpos < 130000000]


ascat.bc.a$SNPpos[startsWith(rownames(ascat.bc.a$SNPpos), "8_12395"), ]


ascat.bc.a$Tumor_LogR[startsWith(rownames(ascat.bc.a$Tumor_LogR), "8_12395"), , drop = FALSE]
ascat.bc.a$Tumor_LogR[startsWith(rownames(ascat.bc.a$Tumor_LogR), "8_127738"), , drop = FALSE] # MYC

ascat.bc.a$Tumor_LogR_segmented[startsWith(rownames(ascat.bc.a$Tumor_LogR_segmented), "8_12395"), , drop = FALSE]
ascat.bc.a$Tumor_LogR_segmented[startsWith(rownames(ascat.bc.a$Tumor_LogR_segmented), "8_127738"), , drop = FALSE] # MYC


ascat.bc.a$Tumor_BAF[startsWith(rownames(ascat.bc.a$Tumor_BAF), "8_12395"), , drop = FALSE]
ascat.bc.a$Tumor_BAF[startsWith(rownames(ascat.bc.a$Tumor_BAF), "8_127738"), , drop = FALSE] # MYC

ascat.bc.a$Tumor_BAF_segmented[[1]][startsWith(rownames(ascat.bc.a$Tumor_BAF_segmented[[1]]), "8_12395"), ]
ascat.bc.a$Tumor_BAF_segmented[[1]][startsWith(rownames(ascat.bc.a$Tumor_BAF_segmented[[1]]), "8_12773"), ] # MYC
