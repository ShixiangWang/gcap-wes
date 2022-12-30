get_highlights = function(fCNA, group = "type",
                          gene_freq = 16, grp_freq = 3, 
                          #min_mean_cn = 16,
                          only_oncogenes = FALSE,
                          return_dt = FALSE) {
  library(data.table)
  library(IDConverter)
  
  data = merge(fCNA$data, fCNA$sample_summary[, unique(c("sample", group)), with = FALSE], by = "sample")
  data
  data[, band2 := sub("(\\..+)", "", band)]
  setnames(data, group, "group")
  dt_band = data[, list(N = length(unique(band[gene_class == "circular"]))), by = list(group, band2)]
  
  dt_band_freq = dt_band[, list(freq = sum(N)), by = list(band2)][order(band2)]
  dt_grp_freq = dt_band[N > 0, list(freq = length(group)), by = list(band2)][order(band2)]
  
  band2_sel = intersect(
    dt_band_freq[freq >= gene_freq]$band2,
    dt_grp_freq[freq >= grp_freq]$band2
  )
  
  dt_gene = data[band2 %in% band2_sel & gene_class == "circular", 
                 list(N = .N, CN = mean(total_cn, na.rm = TRUE)),
                 by = list(band2, gene_id)][order(band2, -N)]
  dt_gene
  options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
  dt_gene[, gene_name := convert_hm_genes(gene_id)]
  dt_gene
  if (only_oncogenes) {
    oncogenes = unique(gcap::oncogenes$OncogeneName)
    dt_gene = dt_gene[gene_name %in% oncogenes]
  }
  
  #dt_gene2 = dt_gene2[N >= gene_freq & CN >= min_mean_cn]
  dt_gene2 = dt_gene[N >= gene_freq]
  if (return_dt) {
    return(dt_gene2)
  }
  
  dt_gene2 = dt_gene2[!is.na(gene_name), head(.SD, 3), by = list(band2)]
  dt_gene2 = merge(dt_gene2, dt_grp_freq, by = "band2", all.x = TRUE)
  dt_anno = dt_gene2[ , list(id = paste0(paste(gene_name, collapse = ","), "(", freq[1], ")")), by = list(band2)]
  rv = ifelse(grepl(",", dt_anno$id), sub(",", " ,", dt_anno$id), sub("\\(", " \\(", dt_anno$id))
  rv
}

get_highlights_onco = function(fCNA, group = "type", gene_freq = 16, grp_freq = 3,
                               return_simple = FALSE) {
  library(data.table)
  library(IDConverter)
  # https://www.nature.com/articles/s41568-020-0290-x#additional-information
  # https://www.intogen.org/download
  # intogen = data.table::fread("data/Compendium_Cancer_Genes.tsv")
  # 223 act genes
  # oncogenes = unique(intogen[ROLE == "Act"]$SYMBOL)
  oncogenes = unique(gcap::oncogenes$OncogeneName)
  
  data = merge(fCNA$data, fCNA$sample_summary[, unique(c("sample", group)), with = FALSE], by = "sample")
  data
  setnames(data, group, "group")
  options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
  data[, gene_name := convert_hm_genes(gene_id)]
  
  dt_gene_freq = data[!is.na(gene_name) & gene_name %in% oncogenes & gene_class == "circular", .N, by = list(gene_name)][order(-N)]
  dt_grp_freq = data[!is.na(gene_name) &gene_name %in% oncogenes & gene_class == "circular", .N, by = list(group, gene_name)]
  dt_grp_freq = dt_grp_freq[N > 0, list(N = length(group)), by = list(gene_name)][order(-N)]
  
  
  rv = dt_grp_freq[gene_name %in% dt_gene_freq[N >= gene_freq]$gene_name][N >= grp_freq]
  if (return_simple) {
    rv$gene_name
  } else {
    rv[order(-N), list(paste0(gene_name, " (", N, ")"))][[1]]
  }
}


