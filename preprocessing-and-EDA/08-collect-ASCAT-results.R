setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(dplyr)

# The paths below points to HPC in old repo, which can not be used now
# Directly use the output file if want to reproduce result.
files <- list.files("~/proj/ecDNA/raw/ascat/result", pattern = ".rds", full.names = TRUE)
# Add samples with circular amplicon records less than 4
files2 = list.files("~/proj/ecDNA/data/outval_result", pattern = ".rds", full.names = TRUE)
# Add samples without circular amplicon
files3 = list.files("~/proj/ecDNA/data/nonEC_result", pattern = ".rds", full.names = TRUE)
files = c(files, files2, files3) # 507 paired results

# The ASCAT result files for WES input have been archived in 51 server
# Path: /data3/wsx_data/tcga_wes_ascat/

read_ascat_cn <- function(x) {
  message("reading ", x)
  x2 <- x
  x <- readRDS(x)
  
  if (object.size(x) < 2000) {
    warning(x2, " was failed in ASCAT calling.", immediate. = TRUE)
    return(invisible(NULL))
  }
  x <- x[c("segments", "aberrantcellfraction", "ploidy")]
  names(x) <- c("data", "purity", "ploidy")
  colnames(x$data) <- c("sample", "chromosome", "start", "end", "major_cn", "minor_cn")
  x$data <- x$data %>%
    dplyr::mutate(total_cn = .data$major_cn + .data$minor_cn) %>%
    dplyr::select(c("chromosome", "start", "end", "total_cn", "minor_cn", "sample"))
  x$data$source <- basename(x2) # track the source file
  return(x)
}

read_ascat_cn_list <- function(x_list) {
  message("reading file list:")
  x_list <- purrr::transpose(purrr::compact(lapply(x_list, read_ascat_cn)))
  y_list <- list()
  message("transforming data")
  y_list$data <- data.table::rbindlist(x_list$data)
  y_list$purity <- purrr::reduce(x_list$purity, c)
  y_list$ploidy <- purrr::reduce(x_list$ploidy, c)
  message("done")
  return(y_list)
}

tcga_wes_cn <- read_ascat_cn_list(files)

saveRDS(tcga_wes_cn, file = "data/tcga_ascat_wes_tumor_cn.rds")

# 500 successful pairs
