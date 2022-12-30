# This is outdated, as I found the data download had something wrong in ID mapping https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/523
# In future, we directly use ASCAT results from Steele, Nature, 2022.

# UPDATE: this bug shall be fixed by the author in latest TCGAbiolinks version.
setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(TCGAbiolinks)
library(stringr)

all_projs <- TCGAbiolinks:::getGDCprojects()$project_id

query <- GDCquery(
  project = all_projs %>% stringr::str_subset("TCGA"),
  data.category = "Copy Number Variation",
  data.type = "Allele-specific Copy Number Segment",
  platform = "Affymetrix SNP 6.0",
  workflow.type = "ASCAT2",
  experimental.strategy = "Genotyping Array"
)

View(getResults(query))

# Some data may be obtained from different tumor-normal pair for a same case
# Remove duplicated data
sum(duplicated(query$results[[1]]$cases))

query$results[[1]]$cases <- make.unique(query$results[[1]]$cases)

GDCdownload(query,
  method = "api",
  # .per.chunk = 10,
  directory = "/data3/wsx_data/GDCdata"
)
data <- GDCprepare(query, directory = "/data3/wsx_data/GDCdata")

saveRDS(data, file = "/data3/wsx_data/GDCdata/gdc_all_tcga_ascat_array.rds")
