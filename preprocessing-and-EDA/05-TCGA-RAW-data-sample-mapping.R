setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(tidyverse)

df <- read_tsv("../gdc-manifest/gdc_sample_sheet.2021-07-23.tsv")

table(df$`Sample Type`)

df <- df %>%
  select(-c(3:4)) %>%
  setNames(c("file_id", "file_name", "project", "case", "sample", "type")) %>%
  mutate(sample = substr(sample, 1, 15), type2 = ifelse(grepl("Normal", type), "normal", "tumor"), center = ifelse(startsWith(
    file_name,
    "C"
  ), stringr::str_split(file_name, pattern = "\\.") %>%
    purrr::map_chr(~ .x[[1]]), "None"))

test <- dplyr::filter(df, case == "TCGA-BL-A0C8")
test %>%
  pivot_wider(id_cols = "case", names_from = "type2", values_from = "sample", values_fn = list) -> z2
unique(expand.grid(z2$tumor[[1]], z2$normal[[1]])) %>%
  setNames(c("tumor", "normal"))

build_df <- function(df) {
  r <- pivot_wider(df, id_cols = "case", names_from = "type2", values_from = "sample", values_fn = list)
  if (!"tumor" %in% colnames(r)) {
    r <- dplyr::tibble(tumor = rep(NA_character_, length(unique(r$normal[[1]]))), normal = unique(r$normal[[1]]))
  } else if (!"normal" %in% colnames(r)) {
    r <- dplyr::tibble(tumor = unique(r$tumor[[1]]), normal = rep(NA_character_, length(unique(r$tumor[[1]]))))
  } else {
    r <- unique(expand.grid(r$tumor[[1]], r$normal[[1]])) %>%
      setNames(c("tumor", "normal"))
  }

  nrows <- cumprod(table(df$type2))[2]
  nr <- nrow(r)
  idxT <- df$sample %in% r$tumor
  idxN <- df$sample %in% r$normal

  if (!is.na(nrows) && nrows > 1) {
    bind_cols(r[rep(seq_len(nr), nrows / nr), ], df[idxT, c("file_id", "file_name", "center")][rep(seq_len(sum(idxT)), nrows / sum(idxT)), ] %>%
      setNames(c("file_id_tumor", "file_name_tumor", "center")), df[idxN, c("file_id", "file_name")][rep(
      seq_len(sum(idxN)),
      nrows / sum(idxN)
    ), ] %>%
      setNames(c("file_id_normal", "file_name_normal"))) %>%
      as_tibble()
  } else {
    bind_cols(r, if (any(idxT)) {
      df[idxT, c("file_id", "file_name", "center")] %>%
        setNames(c("file_id_tumor", "file_name_tumor", "center"))
    } else {
      NULL
    }, if (any(idxN)) {
      df[idxN, c("file_id", "file_name")] %>%
        setNames(c("file_id_normal", "file_name_normal"))
    }) %>%
      as_tibble()
  }
}

# debug(build_df) build_df(z[[2]])

tcga_sample_info <- df %>%
  group_split(center, case) %>%
  map_df(build_df) %>%
  filter(!is.na(tumor)) %>%
  left_join(df %>%
    filter(type2 == "tumor") %>%
    select(-file_id, -file_name, -type2) %>%
    unique(), by = c(tumor = "sample", center = "center")) %>%
  mutate(pair_id_uniq = make.unique(tumor)) %>%
  select(case, pair_id_uniq, tumor, file_name_tumor, normal, file_name_normal, type, project, center, file_id_tumor, file_id_normal)


# There are two cases have only normal samples
setdiff(unique(df$case), unique(tcga_sample_info$case))

length(unique(tcga_sample_info$case))

write_csv(tcga_sample_info, file = "data/tcga_wes_pair_info.csv")


# ---------------------------------------------------------------------------------
# Tumor normal pairing for tumor with 1-3 ecDNA regions

library(tidyverse)

df <- read_tsv("../gdc-manifest/gdc_sample_sheet.2021-10-15.tsv")

table(df$`Sample Type`)

df <- df %>%
  select(-c(3:4)) %>%
  setNames(c("file_id", "file_name", "project", "case", "sample", "type")) %>%
  mutate(type2 = ifelse(grepl("Normal", type), "normal", "tumor"))

table(df$type)
table(df$type2)

# sum(grepl("hg19", df$file_name))
# Although many file names have been laebelled as hg19
# but the GDC portal indicates their reference genome is hg38
# e.g., https://portal.gdc.cancer.gov/files/04a38123-6191-4848-8b0e-3a6b1749e2a7

build_df <- function(df) {
  message("handling sample ", paste(df$sample, collapse = ","))
  
  if (any(grepl("Solid", df$type)) & any(grepl("Blood", df$type))) {
    # Use blood normal as reference
    df <- df %>% filter(!grepl("Solid", df$type))
  }
  
  r <- pivot_wider(df,
                   id_cols = "case",
                   names_from = "type2", values_from = "sample", values_fn = list
  )
  if (!"tumor" %in% colnames(r)) {
    r <- dplyr::tibble(
      tumor = rep(NA_character_, length(unique(r$normal[[1]]))),
      normal = unique(r$normal[[1]])
    )
  } else if (!"normal" %in% colnames(r)) {
    r <- dplyr::tibble(
      tumor = unique(r$tumor[[1]]),
      normal = rep(NA_character_, length(unique(r$tumor[[1]])))
    )
  } else {
    r <- expand.grid(r$tumor[[1]], r$normal[[1]]) %>%
      setNames(c("tumor", "normal"))
  }
  
  nrows <- cumprod(table(df$type2))[2]
  nr <- nrow(r)
  idxT <- df$sample %in% r$tumor
  idxN <- df$sample %in% r$normal
  
  if (!is.na(nrows) && nrows > 1) {
    bind_cols(
      r[rep(seq_len(nr), nrows / nr), ],
      df[idxT, c("file_id", "file_name")][rep(seq_len(sum(idxT)), nrows / sum(idxT)), ] %>% # nolint
        setNames(c("file_id_tumor", "file_name_tumor")),
      df[idxN, c("file_id", "file_name")][rep(seq_len(sum(idxN)), nrows / sum(idxN)), ] %>% # nolint
        setNames(c("file_id_normal", "file_name_normal"))
    ) %>%
      as_tibble()
  } else {
    bind_cols(
      r,
      if (any(idxT)) {
        df[idxT, c("file_id", "file_name")] %>%
          setNames(c("file_id_tumor", "file_name_tumor"))
      } else {
        NULL
      },
      if (any(idxN)) {
        df[idxN, c("file_id", "file_name")] %>%
          setNames(c("file_id_normal", "file_name_normal"))
      }
    ) %>%
      as_tibble()
  }
}

pair_info <- df %>%
  group_split(case) %>%
  map_df(build_df) %>%
  filter(!is.na(tumor)) %>%
  left_join(df %>%
              filter(type2 == "tumor") %>%
              select(-file_id, -file_name, -type2) %>%
              unique(),
            by = c("tumor" = "sample")
  ) %>%
  mutate(pair_id_uniq = make.unique(tumor)) %>%
  select(
    case, pair_id_uniq,
    tumor, file_name_tumor,
    normal, file_name_normal,
    type, project, file_id_tumor, file_id_normal
  )

# Cases with only normal data
df %>%
  filter(case %in% setdiff(unique(df$case), unique(pair_info$case)))

length(unique(pair_info$case))
table(pair_info$type) # Check ids of metastatic samples

# NOTE in our analyses, we don't need handle them all,
# because only should be a sample per case has ecDNA record
# I did extra uncessary work in our train data set and exploration step
# so I had to remove incorrect data and redo. 

load("data/pancan_amplicon_list_and_summary.RData")

head(pair_info$tumor)
head(pair_info$normal)

pair_info2 = pair_info %>% 
  dplyr::mutate(tumor = substr(tumor, 1, 15),
                normal = substr(normal, 1, 15)) %>% 
  dplyr::filter(tumor %in% data_summary_tcga$sample_barcode)

write_csv(pair_info2, file = "data/outer_test_wes_pair_info.csv") # NOTE: this is not used as outer test as the original exploration now
                                                                  # All WES data have been merged together

# ---------------------------------------------------------------------------------
# Non-Ec

# Tumor normal pairing for tumor without regions

library(tidyverse)

df <- read_tsv("../gdc-manifest/gdc_sample_sheet.2021-12-23-2.tsv")

table(df$`Sample Type`)

df <- df %>%
  select(-c(3:4)) %>%
  setNames(c("file_id", "file_name", "project", "case", "sample", "type")) %>%
  mutate(type2 = ifelse(grepl("Normal", type), "normal", "tumor"))

table(df$type)
table(df$type2)

# sum(grepl("hg19", df$file_name))
# Although many file names have been laebelled as hg19
# but the GDC portal indicates their reference genome is hg38
# e.g., https://portal.gdc.cancer.gov/files/04a38123-6191-4848-8b0e-3a6b1749e2a7

build_df <- function(df) {
  message("handling sample ", paste(df$sample, collapse = ","))
  
  if (any(grepl("Solid", df$type)) & any(grepl("Blood", df$type))) {
    # Use blood normal as reference
    df <- df %>% filter(!grepl("Solid", df$type))
  }
  
  r <- pivot_wider(df,
                   id_cols = "case",
                   names_from = "type2", values_from = "sample", values_fn = list
  )
  if (!"tumor" %in% colnames(r)) {
    r <- dplyr::tibble(
      tumor = rep(NA_character_, length(unique(r$normal[[1]]))),
      normal = unique(r$normal[[1]])
    )
  } else if (!"normal" %in% colnames(r)) {
    r <- dplyr::tibble(
      tumor = unique(r$tumor[[1]]),
      normal = rep(NA_character_, length(unique(r$tumor[[1]])))
    )
  } else {
    r <- expand.grid(r$tumor[[1]], r$normal[[1]]) %>%
      setNames(c("tumor", "normal"))
  }
  
  nrows <- cumprod(table(df$type2))[2]
  nr <- nrow(r)
  idxT <- df$sample %in% r$tumor
  idxN <- df$sample %in% r$normal
  
  if (!is.na(nrows) && nrows > 1) {
    bind_cols(
      r[rep(seq_len(nr), nrows / nr), ],
      df[idxT, c("file_id", "file_name")][rep(seq_len(sum(idxT)), nrows / sum(idxT)), ] %>% # nolint
        setNames(c("file_id_tumor", "file_name_tumor")),
      df[idxN, c("file_id", "file_name")][rep(seq_len(sum(idxN)), nrows / sum(idxN)), ] %>% # nolint
        setNames(c("file_id_normal", "file_name_normal"))
    ) %>%
      as_tibble()
  } else {
    bind_cols(
      r,
      if (any(idxT)) {
        df[idxT, c("file_id", "file_name")] %>%
          setNames(c("file_id_tumor", "file_name_tumor"))
      } else {
        NULL
      },
      if (any(idxN)) {
        df[idxN, c("file_id", "file_name")] %>%
          setNames(c("file_id_normal", "file_name_normal"))
      }
    ) %>%
      as_tibble()
  }
}

pair_info <- df %>%
  group_split(case) %>%
  map_df(build_df) %>%
  filter(!is.na(tumor)) %>%
  left_join(df %>%
              filter(type2 == "tumor") %>%
              select(-file_id, -file_name, -type2) %>%
              unique(),
            by = c("tumor" = "sample")
  ) %>%
  mutate(pair_id_uniq = make.unique(tumor)) %>%
  select(
    case, pair_id_uniq,
    tumor, file_name_tumor,
    normal, file_name_normal,
    type, project, file_id_tumor, file_id_normal
  )

# Cases with only normal data
df %>%
  filter(case %in% setdiff(unique(df$case), unique(pair_info$case)))

length(unique(pair_info$case))
table(pair_info$type) # Check ids of metastatic samples

# NOTE in our analyses, we don't need handle them all,
# because only should be a sample per case has ecDNA record

load("data/pancan_amplicon_list_and_summary.RData")

head(pair_info$tumor)
head(pair_info$normal)

pair_info2 = pair_info %>% 
  dplyr::mutate(tumor = substr(tumor, 1, 15),
                normal = substr(normal, 1, 15)) %>% 
  dplyr::filter(tumor %in% data_summary_tcga$sample_barcode)

write_csv(pair_info2, file = "data/nonEC_wes_pair_info.csv")
