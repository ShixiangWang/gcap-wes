#!/usr/bin/env Rscript
# Train XGB with different folds and finally selected hyper-parameters

args <- commandArgs(trailingOnly = TRUE)
idx = as.integer(args[1])
#idx <- 1 # as.integer(args[1])

out_path <- "data"
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

nthread <- 38

library(xgboost)

## Use the idx to determine the target

if (idx == 1L) {
  rm_type <- TRUE
  rm_cli <- FALSE
  rm_cns <- FALSE
} else if (idx == 2L) {
  rm_type <- TRUE
  rm_cli <- TRUE
  rm_cns <- TRUE
} else if (idx == 3L) {
  rm_type <- FALSE
  rm_cli <- FALSE
  rm_cns <- FALSE
}

## Preprocess data
source("../lib/get_data.R")
folds = c(3, 5, 10, 20)
tm = if (idx == 1) {
  readRDS("../modeling/data/xgb_v4/XGB_NF32.rds")
} else if (idx == 2) {
  readRDS("../modeling/data/xgb_v4/XGB_NF11.rds")
} else {
  readRDS("../modeling/data/xgb_v4/XGB_NF56.rds")
}
params = tm$params
params$nthread = nthread

result = list()

for (i in folds) {
  data_list <- get_model_data2(rm_type = rm_type,
                               rm_cli = rm_cli, 
                               rm_cns = rm_cns, 
                               Nfolds = i,
                               version = "v4", fully_random = TRUE)
  names(data_list)[1] = "data_train"
  Features <- colnames(data_list$data_train)
  Features <- Features[-length(Features)]
  message("Features used:")
  print(Features)
  name_addon <- paste0("NF", length(Features))
  
  set.seed(2021)
  ## XGBOOST Cross Validation
  library(parallel)
  
  data = data_list$data_train[, colnames(data_list$data_train) != "y"]
  label = data_list$data_train[, "y"]
  cv_list = data_list$cv_list
  rm(data_list); gc()
  
  rv = list()
  for (j in 1:3) {
    message("Iteration #", j)
    cv_fit = tryCatch(
      {
        xgb.cv(
          params = params,
          data = data,
          label = label,
          nrounds = tm$niter,
          showsd = TRUE,
          print_every_n = 1,
          #early_stopping_rounds = 10,
          metrics = "aucpr",
          verbose = TRUE,
          folds = cv_list[[j]]  # See doc for detail
        )
      },
      error = function(e) {
        print(e$message)
        message("Model CV is failed for Nfold: ", i)
        NULL
      }
    )
    print(cv_fit$evaluation_log)
    rv[[j]] = cv_fit[c("evaluation_log", "params", "best_ntreelimit")]
    rm(cv_fit); gc()
  }
  
  result[[as.character(i)]] = rv
  message("Done for folds: ", i)
}

saveRDS(result, file = sprintf("data/%s_with_diff_folds.rds", name_addon))

