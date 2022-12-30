#!/usr/bin/env Rscript
# 与11-*.R编写的代码不同，这里更改了源数据，建模保证了数据不变
# 保留纯标签 circle, HR, linear数据，并对nofocal进行采样equal-size采样
PROJ_DIR = "/data3/wsx_data/gcap-analysis/modeling/"
setwd(PROJ_DIR)

args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])

out_path <- "data/xgb_tune_downsample"
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

nthread <- 20
N <- 1000 # Number of parameter pairs

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
data_list <- get_model_data2(rm_type = rm_type, rm_cli = rm_cli, rm_cns = rm_cns, 
                             version = "subsample", fully_random = TRUE)

Features <- colnames(data_list$data)
Features <- Features[-length(Features)]
message("Features used:")
print(Features)
name_addon <- paste0("NF", length(Features))

set.seed(2021)
## XGBOOST Cross Validation
library(parallel)
# random search for hyper-parameters controlling overfitting

message("generating search space...")

# https://xgboost.readthedocs.io/en/stable/tutorials/param_tuning.html#control-overfitting
# https://xgboost.readthedocs.io/en/stable/tutorials/param_tuning.html#handle-imbalanced-dataset
search_space <- data.frame(
  eta = sample(c(0.01, 0.1, 0.3), N, replace = TRUE, prob = c(0.2, 0.6, 0.2)),
  max_delta_step = sample(c(0, 1, 10), N, replace = TRUE, prob = c(0.2, 0.6, 0.2)),
  max_depth = sample(2:6, N, replace = TRUE),
  min_child_weight = sample(c(1, 2, 5, 10, 20, 100), N, replace = TRUE, prob = c(0.1, 0.3, 0.2, 0.2, 0.1, 0.1)),
  alpha = sample(c(0, 0.5, 1), N, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
  lambda = sample(c(1, 0.5, 0), N, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
  gamma = sample(c(0, 1, 10), N, replace = TRUE, prob = c(0.2, 0.4, 0.4)),
  subsample = sample(c(0.5, 0.6, 0.7, 0.8, 0.9), N, replace = TRUE, prob = rep(0.2, 5)),
  colsample_bytree = sample(c(0.6, 0.7, 0.8, 0.9, 1), N, replace = TRUE, prob = rep(0.2, 5))
)

message(nrow(search_space), " hyper-parameter pairs in total")

data = data_list$data[, colnames(data_list$data) != "y"]
label = data_list$data[, "y"]
cv_list = data_list$cv_list
rm(data_list); gc()

message("Repeating cross validation 3 times")

ss = 1:nrow(search_space)

for (i in ss) {
  message("handling parameter pair ", i, "...")
  params <- list(
    max_delta_step = search_space$max_delta_step[[i]],
    booster = "gbtree",
    objective = "binary:logistic",
    nthread = nthread,
    eta = search_space$eta[[i]],
    alpha = search_space$alpha[[i]],
    lambda = search_space$lambda[[i]],
    gamma = search_space$gamma[[i]],
    max_depth = search_space$max_depth[[i]],
    min_child_weight = search_space$min_child_weight[[i]],
    subsample = search_space$subsample[[i]],
    colsample_bytree = search_space$colsample_bytree[[i]]
  )
  print(params)
  
  fp <- file.path(
    out_path,
    paste0("xgb_cv_result_", name_addon, "_", i, ".rds")
  )
  message("model evaluation will be saved to ", fp)
  
  if (file.exists(fp)) {
    message("The hyperparameters have been done before, just skip it.")
  } else {
    rv = list()
    for (j in 1:1) {
      message("Iteration #", j)
      cv_fit = tryCatch(
        {
          xgb.cv(
            params = params,
            data = data,
            label = label,
            nrounds = 1000,
            showsd = TRUE,
            print_every_n = 5,
            early_stopping_rounds = 10,
            metrics = "aucpr",
            verbose = TRUE,
            folds = cv_list[[2]]
          )
        },
        error = function(e) {
          print(e$message)
          message("Model CV is failed for parameter pair id: ", i)
          NULL
        }
      )
      print(cv_fit$evaluation_log)
      rv[[j]] = cv_fit[c("evaluation_log", "params", "best_ntreelimit")]
      rm(cv_fit); gc()
    }
    saveRDS(rv, file = fp)
  }
}

message("Done")
