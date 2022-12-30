#!/usr/bin/env Rscript
# Train XGB32 with different folds and finally selected hyper-parameters

PROJ_DIR = "/data3/wsx_data/gcap-analysis/modeling/"
setwd(PROJ_DIR)

#args <- commandArgs(trailingOnly = TRUE)
idx <- 1 # as.integer(args[1])

out_path <- "data/xgb_v3/"
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

nthread <- 10

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
tm = readRDS("data/xgb_v3/XGB_NF32.rds")
params = tm$params
params$nthread = nthread

result = list()

for (i in folds) {
  data_list <- get_model_data(rm_type = rm_type, rm_cli = rm_cli, rm_cns = rm_cns, Nfolds = i)
  
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
          nrounds = 1000,
          showsd = TRUE,
          print_every_n = 1,
          early_stopping_rounds = 10,
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

saveRDS(result, file = "data/xgb_v3/XGB32_with_diff_folds.rds")

str(result$`3`[[1]], max.level = 1)
result$`3`[[1]]$evaluation_log

result_df = rbindlist(
  lapply(result, function(ls) {
    purrr::map_df(ls, ~.$evaluation_log[.$best_ntreelimit])
  }), idcol = "folds"
)

result_df
saveRDS(result_df, file = "data/xgb_v3/XGB32_with_diff_folds.rds")

library(ggpubr)
p1 = ggpubr::ggbarplot(result_df[, .(folds, aucpr = round(train_aucpr_mean, 3))],
                  x = "folds", y = "aucpr", 
                  add = "mean_sd",
                  fill = "steelblue", label = TRUE, lab.vjust = -0.5, lab.nb.digits = 3) +
  labs(y = "Train auPR")

p2 = ggpubr::ggbarplot(result_df[, .(folds, aucpr = round(test_aucpr_mean, 3))],
                  x = "folds", y = "aucpr", 
                  add = "mean_sd",
                  fill = "steelblue", label = TRUE, lab.vjust = -0.5, lab.nb.digits = 3) +
  labs(y = "Test auPR")

p3 = cowplot::plot_grid(p1, p2, align = "hv", ncol = 2)
cowplot::save_plot("plots/XGB32_auPR_with_different_folds.pdf", ncol = 2, plot = p3, base_height = 3.5, base_width = 2.4)
