PROJ_DIR = "~/gcap-analysis/modeling/"
setwd(PROJ_DIR)

library(xgboost)
library(dplyr)

model = c("NF11", "NF32", "NF56")
# XGB models ---------------------------------

for (m in model) {
  fl = list.files("data/xgb_v4_tune/", pattern = paste0("xgb_cv_result_", m), full.names = T)
  fl = fl[file.info(fl)$size > 200]
  data_list = purrr::map(fl, readRDS)
  perf_list = purrr::map(data_list, function(l) {
    purrr::map_df(l, ~.$evaluation_log[.$best_ntreelimit, ])
  })
  
  perf_df = purrr::map_df(perf_list, function(l) {
    data.frame(iter_mean = mean(l$iter),
               iter_max = max(l$iter),
               train_aucpr_mean = mean(l$train_aucpr_mean),
               train_aucpr_sd = mean(l$train_aucpr_std),
               test_aucpr_mean = mean(l$test_aucpr_mean),
               test_aucpr_sd = mean(l$test_aucpr_std),
               aucpr_diff = mean(l$test_aucpr_mean - l$train_aucpr_mean))
  }, .id = "cv") %>%
    mutate(file = fl, params = purrr::map(data_list, ~.[[1]]$params)) %>% 
    #filter(test_aucpr_mean > 0.7) %>% 
    #arrange(abs(aucpr_diff)) %>% 
    arrange(desc(round(test_aucpr_mean, 3)), round(test_aucpr_sd, 3), desc(iter_mean)) %>% 
    head(20)
  
  # Choose top 1 for minimal aucpr_diff
  # perf_df = purrr::map_df(perf_list, function(l) {
  #   data.frame(iter_mean = mean(l$iter),
  #              iter_max = max(l$iter),
  #              train_aucpr_mean = mean(l$train_aucpr_mean),
  #              train_aucpr_sd = mean(l$train_aucpr_std),
  #              test_aucpr_mean = mean(l$test_aucpr_mean),
  #              test_aucpr_sd = mean(l$test_aucpr_std),
  #              aucpr_diff = mean(l$test_aucpr_mean - l$train_aucpr_mean))
  # }, .id = "cv") %>%
  #   mutate(file = fl, params = purrr::map(data_list, ~.[[1]]$params)) %>% 
  #   filter(test_aucpr_mean > 0.7) %>% 
  #   arrange(abs(aucpr_diff)) %>% 
  #   head(20) %>% 
  #   arrange(desc(round(test_aucpr_mean, 3)), round(test_aucpr_sd, 3), desc(iter_mean))
  
  print(perf_df)
  
  # Pick the first set
  pl = perf_df$params[[1]]
  print(pl)
  
  # Set parameters as pl
  params <- pl
  
  source("../lib/get_data.R")
  data <- get_model_data2(rm_type = ifelse(m == "NF56", FALSE, TRUE), 
                          rm_cli = ifelse(m == "NF11", TRUE, FALSE), 
                          rm_cns = ifelse(m == "NF11", TRUE, FALSE), 
                          version = "v4", fully_random = TRUE)
  dtrain = xgb.DMatrix(data$data[, !colnames(data$data) %in% "y"],
                       label = data$data[, "y"])
  
  head(perf_df, 3)
  
  set.seed(2021)
  bst = xgb.train(
    params = params,
    data = dtrain,
    watchlist = list(train = dtrain),
    nrounds = perf_df$iter_max[1],
    print_every_n = 1
  )
  
  # aurpc
  gcap::get_auc(predict(bst, dtrain), xgboost::getinfo(dtrain, "label"), curve = TRUE) %>% plot()
  
  # auc
  gcap::get_auc(predict(bst, dtrain), xgboost::getinfo(dtrain, "label"), type = "roc", curve = TRUE) %>% plot()
  
  dir.create("data/xgb_v4")
  saveRDS(bst, file = sprintf("data/xgb_v4/XGB_%s.rds", m))
  saveRDS(perf_df, file = sprintf("data/xgb_v4/XGB_%s_perf_top20.rds", m))
}

#saveRDS(bst, file = sprintf("data/xgb_v4/XGB_%s.rds", "NF32_mindiff"))
#saveRDS(perf_df, file = sprintf("data/xgb_v4/XGB_%s_perf_top20.rds", "NF32_mindiff"))
