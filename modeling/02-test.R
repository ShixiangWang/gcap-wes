library(xgboost)

source("03-modeling-renew/get_data_renew.R")
data_list = get_model_data(rm_type = TRUE, rm_cli = FALSE, rm_cns = FALSE)

dtrain = xgb.DMatrix(data_list$data_train[, !colnames(data_list$data_train) %in% "y"],
                     label = data_list$data_train[, "y"])
dtest = xgb.DMatrix(data_list$data_test[, !colnames(data_list$data_test) %in% "y"],
                    label = data_list$data_test[, "y"])

params <- list(
  max_delta_step = 1,
  eval_metric = "aucpr",
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.3,
  gamma = 0,
  max_depth = 3,
  min_child_weight = 5,
  colsample_bytree = 1,
  nthread = 20
)

set.seed(2021)
bst_cv = xgb.cv(
  params = params, 
  data = dtrain, 
  nrounds = 100,
  nfold = 5,
  verbose = 2, 
  print_every_n = 1,
  early_stopping_rounds = 5
)
#>
# 模型本身没有过拟合，
# 但实际上，上述方式造成了样本信息的相互泄露
# 数据的生物学意义和数值表征存在pitfall
#<

# It's a pitfall
saveRDS(bst_cv$evaluation_log, file = "data/renew-test-log.rds")

params <- list(
  max_delta_step = 6,
  eval_metric = "aucpr",
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.177,
  gamma = 0.1,
  max_depth = 4,
  min_child_weight = 5,
  colsample_bytree = 0.5,
  lambda = 0.1,
  nthread = 20
)
set.seed(2021)
bst_cv = xgb.cv(
  params = params, 
  data = dtrain, 
  nrounds = 50,
  folds = data_list$cv_list[[1]],
  verbose = 2, 
  print_every_n = 1,
  early_stopping_rounds = 5
)

params <- list(
  max_delta_step = 1,
  eval_metric = "aucpr",
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.3,
  gamma = 0,
  max_depth = 7,
  min_child_weight = 1000,
  colsample_bytree = 1,
  nthread = 20
)

set.seed(2021)
bst = xgb.train(
  params = params,
  data = dtrain,
  watchlist = list(train = dtrain, test = dtest),
  nrounds = 100,
  print_every_n = 1
)

gcap::get_auc(predict(bst, dtest), getinfo(dtest, "label"))$auc.integral


# Test repeats ------------------------------------------------------------

data = readRDS("data/train_data_v3_raw.rds")
data_repeats = readRDS("/data3/wsx_data/rmsk_gene_test.rds")
rmsk_gene_sum = data_repeats[, .N, by = .(gene_id, V4)]

rmsk_gene_sum$gene_id = factor(rmsk_gene_sum$gene_id, levels = unique(substr(data$gene_id, 1, 15)))
rmsk_gene_sum = data.table::dcast(rmsk_gene_sum, gene_id ~ V4, fill = 0L, value.var = "N", drop = FALSE)

data2 = merge(data[, .(sample, gene_id, y)], unique(rmsk_gene_sum), by = "gene_id")
#data2 = data[, .(sample, CN4, CN5, CN6, CN7, CN8, y)]

set.seed(2021)
id_list = unique(data2$sample)
id_list2 = list(train = sample(id_list, 200))
id_list2$test = setdiff(id_list, id_list2$train)

library(xgboost)

dtrain = xgb.DMatrix(as.matrix(data2[sample %in% id_list2$train, -c("gene_id", "sample", "y")]),
                     label = data2[sample %in% id_list2$train]$y)
dtest = xgb.DMatrix(as.matrix(data2[sample %in% id_list2$test, -c("gene_id", "sample", "y")]),
                    label = data2[sample %in% id_list2$test]$y)

params <- list(
  max_delta_step = 1,
  eval_metric = "aucpr",
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.3,
  gamma = 0,
  max_depth = 7,
  min_child_weight = 1000,
  colsample_bytree = 1,
  nthread = 20
)

bst = xgb.train(
  params = params,
  data = dtrain,
  watchlist = list(train = dtrain, test = dtest),
  nrounds = 100,
  print_every_n = 1
)

gcap::get_auc(predict(bst, dtest), getinfo(dtest, "label"))$auc.integral

# No association