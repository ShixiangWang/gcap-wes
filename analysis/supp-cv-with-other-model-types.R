# https://github.com/topepo/caret/issues/1193#issuecomment-776220846

folds = c(3, 5, 10, 20)

all_result = list()
source("../lib/get_data.R")

library(foreach)
library(doParallel)

doParallel::registerDoParallel(cores = 4)

all_result = foreach (f = folds, .verbose = TRUE) %dopar% {
  cat("handling fold ", f, "\n", file = "progress", append = TRUE)
  data_list <- get_model_data2(rm_type = TRUE, rm_cli = TRUE, rm_cns = TRUE, 
                               version = "v4", fully_random = TRUE, Nfolds = f)
  data = as.data.frame(data_list$data)
  colnames(data)
  cv_list = data_list$cv_list
  
  auc_rv = list()
  for (j in 1:3) {
    auc_list = c()
    for (i in cv_list[[j]]) {
      mod_fit <- glm(y ~ ., data = data[-i, ], family = binomial())
      
      y_pred <- predict(mod_fit, type = "response")
      y_true = data[-i, "y"]
      cat("train: \n", file = "progress", append = TRUE)
      cat(gcap::get_auc(y_pred, y_true, type = "pr")$auc.integral, "\n", file = "progress", append = TRUE)
      
      cat("test: \n", file = "progress", append = TRUE)
      y_pred <- predict(mod_fit, data[i,], type = "response")
      y_true = data[i, "y"]
      mod_aucpr = gcap::get_auc(y_pred, y_true, type = "pr")$auc.integral
      cat(mod_aucpr, "\n", file = "progress", append = TRUE)
      auc_list = c(auc_list, mod_aucpr)
    }
    auc_rv[[j]] = auc_list
  }
  #all_result[[as.character(f)]] = auc_rv
}

all_result = list(
  "3" = list(
    c(0.7293068, 0.6916844, 0.7169736),
    c(0.6600780, 0.7612615, 0.6912391),
    c(0.7088060, 0.6640384, 0.7354969)
  ),
  "5" = list(
    c(0.7642400, 0.6992629, 0.7223684, 0.7500590, 0.6300367),
    c(0.7336811, 0.7784011, 0.6820034, 0.7167809, 0.6537468),
    c(0.6947834, 0.6575073, 0.7759784, 0.7349271, 0.7119696)
  ),
  "10" = list(
    c(0.8327682, 0.6443673, 0.7027079, 0.7571767, 0.6467179, 
      0.7801534, 0.4829009, 0.6622523, 0.8554683, 0.7549624),
    c(0.6917729, 0.8741452, 0.8174191, 0.6120799, 0.7363005,
      0.8637159, 0.7140282, 0.6876676, 0.7340060, 0.6175352),
    c(0.7252827, 0.7895338, 0.8536189, 0.6539749, 0.6929545,
      0.7616684, 0.7550123, 0.6994819, 0.8070102, 0.5672349)
  ),
  "20" = list(
    c(0.6917046, 0.5362192, 0.9136018, 0.4401381, 0.7618686, 
      0.6698952, 0.7550165, 0.8124759, 0.7375635, 0.8284697,
      0.7856090, 0.7218766, 0.7880894, 0.7456920, 0.7911274,
      0.8067173, 0.6917731, 0.7873469, 0.6888845, 0.6739047),
    c(0.7726155, 0.8107380, 0.7046723, 0.8111980, 0.6350234,
      0.5331211, 0.7180136, 0.5574510, 0.6289307, 0.8413848,
      0.9261895, 0.8060359, 0.6282646, 0.6962094, 0.7135926,
      0.8239483, 0.6150224, 0.6770875, 0.5189798, 0.9008599),
    c(0.6261730, 0.6318065, 0.8954325, 0.7314405, 0.5971092,
      0.8844786, 0.4912823, 0.7513299, 0.6861764, 0.8216095)
  )
)

saveRDS(all_result, file = "data/model_perf_logistic.rds")


# library(caret)
# folds <- createFolds(mtcars$mpg, k = 5, returnTrain = TRUE)
# train_res <- train(mpg ~ ., data = mtcars, method = "lm", 
#                    trControl = trainControl(method = "cv", index = folds,
#                                             savePredictions = TRUE))
# train_res
# randomForest(iris[,-5],y=iris[,5],ntree=ntree,do.trace=100) 
