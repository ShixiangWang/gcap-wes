wd = "/data3/wsx_data/gcap-analysis"
stopifnot(dir.exists(wd))

library(dplyr)
library(cli)

get_perf <- function(fit, idx) {
  old_wd = getwd()
  setwd(wd)
  on.exit(setwd(old_wd))
  ## Use the idx to determine the target
  if (idx == 1L) {
    rm_type = TRUE
    rm_cli = FALSE
    rm_cns = FALSE
  } else if (idx == 2L) {
    rm_type = TRUE
    rm_cli = TRUE
    rm_cns = TRUE
  } else if (idx == 3L) {
    rm_type = FALSE
    rm_cli = FALSE
    rm_cns = FALSE
  }
  
  library(stats)
  suppressMessages(library(xgboost))
  suppressMessages(library(caret))
  suppressMessages(library(purrr))
  stopifnot(requireNamespace("gcap", quietly = TRUE))
  
  ## Preprocess data
  source("lib/get_data.R")
  data_list = get_model_data(rm_type = rm_type, rm_cli = rm_cli, rm_cns = rm_cns)
  
  dtrain = xgb.DMatrix(data_list$data_train[, !colnames(data_list$data_train) %in% "y"],
                       label = data_list$data_train[, "y"])
  dtest = xgb.DMatrix(data_list$data_test[, !colnames(data_list$data_test) %in% "y"],
                      label = data_list$data_test[, "y"])
  
  Features = colnames(data_list$data_train)
  Features = Features[-length(Features)]
  cli_alert(c("Model: ", length(Features), " features"))
  
  rm(data_list); invisible(gc())
  
  cli_alert("Eval probs in datasets")
  train = predict(fit, dtrain, ntreelimit = fit$best_ntreelimit)
  test = predict(fit, dtest, ntreelimit = fit$best_ntreelimit)
  
  eval_overall_perf = function(x, y) {
    rv = gcap::get_auc(x, getinfo(y, "label"), curve = TRUE, type = "roc")
    rv2 = gcap::get_auc(x, getinfo(y, "label"), curve = TRUE)
    
    list(
      ROC = list(
        data = rv$curve,
        auc = round(rv$auc, 3),
        type = "auROC"
      ),
      PR = list(
        data = rv2$curve,
        auc = round(rv2$auc.integral, 3),
        type = "auPR"
      )
    )
  }
  
  cli_alert("Eval overvall performance")
  overall_perf = list(
    train = eval_overall_perf(train, dtrain),
    test = eval_overall_perf(test, dtest)
  )
  
  cli_alert("Eval performance using different prob as classifier threshold")
  df <- map_df(c(0.01, 0.05, seq(0.1, 0.95, 0.05)), function(cutoff) {
    cli_alert("> cutoff: {cutoff}")

    pred_train2 <- ifelse(train > cutoff, 1L, 0L)
    pred_test2 <- ifelse(test > cutoff, 1L, 0L)
    
    #print(summary(pred_train2))
    #print(summary(pred_test2))
    
    if (all(pred_train2 == 0) | all(pred_train2 == 1)) {
      cli_alert("No proper data, skip")
      return(data.frame())
    }
    
    # Sensitivity is also known as recall
    cfmat_train <- confusionMatrix(table(Prediction = pred_train2, Reference = getinfo(dtrain, "label")), positive = "1")
    out1 <- cfmat_train[["byClass"]][c("Sensitivity", "Specificity", "Precision", "F1")]
    
    cfmat_test <- confusionMatrix(table(Prediction = pred_test2, Reference = getinfo(dtest, "label")), positive = "1")
    out2 <- cfmat_test[["byClass"]][c("Sensitivity", "Specificity", "Precision", "F1")]
    
    data.frame(cutoff = cutoff,
               data_type = rep(c("train", "test"), each = 4),
               measure = rep(c("Sensitivity", "Specificity", "Precision", "F1"), 2),
               score = c(as.numeric(out1), as.numeric(out2)))
  })
  cli_alert_info("Performance measures under different prob cutoff:")
  print(df)
  
  cli_alert_success("Done")
  list(
    data = list(train = data.frame(y = getinfo(dtrain, "label"), y_hat = train),
                test = data.frame(y = getinfo(dtest, "label"), y_hat = test)),
    perf_overall = overall_perf,
    perf_stats = df
  )
}


plot_perf <- function(in_list, type = "auROC") {
  suppressMessages(library(ggplot2))
  
  if (type == "auROC") aucs = c(in_list$train$ROC$auc, in_list$test$ROC$auc)
  else aucs = c(in_list$train$PR$auc, in_list$test$PR$auc)
  
  types = c(paste0("Train (AUC: ", aucs[1], ")"), 
            paste0("Test  (AUC: ", aucs[2], ")"))
  
  idx = if (type == "auROC") 1 else 2
  data = rbind(cbind(data.frame(in_list[[1]][[idx]]$data), type = types[1]),
               cbind(data.frame(in_list[[2]][[idx]]$data), type = types[2]))
  
  ggplot(data, aes(x = X1, y = X2, color = type)) +
    geom_line() +
    cowplot::theme_cowplot() +
    labs(x = ifelse(type == "auPR", "Sensitivity (Recall)",
                    "False positive rate"
    ), y = ifelse(type == "auPR", "Precision", "True positive rate"), color = NULL) +
    scale_color_manual(values = c(
      "red",
      "blue"
    ) %>%
      setNames(types)) +
    theme(legend.position = if (type == "PR") {
      c(0.05, 0.2)
    } else {
      c(0.3, 0.2)
    }) +
    expand_limits(x = 0, y = 0)
}
