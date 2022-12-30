seq2 = function(x, ...) {
  lapply(x, function(i) {
    seq(i, ...)
  })
}

seq_random = function(x, n) {
  rl = as.integer(unlist(lapply(sigminer:::chunk2(seq_len(x), ceiling(x/n)), sample)))
  lapply(seq2(seq_len(n), x, n), function(x) rl[x])
}

convert_idx = function(ix, df, ids) {
  which(ids %in% df$sample[ix])
}

get_model_data = function(rm_type = FALSE, rm_cli = FALSE, rm_cns = FALSE,
                          Nfolds = 5, 
                          cv_repeats = Nfolds,
                          fully_random = FALSE,
                          only_samples = FALSE,
                          version = "v3") {
  cwd = getwd()
  towd = "/data3/wsx_data/gcap-analysis/modeling/"
  if (dir.exists(towd)) {
    setwd(towd)
  }
  on.exit(setwd(cwd))
  
  load(sprintf("data/train_data_extra_info_%s.RData", version)) # data_sample_record, data_gene_load
  data_gene_load = data.table::as.data.table(data_gene_load)
  
  # Split train (80%) and test (20%)
  idx_list = seq2(seq_len(5), nrow(data_gene_load), 5)
  
  # data split
  idx_tst = idx_list[[2]]  # Use the second split as test
  idx_trn = sort(setdiff(seq_len(nrow(data_gene_load)), idx_tst))
  
  if (only_samples) return(list(
    train = data_gene_load[idx_trn]$sample,
    test = data_gene_load[idx_tst]$sample
  ))
  
  mat = readRDS(sprintf("data/train_data_%s.rds", version))
  
  if (rm_type) {
    mat = mat[, !startsWith(colnames(mat), "type_")]
  }
  if (rm_cli) {
    mat = mat[, !colnames(mat) %in% c("age", "gender")]
  }
  if (rm_cns) {
    mat = mat[, !startsWith(colnames(mat), "CN")]
  }
  
  ids = inverse.rle(data_sample_record)
  
  # Now we convert sample index to data row index
  data_train = mat[convert_idx(idx_trn, data_gene_load, ids), ]
  data_test  = mat[convert_idx(idx_tst, data_gene_load, ids), ]
  
  # In train set, generate N folds for cross validation with N times
  data_load_cv = data_gene_load[idx_trn, ]
  ids_cv = ids[ids %in% data_load_cv$sample]
  
  cv_list = list()
  for (i in seq_len(cv_repeats)) {
    # 这里默认生成的list有点需要注意
    # 如果为3，则生成3次3-fold交叉验证的test fold index vector
    # 如果为10，则生成10次10-fold交叉验证的test fold index vector
    if (!fully_random) {
      # Randomly sampling in specific chunks
      cv_samps = seq_random(nrow(data_load_cv), Nfolds)
    } else {
      # Fully random
      cv_samps = sigminer:::chunk2(sample(seq_len(nrow(data_load_cv)), nrow(data_load_cv)), Nfolds)
      names(cv_samps) = NULL
    }
    
    cv_list[[i]] = lapply(cv_samps, function(fold) {
      convert_idx(fold, data_load_cv, ids_cv)
    })
  }
  
  list(
    data_train = data_train,
    data_test = data_test,
    cv_list = cv_list
  )
}

# This version uses all data to generate CV folds
get_model_data2 = function(rm_type = FALSE, rm_cli = FALSE, rm_cns = FALSE,
                           Nfolds = 10, 
                           cv_repeats = Nfolds,
                           fully_random = FALSE,
                           only_samples = FALSE,
                           version = "v4", seed = 2021) {
  cwd = getwd()
  towd = "/data3/wsx_data/gcap-analysis/modeling/"
  if (dir.exists(towd)) {
    setwd(towd)
  } else if (dir.exists("~/gcap-analysis/modeling/")) {
    setwd("~/gcap-analysis/modeling/")
  }
  on.exit(setwd(cwd))
  
  load(sprintf("data/train_data_extra_info_%s.RData", version)) # data_sample_record, data_gene_load
  data_gene_load = data.table::as.data.table(data_gene_load)
  
  if (only_samples) return(data_gene_load$sample)
  
  mat = readRDS(sprintf("data/train_data_%s.rds", version))
  
  if (rm_type) {
    mat = mat[, !startsWith(colnames(mat), "type_")]
  }
  if (rm_cli) {
    mat = mat[, !colnames(mat) %in% c("age", "gender")]
  }
  if (rm_cns) {
    mat = mat[, !startsWith(colnames(mat), "CN")]
  }
  
  ids = inverse.rle(data_sample_record)
  
  # Generate N folds for cross validation with N times
  data_load_cv = data_gene_load
  ids_cv = ids
  
  cv_list = list()
  for (i in seq_len(cv_repeats)) {
    # 这里默认生成的list有点需要注意
    # 如果为3，则生成3次3-fold交叉验证的test fold index vector
    # 如果为10，则生成10次10-fold交叉验证的test fold index vector
    if (!fully_random) {
      # Randomly sampling in specific chunks
      set.seed(seed)
      cv_samps = seq_random(nrow(data_load_cv), Nfolds)
    } else {
      # Fully random
      cv_samps = sigminer:::chunk2(sample(seq_len(nrow(data_load_cv)), nrow(data_load_cv)), Nfolds)
      names(cv_samps) = NULL
    }
    
    cv_list[[i]] = lapply(cv_samps, function(fold) {
      convert_idx(fold, data_load_cv, ids_cv)
    })
  }
  
  list(
    data = mat,
    cv_list = cv_list
  )
}

