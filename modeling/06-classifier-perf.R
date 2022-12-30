# Performance analysis for classifying sample into circular/others
# Focus on XGB32
setwd(PROJ_DIR)

library(gcap)
library(data.table)
source("lib/get_data.R")
data_split = get_model_data(only_samples = TRUE)
data_split = data.table(
  sample = c(data_split$train, data_split$test),
  set = c(rep("train", length(data_split$train)), rep("test", length(data_split$test)))
)

# True labels
data_label = readxl::read_excel("preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3)
data_label = data_label %>% dplyr::filter(sample_barcode %in% data_split$sample) %>% 
  dplyr::select(sample = sample_barcode, y_label = sample_classification) %>% 
  data.table::as.data.table()
data_label = merge(data_label, data_split, by = "sample")
data_label[, y_true_circular := fcase(y_label %in% c("Circular"), 1L, default = 0L)]
data_label[, y_true_noncircular := fcase(!y_label %in% c("No-fSCNA", "Circular"), 1L, default = 0L)]
data_label[, y_true_amplicon := fcase(!y_label %in% c("No-fSCNA"), 1L, default = 0L)]
data_label[, y_label := fcase(y_label %in% c("Circular"), "circular",
                              y_label %in% c("No-fSCNA"), "nofocal",
                              default = "noncircular")]

table(data_label$set, data_label$y_label)
#           circular nofocal noncircular
# test        65       5           7
# train      261      25          23

data = readRDS("modeling/data/train_data_v3_raw.rds")
colnames(data)[3] <- "total_cn"
data$prob = gcap.runPrediction(data, model = "XGB32")

get_perf = function(df) {
  library(ModelMetrics)
  ms = c("auc", "accuracy", "precision", "recall", "specificity", "f1Score")
  #ms = c("accuracy", "auROC", "auPR")
  typs = c("circular", "noncircular", "amplicon")
  allcb = expand.grid(ms = ms, typs = typs, stringsAsFactors = FALSE)
  
  rs = list()
  for (i in 1:nrow(allcb)) {
    #message(allcb$ms[i], "-", allcb$typs[i])
    if (allcb$ms[i] == "accuracy") {
      rs[[paste(allcb$ms[i], allcb$typs[i], sep = "_")]] = mean(df[[paste0("y_true_", allcb$typs[i])]] == df[[paste0("y_pred_", allcb$typs[i])]], na.rm = TRUE)
    } else {
      rs[[paste(allcb$ms[i], allcb$typs[i], sep = "_")]] = do.call(
        allcb$ms[i],
        args = list(df[[paste0("y_true_", allcb$typs[i])]], df[[paste0("y_pred_", allcb$typs[i])]]))
      # rs[[paste(allcb$ms[i], allcb$typs[i], sep = "_")]] = if (allcb$ms[i] == "auPR") {
      #   gcap::get_auc(df[[paste0("y_pred_", allcb$typs[i])]], df[[paste0("y_true_", allcb$typs[i])]], type = "pr")$auc.integral
      # } else {
      #   gcap::get_auc(df[[paste0("y_pred_", allcb$typs[i])]], df[[paste0("y_true_", allcb$typs[i])]], type = "roc")$auc
      # }
      
    }
  }
  rs[["accuracy"]] = mean(df$y_pred == df$y_label, na.rm = TRUE)
  as.data.frame(rs)
}


get_all_perf = function(i) {
  message("Running result for #", i)
  
  data.table::setDTthreads(1)
  data_pred = gcap.runScoring(data,
                              tightness = search_space$tightness[i], 
                              gap_cn = search_space$gap_cn[i])$fCNA$sample_summary[, .(sample, y_pred = class)]
  data_pred[, y_pred_circular := fcase(y_pred %in% c("circular", "possibly_circular"), 1L, default = 0L)]
  data_pred[, y_pred_noncircular := fcase(y_pred %in% c("noncircular"), 1L, default = 0L)]
  data_pred[, y_pred_amplicon := fcase(!y_pred %in% c("nofocal"), 1L, default = 0L)]
  data_pred[, y_pred := fcase(y_pred %in% c("circular", "possibly_circular"), "circular",
                              y_pred == "noncircular", "noncircular",
                              default = "nofocal")]
  
  data_y = dplyr::full_join(data_label, data_pred, by = "sample")
  print(table(data_y$y_label, data_y$y_pred))
  
  rv = rbind(
    cbind(
      data.frame(
        set = "train",
        tightness = search_space$tightness[i],
        gap_cn = search_space$gap_cn[i]
      ),
      get_perf(data_y[set == "train"])),
    cbind(
      data.frame(
        set = "test",
        tightness = search_space$tightness[i],
        gap_cn = search_space$gap_cn[i]
      ),
      get_perf(data_y[set == "test"]))
  )
  print(rv)
  rv
}

search_space = expand.grid(
  tightness = 0:5,
  gap_cn = 2:6, 
  stringsAsFactors = FALSE)
result_grid = lapply(1:nrow(search_space), get_all_perf)

library(future)
plan(multisession, gc = TRUE, workers = 6)
options(future.globals.maxSize = 5e9)

result_grid = future.apply::future_lapply(1:nrow(search_space), get_all_perf)
result_grid = rbindlist(result_grid)
head(result_grid)

result_grid$set = factor(result_grid$set, c("train", "test"))

saveRDS(result_grid, file = "modeling/data/classifier_performance.rds")
 
# Overall accuracy
p = ggplot(result_grid, aes(tightness, gap_cn, fill = accuracy)) +
  geom_tile() +
  geom_text(aes(label = round(accuracy, 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_overall_accuracy.pdf", plot = p, width = 8, height = 4)

p = ggplot(result_grid, aes(tightness, gap_cn, fill = auc_circular)) +
  geom_tile() +
  geom_text(aes(label = round(auc_circular, 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_circular_auc.pdf", plot = p, width = 8, height = 4)

p = ggplot(result_grid, aes(tightness, gap_cn, fill = accuracy_circular )) +
  geom_tile() +
  geom_text(aes(label = round(accuracy_circular , 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_circular_accuracy.pdf", plot = p, width = 8, height = 4)

p = ggplot(result_grid, aes(tightness, gap_cn, fill = precision_circular  )) +
  geom_tile() +
  geom_text(aes(label = round(precision_circular  , 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_circular_precision.pdf", plot = p, width = 8, height = 4)

p = ggplot(result_grid, aes(tightness, gap_cn, fill = recall_circular)) +
  geom_tile() +
  geom_text(aes(label = round(recall_circular, 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_circular_recall.pdf", plot = p, width = 8, height = 4)


p = ggplot(result_grid, aes(tightness, gap_cn, fill = specificity_circular)) +
  geom_tile() +
  geom_text(aes(label = round(specificity_circular, 3))) +
  scale_fill_viridis_c() +
  facet_wrap(~set) +
  theme_bw()
p
ggsave("modeling/plots/sample_class_circular_specificity.pdf", plot = p, width = 8, height = 4)


library(dplyr)
p = result_grid %>% 
  dplyr::filter(tightness == 1, gap_cn == 4) %>% 
  dplyr::select(set, dplyr::ends_with("_circular")) %>% 
  tidyr::pivot_longer(cols = -"set", names_to = "type", values_to = "score") %>% 
  dplyr::mutate(type = stringr::str_remove(type, "_circular")) %>% 
  dplyr::filter(type %in% c("auc", "f1Score", "accuracy", "specificity")) %>% 
  dplyr::mutate(type = factor(type, c("auc", "f1Score", "accuracy", "specificity"))) %>% 
  ggplot(aes(x = type, y = score, fill = set)) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.6) +
  geom_text(aes(label = round(score, 2), group = set), position=position_dodge(1), vjust=-0.5) +
  cowplot::theme_cowplot() +
  labs(x = NULL)
p  
ggsave("modeling/plots/sample_class_circular_performance2.pdf", plot = p, width = 5, height = 4)


