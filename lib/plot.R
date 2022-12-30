plot_clean_stat_box = function(data, x, y, ylab = "value", colors = c("grey", "blue", "red"),
                               centrality.plotting = TRUE, type = "np", centrality.type = type) {
  library(ggplot2)
  library(ggstatsplot)
  ggbetweenstats(
    data = data,
    x = {{x}},
    y = {{y}}, 
    type = type,
    centrality.type = centrality.type,
    point.args = list(alpha = 0),
    centrality.point.args = list(alpha = 0),
    centrality.plotting = centrality.plotting,
    centrality.label.args = list(size = 2, nudge_x = 0.4, segment.linetype = 4,
         min.segment.length = 0),
    point.path = FALSE,
    violin.args = list(width = 0.5, alpha = 0.2, mapping = aes(fill = {{x}})),
    results.subtitle = FALSE,
    bf.message = F
  ) -> p
  p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
    scale_fill_manual(values = colors) + labs(x = NULL, y = ylab, caption = NULL)
  p
}

# plot_clean_stat_box(data_tcga, x = class, y = cna_burden, ylab = "Copy number alteration burden")



