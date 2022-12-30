enrichBP <- function(fCNA, freq = c(5, 0), usebreak = FALSE) {
    library(clusterProfiler)
    library(gcap)
    library(org.Hs.eg.db)
    library(dplyr)

    geneList_linear <- fCNA$gene_summary[noncircular > freq[1]]$gene_id
    message("Linear genes: ", length(geneList_linear))
    geneList_circular <- fCNA$gene_summary[(circular + possibly_circular) > freq[2]]$gene_id
    message("Circular genes: ", length(geneList_circular))

    if(usebreak) return()

    geneList <- list(
        linear = unique(bitr(geneList_linear, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID),
        circular = unique(bitr(geneList_circular, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
    )

    ck_go <- compareCluster(geneCluster = geneList, fun = enrichGO, ont = "BP", OrgDb = org.Hs.eg.db, readable = TRUE)
    ck_go <- simplify(ck_go)

    circular_and_linear_gene_cluster_GO_BP_comparison <- ck_go@compareClusterResult
    print(head(circular_and_linear_gene_cluster_GO_BP_comparison))
    z <- circular_and_linear_gene_cluster_GO_BP_comparison %>%
        group_by(Cluster) %>%
        summarise(ls = stringr::str_split(geneID, pattern = "/") %>%
            unlist()) %>%
        group_split() %>%
        purrr::map(~ unique(.$ls))
    list(
        result = ck_go,
        tidy_result = z
    )
}
