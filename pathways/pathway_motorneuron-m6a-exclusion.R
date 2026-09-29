polo_gse <- readRDS("~/ALS/genes/GSEA/clean/polo_GSEA_cleaned.RDS")
swift_gse <- readRDS("~/ALS/genes/GSEA/clean/swift_GSEA_cleaned.RDS")
#M6a_gse <- readRDS("~/ALS/genes/GSEA/clean/M6a_GSEA_cleaned.RDS")
fus_gse <-readRDS("~/ALS/genes/GSEA/clean/fus_GSEA_cleaned.RDS")
multi_gse <- readRDS("~/ALS/genes/GSEA/clean/multi_GSEA_cleaned.RDS")
#install.packages("dplyr")
library(dplyr)
polo_sig <- polo_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D1")
swift_sig <- swift_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D2")
fus_sig <- fus_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D3")
multi_sig <- multi_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D4")
all_go <- bind_rows(multi_sig,fus_sig,swift_sig,polo_sig)
go_ids <- unique(all_go$ID)
library(GOSemSim)
library(org.Hs.eg.db)
hsGO <- godata('org.Hs.eg.db',ont="BP")

sim_matrix <- mgoSim(go_ids,
                     go_ids,
                     semData = hsGO,
                     measure = "Wang",
                     combine = NULL)



# Option 1: Use a real distance + valid hclust method
dist_mat <- dist(1 - sim_matrix, method = "euclidean")  # or "manhattan"
hc2 <- hclust(dist_mat, method = "average")   # usually no error
clusters2 <- cutree(hc2, h = 0.8)
cluster_df <- data.frame(
    ID = names(clusters2),
    cluster = clusters2
)

all_go_clustered <- left_join(all_go,cluster_df, by = "ID")
all_go_clustered %>%
    group_by(cluster) %>%
    summarise(
        n_datasets = n_distinct(dataset),
        mean_NES = mean(NES)) %>%
    arrange(desc(n_datasets))

cluster_summary <- all_go_clustered %>%
    group_by(cluster,dataset) %>%
    summarise(
        mean_NES = mean(NES),
        min_padj = min(p.adjust),
        n_terms = n(),
        .groups = "drop"
    )

library(dplyr)
cluster_labels <- all_go_clustered %>%
    group_by(cluster) %>%
    slice_min(p.adjust, n = 1 , with_ties = FALSE) %>%
    dplyr::select(cluster, ID,Description)

cluster_summary_labeled <- cluster_summary %>% left_join(cluster_labels, by= "cluster")

library(dplyr)
cluster_rank <- cluster_summary_labeled %>%
    group_by(Description) %>%
    summarise(n_datasets = n_distinct(dataset),
              mean_abs_NES = mean(abs(mean_NES))
    ) %>%
    arrange(desc(n_datasets),desc(mean_abs_NES))
# 15 top clusters that apeared in 3 or more than 3 datasets (reproducible)
top_clusters <- cluster_rank %>%
    filter(n_datasets >= 3) %>%
    slice_head( n = 15)
# all clusters that appeared in 3 or more than 3 datasets.
top_clusters_3 <- cluster_rank %>%
    filter(n_datasets >= 3)
saveRDS(top_clusters_3,"top_clusters_3.RDS")
saveRDS(top_clusters,"top_clusters.RDS")
saveRDS(cluster_summary_labeled,"cluster_summary_labeled.RDS")
