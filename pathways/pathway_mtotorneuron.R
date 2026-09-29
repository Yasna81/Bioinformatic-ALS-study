polo_gse <- readRDS("~/ALS/genes/GSEA/clean/polo_GSEA_cleaned.RDS")
swift_gse <- readRDS("~/ALS/genes/GSEA/clean/swift_GSEA_cleaned.RDS")
M6a_gse <- readRDS("~/ALS/genes/GSEA/clean/M6a_GSEA_cleaned.RDS")
fus_gse <-readRDS("~/ALS/genes/GSEA/clean/fus_GSEA_cleaned.RDS")
multi_gse <- readRDS("~/ALS/genes/GSEA/clean/multi_GSEA_cleaned.RDS")
polo_sig <- polo_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D1")
swift_sig <- swift_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D2")
M6a_sig <- M6a_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D3")
fus_sig <- fus_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D4")
multi_sig <- multi_gse %>%
    filter(p.adjust < 0.05) %>%
    mutate(dataset = "D5")
all_go <- bind_rows(multi_sig,fus_sig,M6a_sig,swift_sig,polo_sig)
go_ids <- unique(all_go$ID)
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db',ont="BP")
sim_matrix <- mgoSim(go_ids,
                     go_ids,
                     semData = hsGO,
                     measure = "Wang",
                     combine = NULL)

dist_matrix <- as.dist(1 - sim_matrix)
hc <- hclust(dist_matrix, method = "average")
clusters <- cutree(hc,h= 0.8)
cluster_df <- data.frame(
    ID = names(clusters),
    cluster = clusters
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
    select(cluster, ID,Description)
cluster_summary_labeled <- cluster_summary %>% left_join(cluster_labels, by= "cluster")
# all significant pathwys with go terms
library(ggplot2)
rep_patways <-ggplot(cluster_summary_labeled,
       aes(x= dataset,
           y = Description,
           size = n_terms,
           color=mean_NES)) +
    geom_point() +
    scale_color_gradient2(low = "blue",
                          mid = "white",
                          high = "red") +
    theme_minimal()
ggsave("all_mt_rep_pathways.jpeg",rep_patways, width = 10,height = 8,dpi = 300)

library(dplyr)
cluster_rank <- cluster_summary_labeled %>%
    group_by(Description) %>%
    summarise(n_datasets = n_distinct(dataset),
              mean_abs_NES = mean(abs(mean_NES))
              ) %>%
    arrange(desc(n_datasets),desc(mean_abs_NES))
# 15 top clusters that apeared in 4 or more than 4 datasets (reproducible)
top_clusters <- cluster_rank %>%
    filter(n_datasets >= 4) %>%
    slice_head( n = 15)
top_names <- top_clusters$Description
plot_df <- cluster_summary_labeled %>%
    filter(Description %in% top_names)
#heatmap of 15
library(pheatmap)
library(tidyr)
head_df <- plot_df %>%
    select(Description,dataset, mean_NES) %>%
    pivot_wider(names_from = dataset,
                values_from = mean_NES)
mat <- as.matrix(head_df[,-1])
rownames(mat) <- head_df$Description
library(RColorBrewer)
library(viridis)
my_colors <- viridis(100)
heatmap_mt <- pheatmap(mat,
         color = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(100),
         border_color = "white",
         labels_col = c("GSE272827","GSE229095","GSE242766","GSE94888","GSE276214"),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         angle_col = 45,
         breaks = seq(-2.5,2.5,length.out= 100),
         fontsize_row = 13,
         fontsize_col = 12,
         row_names_side = "left",
         width = 10,
         height =  6,
         dpi= 300,
         filename = "plot.jpeg")

# info

cluster_rank %>%
    count(n_datasets)

library(ggplot2)
library(forcats) # for ordering

# Ensure Description is ordered by mean NES or reproducibility
plot_df <- plot_df %>%
    mutate(Description = fct_reorder(Description, mean_NES, .fun = mean, .desc = TRUE))

# Bubble (Emma-style) plot
dot <- ggplot(plot_df, aes(x = dataset,
                    y = Description,
                    size = n_terms, # size = number of GO terms in cluster
                    color = mean_NES)) + # color = NES
    geom_point(alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(3,10)) +
    theme_minimal(base_size = 11) +
    theme(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(
        x = "Dataset",
        y = "Biological Theme",
        color = "Mean NES",
        size = "Number of GO terms"
    ) +
    scale_x_discrete(labels = c("D1"="GSE272827",
                                "D2" = "GSE229095",
                                "D3"="GSE242766",
                                "D4" ="GSE94888",
                                "D5"= "GSE276214"))

ggsave("15-mt-dotplot.jpeg",dot, width = 10,height = 8,dpi = 300)

saveRDS(top_clusters,"top_clusters.RDS")
#reprducabiliy
saveRDS(sim_matrix,"sim_matrix.RDS")
saveRDS(all_go, "all_go.RDS")
saveRDS(cluster_summary,"cluster_summary_2.RDS")
saveRDS(cluster_summary_labeled,"cluster_summary_rank2.RDS")
