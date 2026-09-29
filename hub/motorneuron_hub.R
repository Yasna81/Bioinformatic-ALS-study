polo_gse <- readRDS("~/ALS/genes/GSEA/gsea_polo.RDS")
swift_gse <- readRDS("~/ALS/genes/GSEA/gsea_swift.RDS")
m6a_gse <- readRDS("~/ALS/genes/GSEA/gsea_M6A.RDS")
fus_gse <- readRDS("~/ALS/genes/GSEA/gsea_fus.RDS")
multi_gse <- readRDS("~/ALS/genes/GSEA/gsea_multi.RDS")
library(dplyr)
polo_gse_sig <- as.data.frame(polo_gse) %>% filter(p.adjust < 0.05)
swift_gse_sig <- as.data.frame(swift_gse) %>% filter(p.adjust < 0.05)
m6a_gse_sig <- as.data.frame(m6a_gse) %>% filter(p.adjust < 0.05)
fus_gse_sig <- as.data.frame(fus_gse) %>% filter(p.adjust < 0.05)
multi_gse_sig <- as.data.frame(multi_gse) %>% filter(p.adjust < 0.05)
sig_list <- list(polo_gse_sig,swift_gse_sig,m6a_gse_sig,fus_gse_sig,multi_gse_sig)
go_ids_list <- lapply(sig_list, function(g) g$ID)
all_ids <- unlist(go_ids_list)
go_freq <- table(all_ids)
shared_go <- names(go_freq[go_freq >=3])
extract_leading <- function(df){
    shared_df <- df %>% filter(ID %in% shared_go)
    genes <- unique(unlist(strsplit(shared_df$core_enrichment,"/")))
    return(genes)
}
        
gene_list <- lapply(sig_list,extract_leading)
all_genes <- unlist(gene_list)
gene_freq <- table(all_genes)
consensus_genes <- names(gene_freq[gene_freq>=3])
length(consensus_genes)
#160
library(STRINGdb)
options(timeout = 2000)
string_db <- STRINGdb$new(
    version = "11.5",
    species = 9606,
    score_threshold = 700,
    input_directory=""
)
library(org.Hs.eg.db)
library(clusterProfiler)
gene_symbols <- bitr(
    consensus_genes,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
)
mapped <- string_db$map(gene_symbols,
                        "SYMBOL",
                        removeUnmappedRows = TRUE)

ppi_2 <- string_db$get_interactions(mapped$STRING_id)

library(igraph)
g <- graph_from_data_frame(ppi_2[,1:2],directed = FALSE)
deg <- degree(g)
bet <- betweenness(g)
close <- closeness(g)
library(dplyr)
hub_table <- data.frame(gene = names(deg),
                        degree = deg,
                        betweenness = bet,
                        closeness = close)
hub_table <- hub_table %>% arrange(desc(degree))
head(hub_table,20)
hub_multi <- slice_head(hub_table, n = 20)
hub_with_symbols <- merge(hub_table, mapped[,c("STRING_id","SYMBOL")],
                          by.x = "gene",by.y = "STRING_id")
hub_ids <- hub_with_symbols$gene


# Load packages
library(igraph)
library(tidygraph)
library(ggraph)
library(viridis)
library(dplyr)

# --------------------------
# 1. Build graph from edge list
# --------------------------
# Assuming ppi has columns: from, to (use gene symbols or STRING IDs)
g <- graph_from_data_frame(ppi_2, directed = FALSE)

# --------------------------
# 2. Add node attributes from hub_astro
# --------------------------
# hub_astro columns: protein_id, symbol, degree, betweenness
node_incidices <- match(V(g)$name, hub_with_symbols$gene)
V(g)$degree <- hub_with_symbols$degree[node_incidices]
V(g)$betweenness <- hub_with_symbols$betweenness[node_incidices]
V(g)$display_lable <- hub_with_symbols$SYMBOL[node_incidices]
# --------------------------
# 3. Convert to tidygraph for ggraph
# --------------------------
graph_tbl <- as_tbl_graph(g)

# --------------------------
# 4. Plot the network
# --------------------------
s <- ggraph(graph_tbl, layout = "nicely") + # Force-directed layout
    # Edges
    geom_edge_link(alpha = 0.3, color = "grey70") +
    # Nodes
    geom_node_point(aes(size = degree, 
                        color = betweenness)) +
    # Labels (optional: repel to avoid overlap)
    geom_node_text(aes(label = display_lable),
                   repel = TRUE,
                   size = 3,
                   max.overlaps = 20,
                   point.padding = unit(0.8,"lines"),
                   fontface = "bold") +
    # Scales
    scale_size_continuous(range = c(3, 10)) + # adjust node sizes
    scale_color_viridis(option = "plasma",trans = "log10",na.value = "grey50") + # betweenness color
    scale_x_continuous(expand= expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    # Clean theme
    theme_void() +
    labs(title = "Core Protein Interaction Hubs in Motorneurons") +
    coord_fixed() +
    labs(color = "Betweenness Centrality", size = "Degree") 

# --------------------------
# 5. Save figure for publication
# --------------------------
ggsave(filename = "PPI_network_motorneuron.jpeg",
       plot = last_plot(),
       device = "jpeg",
       width = 12,
       height = 8,
       dpi = 300,
       bg = "white")
ggsave("PPI_network_multi.png", width = 12, height = 8, units =  "in",dpi= 600, bg = "white")
saveRDS(ppi_2, "ppi_multi.RDS")
saveRDS(mapped , "mapped_MULTI.RDS")
saveRDS(hub_with_symbols, "hub_symbol_mnorun.RDS")
