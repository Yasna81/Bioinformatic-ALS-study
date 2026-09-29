astro_gse <- readRDS("~/ALS/genes/GSEA/gsea_astor.RDS")
astro_gse_sig <- as.data.frame(astro_gse) %>% filter(p.adjust < 0.05)
leading_genes <- unique(unlist(strsplit(astro_gse_sig$core_enrichment, "/")))
library(STRINGdb)
string_db <- STRINGdb$new(
    version="11.5",
    species = 9606,
    score_threshold = 700
)
library(GOSemSim)
library(org.Hs.eg.db)
library(clusterProfiler)
gene_symbols_2 <- bitr(
    leading_genes,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
)
mapped <- string_db$map(gene_symbols_2,
                        "SYMBOL",
                        removeUnmappedRows = TRUE)
library(DESeq2)

ppi <- string_db$get_interactions(mapped$STRING_id)

library(igraph)
g <- graph_from_data_frame(ppi[,1:2],directed = FALSE)
deg <- degree(g)
bet <- betweenness(g)
close <- closeness(g)
library(dplyr)
hub_table_astro <- data.frame(gene = names(deg),
                        degree = deg,
                        betweenness = bet,
                        closeness = close)
hub_table_astro <- hub_table_astro %>% arrange(desc(degree))
head(hub_table_astro,20)
hub_astro_20 <- slice_head(hub_table_astro, n = 20)
hub_astro_symbols <- merge(hub_astro_20, mapped[,c("STRING_id","SYMBOL")],
                          by.x = "gene",by.y = "STRING_id")
#Load packages
library(igraph)
library(tidygraph)
library(ggraph)
library(viridis)
library(dplyr)

# --------------------------
# 1. Build graph from edge list
# --------------------------
# Assuming ppi has columns: from, to (use gene symbols or STRING IDs)
g <- graph_from_data_frame(ppi, directed = FALSE)

# --------------------------
# 2. Add node attributes from hub_astro
# --------------------------
# hub_astro columns: protein_id, symbol, degree, betweenness
node_incidices <- match(V(g)$name, hub_astro_symbols$gene)
V(g)$degree <- hub_astro_symbols$degree[node_incidices]
V(g)$betweenness <- hub_astro_symbols$betweenness[node_incidices]
V(g)$display_lable <- hub_astro_symbols$SYMBOL[node_incidices]
# --------------------------
# 3. Convert to tidygraph for ggraph
# --------------------------
graph_tbl <- as_tbl_graph(g)

# --------------------------
# 4. Plot the network
# --------------------------

vis_tbl <- graph_tbl %>% 
    activate(nodes) %>% 
    filter(degree > 20) # Lower this number if the plot is still too empty

n <- ggraph(vis_tbl, layout = "stress") + 
    geom_edge_link(alpha = 0.2, color = "grey80") + # Increased alpha since fewer edges
    geom_node_point(aes(size = degree, color = betweenness), alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", trans = "sqrt") + # sqrt is safer than log10
    geom_node_text(aes(label = display_lable), 
                   repel = TRUE, 
                   size = 4, 
                   fontface = "bold",
                   max.overlaps = 100) + 
    theme_void() +
    labs(title = "Core Protein Interaction Hubs in Astrocytes") 
    ggsave(filename = "PPI_network_astor.jpeg",
           plot = last_plot(),
           device = "jpeg",
           width = 10,
           height = 8,
           dpi = 300,
           bg = "white")
# --------------------------
# 5. Save figure for publication
# --------------------------
ggsave("PPI_network_astor.png", width = 12, height = 8, units =  "in",dpi= 600, bg = "white")
saveRDS(ppi_2,"ppi_astro.RDS")
saveRDS(mapped, "mapped_astro.RDS")
saveRDS(hub_astro_symbols,"hub_astro_symbols.RDS")
