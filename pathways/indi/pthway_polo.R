library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
#polokinase GSEA
gene_df_1 <- bitr(polokinas_res_clean$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- polokinas_res_clean %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
)
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1) %>% ungroup()
polokinase_list <- mapped_df1$stat
names(polokinase_list) <- mapped_df1$ENTREZID
polokinase_list <- sort(polokinase_list,decreasing = TRUE)

gsea_polo <- gseGO(
    geneList = polokinase_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)
# save intial data 
saveRDS(gsea_polo, "gsea_polo.RDS")
# now dotplot
h <- dotplot(gsea_polo,showCategory = 15, color = "NES")
x <- h + labs(title = "GSEA of GSE272827-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE272827_polo.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)

polo_go_res <- as.data.frame(gsea_polo)
polo_go_res <- polo_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(polo_go_res, "polo_GSEA_cleaned.RDS")
