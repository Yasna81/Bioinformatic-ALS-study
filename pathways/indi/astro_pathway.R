library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
#polokinase GSEA
gene_df_1 <- bitr(res_astro_merg_clean$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- res_astro_merg_clean %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
    )
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1, with_ties = FALSE) %>% ungroup()
astor_list <- mapped_df1$stat
names(astor_list) <- mapped_df1$ENTREZID
astor_list <- sort(astor_list,decreasing = TRUE)

gsea_astor <- gseGO(
    geneList = astor_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)
# save intial data 
saveRDS(gsea_astor, "gsea_astor.RDS")
# now dotplot
library(ggplot2)
p <- dotplot(gsea_astor,showCategory = 15, color = "NES")
x <- p + labs(title = "GSEA of GSE196219-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE196219_astro.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)

astro_go_res <- as.data.frame(gsea_astor)
astro_go_res <- astro_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(astro_go_res, "astro_GSEA_cleaned.RDS")

