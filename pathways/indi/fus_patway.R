gene_df_1 <- bitr(res_fus_clean$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- res_fus_clean %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
    )
mapped_df1 <- mapped_df1 %>% filter(!is.na(stat))
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1, with_ties = FALSE) %>% ungroup()
fus_list <- mapped_df1$stat
names(fus_list) <- mapped_df1$ENTREZID
fus_list <- sort(fus_list,decreasing = TRUE) 
#I made a typo here , but has cleaned the swift before.
gsea_fus <- gseGO(
    geneList = fus_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)
# save intial data 
saveRDS(gsea_fus, "gsea_fus.RDS")
# now dotplot
library(ggplot2)
d <- dotplot(gsea_fus,showCategory = 15, color = "NES")
x <- d + labs(title = "GSEA of GSE94888-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE94888_fus.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)

fus_go_res <- as.data.frame(gsea_fus)
fus_go_res <- fus_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(fus_go_res, "fus_GSEA_cleaned.RDS")
