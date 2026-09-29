gene_df_1 <- bitr(multiple_unique$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- multiple_unique %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
    )
mapped_df1 <- mapped_df1 %>% filter(!is.na(stat))
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1, with_ties = FALSE) %>% ungroup()
multiple_list <- mapped_df1$stat
names(multiple_list) <- mapped_df1$ENTREZID
multiple_list <- sort(multiple_list,decreasing = TRUE) 
#I made a typo here , but has cleaned the swift before.
gsea_multi <- gseGO(
    geneList = multiple_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)
# save intial data 
saveRDS(gsea_multi, "gsea_multi.RDS")
gsea_multi <- readRDS("~/ALS/genes/GSEA/gsea_multi.RDS")
# now dotplot
k <- dotplot(gsea_multi,showCategory = 15, color = "NES")
x <- k + labs(title = "GSEA of GSE276214-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE276214_multi.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)

multi_go_res <- as.data.frame(gsea_multi)
multi_go_res <- multi_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(multi_go_res, "multi_GSEA_cleaned.RDS")
