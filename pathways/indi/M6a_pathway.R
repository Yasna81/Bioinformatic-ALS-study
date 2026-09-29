gene_df_1 <- bitr(m6a_res_clean$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- m6a_res_clean %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
    )
mapped_df1 <- mapped_df1 %>% filter(!is.na(stat))
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1, with_ties = FALSE) %>% ungroup()
m6a_list <- mapped_df1$stat
names(m6a_list) <- mapped_df1$ENTREZID
m6a_list <- sort(m6a_list,decreasing = TRUE) 
#I made a typo here , but has cleaned the swift before.
gsea_swift <- gseGO(
    geneList = m6a_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)
# save intial data 
saveRDS(gsea_swift, "gsea_M6A.RDS")
# now dotplot
library(ggplot2)
i <- dotplot(gsea_swift,showCategory = 15, color = "NES")
x <- i + labs(title = "GSEA of GSE242766-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE242766_M6a.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)
swift_go_res <- as.data.frame(gsea_swift)
swift_go_res <- swift_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(swift_go_res, "M6a_GSEA_cleaned.RDS")
