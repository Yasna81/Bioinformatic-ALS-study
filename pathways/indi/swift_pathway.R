#swift 
#swift GSEA
gene_df_1 <- bitr(res_df_2_clean$gene_id,
                  fromType = "ENSEMBL",
                  toType= "ENTREZID",
                  OrgDb = org.Hs.eg.db)
mapped_df1 <- res_df_2_clean %>%
    inner_join(gene_df_1, by = c("gene_id"="ENSEMBL")
    )
mapped_df1 <- mapped_df1 %>% filter(!is.na(stat))
mapped_df1 <- mapped_df1 %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n=1, with_ties = FALSE) %>% ungroup()
swift_list <- mapped_df1$stat
names(swift_list) <- mapped_df1$ENTREZID
swift_list <- sort(swift_list,decreasing = TRUE) 

gsea_swift <- gseGO(
    geneList = swift_list,
    OrgDb = org.Hs.eg.db,
    ont= "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
)

# save intial data 
saveRDS(gsea_swift, "gsea_swift.RDS")
# now dotplot
n <- dotplot(gsea_swift,showCategory = 15, color = "NES")
x <- n + labs(title = "GSEA of GSE229095-Top Enriched GO Terms") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14))
ggsave("GSE229095_swift.jpeg",
       plot = x ,
       width = 10,
       height = 6,
       dpi = 300)

swift_go_res <- as.data.frame(gsea_swift)
swift_go_res <- swift_go_res %>% select(ID,Description,NES,p.adjust) %>% arrange(desc(NES))
saveRDS(swift_go_res, "swift_GSEA_cleaned.RDS")
