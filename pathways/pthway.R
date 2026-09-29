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

gsea_polo_keg <- gseKEGG(geneList = polokinase_list,
                         organism = "hsa",
                         minGSSize= 10,
                         pvalueCutoff = 0.05,
                         verbose = FALSE)
dotplot(gsea_polo,showCategory = 15)
dotplot(gsea_polo_keg, showCategory= 15)
