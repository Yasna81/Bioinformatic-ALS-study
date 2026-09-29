polokin_count <- readRDS("counts/polocount_GSE272827.RDS")
write.csv(polokin_count, "polokin_count_GSE272827.csv", row.names = TRUE)
#------------------------------res conversion--------------------------
res_polo <- readRDS("~/ALS/genes/arranged/polokinase_clean.RDS")
res_polo <- res_polo[order(res_polo$padj),]
res_polo$direction <-ifelse(res_polo$padj < 0.05 & res_polo$log2FoldChange > 0, "UP",
                            ifelse(res_polo$padj < 0.05 & res_polo$log2FoldChange < 0 , "Down", "NS"))
write.csv(res_polo, "polokin_res.csv", row.names = FALSE)
#---------------------------------------------------------
res_swift <- readRDS("~/ALS/genes/arranged/swift_clean.RDS")
res_swift <- res_swift[order(res_swift$padj),]
res_swift$direction <-ifelse(res_swift$padj < 0.05 & res_swift$log2FoldChange > 0, "UP",
                            ifelse(res_swift$padj < 0.05 & res_swift$log2FoldChange < 0 , "Down", "NS"))
write.csv(res_swift, "swift_res.csv", row.names = FALSE)
#-------------------------------
res_m6a <- readRDS("~/ALS/genes/arranged/M6a_clean.RDS")
res_fus <- readRDS("~/ALS/genes/arranged/fus_final_clean.RDS")
res_multi<- readRDS("~/ALS/genes/arranged/soma_axon.RDS")
res_astro <- readRDS("~/ALS/genes/arranged/astrores_final_clean.RDS")
