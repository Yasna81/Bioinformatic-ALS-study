multiple_counts <- read.table("~/ALS/data/Als/GSE276214_bfx731.hg38.e81.txt",
                header = TRUE,
                sep = "",
                stringsAsFactors = FALSE,
                check.names = FALSE)

samples_to_keep <- c("Geneid",
                     "N1_WT-axons",#GSM142
                     "wt_soma_JJ_n1",#GSM144
                     "wt_soma_LM_n2",#GSM146
                     "wt_axons_LM_n2",#GSM143
                     "wt_soma_JJ_n3",#GSM145
                     "wt_axons_JJ_n3",#GSM147
                     "axons_525L",
                     "m525L_somaJJ_n1",
                     "m525L_axonsLM_2",
                     "m525L_somaLM_2",
                     "m525L_somaJJ_n3",
                     "m525L_axonsJJn3"
)

multiple <- multiple_counts[,samples_to_keep]
#sum(duplicated(polokin_filter$Geneid)) was 0 so no duplicated gene!
rownames(multiple) <- multiple$Geneid
multiple$Geneid <- NULL
#time to merge 
multiple$WT_N1 <- multiple$`N1_WT-axons` + multiple$wt_soma_JJ_n1
multiple$WT_N2 <- multiple$wt_soma_LM_n2 + multiple$wt_axons_LM_n2
multiple$WT_N3 <- multiple$wt_soma_JJ_n3 + multiple$wt_axons_JJ_n3
multiple$mp525l_N1 <- multiple$axons_525L + multiple$m525L_somaJJ_n1
multiple$mp525l_N2 <- multiple$m525L_axonsLM_2 + multiple$m525L_somaLM_2
multiple$m525l_N3 <- multiple$m525L_somaJJ_n3 + multiple$m525L_axonsJJn3
#omiting previous counts :
multiple_1 <- multiple[,!grepl("axons|soma",colnames(multiple))]
saveRDS(multiple_1,"finalmulti_GSE276214_countmatrix.RDS")
#creating metadata :
sample_info <- data.frame(
    condition = factor(c("Control","Control","Control","Mutation","Mutation","Mutation")),
    row.names = colnames(multiple_1)
)
library(DESeq2)
sample_info$condition <- relevel(sample_info$condition, ref = "Control")
dds <-DESeqDataSetFromMatrix(countData = multiple_1,
                             colData = sample_info,
                             design = ~ condition)
#filter genes with low count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
#
norm_count_multi <- counts(dds,normalized = TRUE)
write.csv(norm_count_multi, "multitcount_GSE276214.csv", row.names = TRUE)
#
res <- results(dds,contrast = c("condition","Mutation","Control"))
res_ordered <- res[order(res$padj),]
res_df_multiple <- as.data.frame(res_ordered)

# did desq2 worked properly?
vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
plotDispEsts(dds)
#desq2 is handling noise :)
plot_colors <- c(rep("skyblue",3),rep("tomato",3))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info$condition,las=2,cex.axis=0.8)
#labeling 
library(org.Hs.eg.db)
res_df_multiple$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df_multiple),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
dup_symbols <- res_df_multiple$symbol[duplicated(res_df_multiple$symbol)]
# no gene id was duplicated only gene symbols were duplicated.
#most of the dups are na-dups. some are real symbols duplicated
####handeling dups
# restoring gene id col so that we can go for RRA
res_df_multiple$gene_id <- rownames(res_df_multiple)
library(dplyr)
multiple_unique <- res_df_multiple %>%
    group_by(symbol) %>%
    filter(!is.na(symbol)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(res_df_multiple %>% filter(is.na(symbol))) %>%
    arrange(desc(stat))
saveRDS(multiple_unique,"soma_axon.RDS")
