counts <- read.delim("~/ALS/data/Als/GSE242766_counts_table.tsv")
samples_to_keep <- c("X",
                     "DOXY._ARS._shMETTL3_FUS.WT_INP_REP1_S3",#GSM0162
                     "DOXY._ARS._shMETTL3_FUS.WT_INP_REP2_S21",#GSM0163
                     "DOXY._ARS._shMETTL3_FUS.WT_INP_REP1_S26",#GSM170
                     "DOXY._ARS._shMETTL3_FUS.WT_INP_REP2_S13",#GSM171
                     "DOXY._ARS._shMETTL3_FUS.P525L_INP_REP1_S11",#GSM0158
                     "DOXY._ARS._shMETTL3_FUS.P525L_INP_REP2_S20",#GSM0159
                     "DOXY._ARS._shMETTL3_FUS.P525L_INP_REP1_S43",#GSM166
                     "DOXY._ARS._shMETTL3_FUS.P525L_INP_REP2_S45"#GSM167
                     
)
M6a_counts <- counts[,samples_to_keep]
#sum(duplicated(M6a_counts$X)) no duplication 
rownames(M6a_counts) <- M6a_counts$X
M6a_counts$X <- NULL
#creating metadata :
sample_info <- data.frame(
    condition = factor(c("Control","Control","Control","Control","Mutation","Mutation","Mutation","Mutation")),
    row.names = colnames(M6a_counts)
)
library(DESeq2)
sample_info$condition <- relevel(sample_info$condition, ref = "Control")
dds <-DESeqDataSetFromMatrix(countData = M6a_counts,
                             colData = sample_info,
                             design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
#
norm_count_m6a <- counts(dds,normalized = TRUE)
write.csv(norm_count_m6a, "m6atcount_GSE242766.csv", row.names = TRUE)
#
res <- results(dds,contrast = c("condition","Mutation","Control"))
res_ordered <- res[order(res$padj),]
res_df <- as.data.frame(res_ordered)
#lets  get symbols
library(org.Hs.eg.db)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
# did desq2 worked properly?
vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
plotDispEsts(dds)
#desq2 is handling noise :)
plot_colors <- c(rep("skyblue",4),rep("tomato",4))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info$condition,las=2,cex.axis=0.8)
# saving initial data 
saveRDS(res_df,"res_M6a.RDS")
#cleaning initial data
m6a_res <- readRDS("genes/res/res_M6a.RDS")
m6a_res$gene_id <- rownames(m6a_res)
any(duplicated(m6a_res$gene_id)) # nothing duplicated in gene ids
any(duplicated(m6a_res$symbol)) #
dup_symbols <- m6a_res$symbol[duplicated(m6a_res$symbol)]
m6a_res[m6a_res$symbol %in% dup_symbols, c("symbol")] #so many gene symbols are duplicated
# arranging and handeking the nas
library(dplyr)
m6a_res_clean <- m6a_res %>%
    group_by(symbol) %>%
    filter(!is.na(symbol)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(m6a_res %>% filter(is.na(symbol))) %>%
    arrange(desc(stat))

saveRDS(m6a_res_clean,"M6a_clean.RDS")
