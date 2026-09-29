control_1 <- read.delim("~/ALS/data/Als/swift/GSM7152417_6_quant.sf")
control_2 <-read.delim("~/ALS/data/Als/swift/GSM8632052_15_quant.sf")
mutation_1 <- read.delim("~/ALS/data/Als/swift/GSM7152418_7_quant.sf")
mutation_2 <- read.delim("~/ALS/data/Als/swift/GSM7152423_12_quant.sf")
#consistency in gene names :
library(data.table)
library(org.Hs.eg.db)
tx2gene <- select(org.Hs.eg.db,
                  keys = keys(org.Hs.eg.db, keytype = "ENSEMBLTRANS"),
                  columns = "ENSEMBL",
                  keytype = "ENSEMBLTRANS")
setDT(tx2gene)
setnames(tx2gene,c("ENST","ENSG"))
#
files <- list(control_1,control_2,mutation_1,mutation_2)
sample_name <- c("Control_1","Control_2","Mutation_1","Mutation_2")
data_list <- lapply(files,function(df) {
    dt <- as.data.table(df)
    return(dt[,.(Name,NumReads)])
})
combined_dt <- Reduce(function(x,y) merge(x,y,by ="Name"), data_list)
setnames(combined_dt, c("ENST", sample_name))
combined_dt[,ENST_base := sub("\\..*","",ENST)]
combined_dt <- merge(combined_dt,tx2gene,by.x = "ENST_base",by.y="ENST")
final_count_dt <- combined_dt[,lapply(.SD,sum),by= ENSG, .SDcols = sample_name]
# no NA
final_count_dt<- na.omit(final_count_dt, cols= "ENSG")
final_counts <- as.data.frame(final_count_dt)
rownames(final_counts) <- final_count_dt$ENSG
final_counts$ENSG <-NULL
final_counts <- round(as.matrix(final_counts))
saveRDS(final_counts ,"final_swift_countmatrix.RDS")
#count matrix ready, now we head to desq2
sample_info_2 <- data.frame(condition= factor(c("Control","Control","Mutation","Mutation"),
                                              levels = c("Control","Mutation"))
                            )
rownames(sample_info_2) <- colnames(final_counts)
sample_info_2$condition <- relevel(sample_info_2$condition, ref = "Control")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData = final_counts,
                             colData = sample_info_2,
                             design = ~ condition)
vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
#kind of fine.
dds <- DESeq(dds)
#check 
norm_count_swift <- counts(dds,normalized = TRUE)
write.csv(norm_count_swift, "swiftcount_GSE229095.csv", row.names = TRUE)
#
res_2 <- results(dds)
res_df_2 <- as.data.frame(res_2)
res_df_2$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(res_df_2),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
# normalization 
plot_colors <- c(rep("skyblue",2),rep("tomato",2))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info_2$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info_2$condition,las=2,cex.axis=0.8)
#up and down :
res_df_2$gene_id <- rownames(res_df_2)

any(duplicated(res_df_2$gene_id)) # nothing duplicated in gene ids
any(duplicated(res_df_2$symbol)) #
dup_symbols <- res_df_2$symbol[duplicated(res_df_2$symbol)]
res_df_2[res_df_2$symbol %in% dup_symbols, c("symbol")] #so many gene symbols are duplicated
# arranging and handeking the nas
library(dplyr)
res_df_2_clean <- res_df_2 %>%
    group_by(symbol) %>%
    filter(!is.na(symbol)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(res_df_2 %>% filter(is.na(symbol))) %>%
    arrange(desc(stat))

saveRDS(res_df_2_clean,"swift_clean.RDS")

