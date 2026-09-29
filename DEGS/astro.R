astro_counts <- read.delim("~/ALS/data/Als/GSE196219_RawCountsTable.txt")
samples_to_keep <- c("Gene.ID",
                     "Gene.Name",
                     "CO3.1.diff.1.S58.AACAGGTT.CTTGGTAT",#GSM2918
                     "CO3.1.diff.2.S59.GGTGAACC.TCCAACGC",#GSM2919
                     "CO3.1.diff.4.S60.CAACAATG.CCGTGAAG",#GSM2920
                     "X3.1.diff.1.S61.AAGATACT.ATGTAAGT",#GSM2912
                     "X3.1.diff.2.S62.GGAGCGTC.GCACGGAC",#GSM2913
                     "X3.1.diff.3.S63.ATGGCATG.GGTACCTT"#GSM2914
)
astro_fil<- astro_counts[,samples_to_keep]
sum(duplicated(astro_fil$Gene.ID)) #was 0 so no duplicated gene!
rownames(astro_fil) <- astro_fil$Gene.ID
#keeping gene name and ids to merge later 
info_frame <- astro_fil[,c("Gene.ID","Gene.Name")]
astro_fil$Gene.ID <- NULL
astro_fil$Gene.Name <- NULL
#creating metadata :
sample_info <- data.frame(
    condition = factor(c("Control","Control","Control","Mutation","Mutation","Mutation")),
    row.names = colnames(astro_fil)
)

sample_info$condition <- relevel(sample_info$condition, ref = "Control")
dds <-DESeqDataSetFromMatrix(countData = astro_fil,
                             colData = sample_info,
                             design = ~ condition)
#filter genes with low count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
#
norm_count_astro <- counts(dds,normalized = TRUE)
write.csv(norm_count_astro, "astro_count_GSE196219.csv", row.names = TRUE)
#

res <- results(dds,contrast = c("condition","Mutation","Control"))
res_ordered <- res[order(res$padj),]
res_df_astr <- as.data.frame(res_ordered)
#desq2
vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
plotDispEsts(dds)
#desq2 is handling noise :)
plot_colors <- c(rep("skyblue",3),rep("tomato",3))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info$condition,las=2,cex.axis=0.8)

#now merging 
library(data.table)
setDT(res_df_astr,keep.rownames = "gene_id")
info_frame$gene_id <- info_frame$Gene.ID
setDT(info_frame)
res_astro_merg <- merge(res_df_astr,info_frame,by="gene_id",all.x = TRUE)
rownames(res_astro_merg) <- res_astro_merg$Gene.ID
res_astro_merg$Gene.ID <- NULL
#initial data
saveRDS(res_astro_merg,"astro_res.RDS")


#cleaning initial data
res_astro_merg$gene_id <- rownames(res_astro_merg)
any(duplicated(res_astro_merg$gene_id)) # nothing duplicated in gene ids
any(duplicated(res_astro_merg$Gene.Name)) #
dup_symbols <- res_astro_merg$Gene.Name[duplicated(res_astro_merg$Gene.Name)]
res_astro_merg[res_astro_merg$Gene.Name %in% dup_symbols, c("Gene.Name")] #so many gene symbols are duplicated
# arranging and handeking the nas
library(dplyr)
res_astro_merg_clean <- res_astro_merg %>%
    group_by(Gene.Name) %>%
    filter(!is.na(Gene.Name)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(res_astro_merg%>% filter(is.na(Gene.Name))) %>%
    arrange(desc(stat))

saveRDS(res_astro_merg_clean,"res_final_clean.RDS")


