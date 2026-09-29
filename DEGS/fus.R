fus_counts <- read.delim("~/ALS/data/Als/GSE94888_raw_counts_GRCh38.p13_NCBI.tsv")
samples_to_keep <- c("GeneID",
                     "GSM2490554",
                     "GSM2490555",
                     "GSM2490556",
                     "GSM2490557",
                     "GSM2490558",
                     "GSM2490559"
)
fus_counts_fil <- fus_counts[,samples_to_keep]

# converting to gene ids .....
library(data.table)
library(org.Hs.eg.db)

# 1. Convert to data.table
setDT(fus_counts_fil)

# 2. Create the translation map (Numeric -> ENSG)
# We use the unique IDs from your column to save time/RAM
mapping <- select(org.Hs.eg.db, 
                  keys = as.character(unique(fus_counts_fil$GeneID)), 
                  columns = "ENSEMBL", 
                  keytype = "ENTREZID")
setDT(mapping)
fus_counts_fil$GeneID <- as.character(fus_counts_fil$GeneID)
# 3. Merge the mapping into your study 3 data
# 'by.x' is your current column name, 'by.y' is the map column name
df3 <- merge(fus_counts_fil, mapping, by.x = "GeneID", by.y = "ENTREZID")

# 4. Remove rows that didn't find a match (NAs) and the old ID column
df3 <- df3[!is.na(ENSEMBL)]
df3[, GeneID := NULL] 

# 5. Collapse: Sum counts if multiple Entrez IDs map to one ENSG
# This is fast and handles the 20,000+ rows easily
sample_names <- colnames(df3)[colnames(df3) != "ENSEMBL"]
final_counts_3 <- df3[, lapply(.SD, sum), by = ENSEMBL, .SDcols = sample_names]
final_df_3 <- as.data.frame(final_counts_3)
rownames(final_df_3) <- final_df_3$ENSEMBL
final_df_3$ENSEMBL <- NULL
final_df_3 <- round(as.matrix(final_df_3))
saveRDS(final_df_3, "final_fus_countmatrix.RDS")

#creating metadata :
sample_info <- data.frame(
    condition = factor(c("Control","Control","Control","Mutation","Mutation","Mutation")),
    row.names = colnames(final_df_3)
)

sample_info$condition <- relevel(sample_info$condition, ref = "Control")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData = final_df_3,
                             colData = sample_info,
                             design = ~ condition)
#filter genes with low count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
#
norm_count_fus <- counts(dds,normalized = TRUE)
write.csv(norm_count_fus, "fusatcount_GSE94888.csv", row.names = TRUE)
#
res <- results(dds,contrast = c("condition","Mutation","Control"))
res_ordered <- res[order(res$padj),]
res_df_fus <- as.data.frame(res_ordered)
#

vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
plotDispEsts(dds)
#desq2 is handling noise :)
plot_colors <- c(rep("skyblue",3),rep("tomato",3))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info$condition,las=2,cex.axis=0.8)
#up and down :
library(org.Hs.eg.db)
res_df_fus$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df_fus),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
#saving initial data
saveRDS(res_df_fus,"res_df_fus.RDS")
#cleaning initial data
res_df_fus$gene_id <- rownames(res_df_fus)
any(duplicated(res_df_fus$gene_id)) # nothing duplicated in gene ids
any(duplicated(res_df_fus$symbol)) #
dup_symbols <- res_df_fus$symbol[duplicated(res_df_fus$symbol)]
res_df_fus[res_df_fus$symbol %in% dup_symbols, c("symbol")] #so many gene symbols are duplicated
# arranging and handeking the nas
library(dplyr)
res_fus_clean <- res_df_fus %>%
    group_by(symbol) %>%
    filter(!is.na(symbol)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(res_df_fus %>% filter(is.na(symbol))) %>%
    arrange(desc(stat))

saveRDS(res_fus_clean,"fus_final_clean.RDS")

