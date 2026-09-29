library(DESeq2)
polokin_27 <- read.delim("~/ALS/data/Als/polokinase/GSE272827_bfx1205.hg38.e92.tsv")
samples_to_keep <- c("Geneid",
                     "MN_WT_27.03.2019",#GSM142
                     "MN_WT_21.03.2019",#GSM144
                     "MN_WT_04.07.2019",#GSM146
                     "MN_P525L_27.03.2019",#GSM143
                     "MN_P525L_21.03.2019",#GSM145
                     "MN_P525L_04.07.2019"#GSM147
                     )
polokin_filter <- polokin_27[,samples_to_keep]
#sum(duplicated(polokin_filter$Geneid)) was 0 so no duplicated gene!
rownames(polokin_filter) <- polokin_filter$Geneid
polokin_filter$Geneid <- NULL
#creating metadata :
sample_info <- data.frame(
    condition = factor(c("Control","Control","Control","Mutation","Mutation","Mutation")),
    row.names = colnames(polokin_filter)
)

sample_info$condition <- relevel(sample_info$condition, ref = "Control")
dds <-DESeqDataSetFromMatrix(countData = polokin_filter,
                             colData = sample_info,
                             design = ~ condition)
#filter genes with low count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
#for accuarey :
norm_count <- counts(dds,normalized = TRUE)
saveRDS(norm_count, "polocount_GSE272827.RDS")
#
res <- results(dds,contrast = c("condition","Mutation","Control"))
res_ordered <- res[order(res$padj),]
res_df <- as.data.frame(res_ordered)

# did desq2 worked properly?
vsd <- vst(dds,blind = FALSE)
plotPCA(vsd,intgroup="condition")
plotDispEsts(dds)
#desq2 is handling noise :)
plot_colors <- c(rep("skyblue",3),rep("tomato",3))
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main="Raw Counts",col= plot_colors,border = "dodgerblue",names=sample_info$condition,las =2,cex.axis=0.8)
boxplot(assay(vsd),main="Normalized (vst) Counts",col=plot_colors,border = "dodgerblue",names= sample_info$condition,las=2,cex.axis=0.8)
#up and down :
#volcano plot :
library(org.Hs.eg.db)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
library(ggplot2)
library(ggrepel)
res_df$label <- ifelse(!is.na(res_df$symbol) & 
                           res_df$padj < 1e-5 &
                           abs(res_df$log2FoldChange) > 3,
                       res_df$symbol, "")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    # 1. Background points (Grey for non-significant)
    geom_point(data = subset(res_df, padj >= 0.05 | abs(log2FoldChange) <= 1), 
               color = "grey80", alpha = 0.4, size = 1) +
    
    # 2. Significant Up-regulated (Red)
    geom_point(data = subset(res_df, padj < 0.05 & log2FoldChange > 1), 
               color = "firebrick", alpha = 0.6, size = 1.5) +
    
    # 3. Significant Down-regulated (Blue)
    geom_point(data = subset(res_df, padj < 0.05 & log2FoldChange < -1), 
               color = "dodgerblue4", alpha = 0.6, size = 1.5) +
    
    # 4. Add the labels using ggrepel (won't overlap)
    geom_text_repel(aes(label = label),
                    size = 3.5,
                    fontface = "italic",
                    box.padding = 0.5, 
                    point.padding = 0.3,
                    max.overlaps = 15,
                    segment.color = 'grey50') +
    
    # 5. Reference lines for the thresholds
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    
    # 6. Theme and Labels
    theme_minimal() +
    labs(title = "Differential Expression: Mutation vs Control",
         subtitle = "Red = Up-regulated, Blue = Down-regulated",
         x = "log2(Fold Change)",
         y = "-log10(Adjusted P-value)") +
    theme(panel.grid.minor = element_blank())
#saving the initial version :
write.csv(res_df,"polokinase_results.csv")
polokinas_res <- read.csv("genes/res/polokinase_results.csv")
#checking if any duplicants in gene ids :
any(duplicated(polokinas_res$X)) # nothing duplicated in gene ids
any(duplicated(polokinas_res$symbol)) #
dup_symbols <- polokinas_res$symbol[duplicated(polokinas_res$symbol)]
polokinas_res[polokinas_res$symbol %in% dup_symbols, c("symbol")] # beside NAs we have some duplicated gene symbols
# arranging and handeking the nas
library(dplyr)
polokinas_res_clean <- polokinas_res %>%
    group_by(symbol) %>%
    filter(!is.na(symbol)) %>%
    slice_max(stat,n=1,with_ties = FALSE) %>%
    ungroup() %>%
    bind_rows(polokinas_res %>% filter(is.na(symbol))) %>%
    arrange(desc(stat))

saveRDS(polokinas_res_clean,"polokinase_clean.RDS")
