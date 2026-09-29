library(ggrepel)
library(ggplot2)

polo <- readRDS("~/ALS/genes/arranged/polokinase_clean.RDS")
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
polo <- polo %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
polo <- polo %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- polo %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
polo_plt <- ggplot(polo, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = symbol),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    coord_cartesian(
        xlim = c(-10,10)
    ) +
    
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE272827",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )

ggsave("GSE272827_polo.jpeg",
       plot = polo_plt,
       width = 10,
       height = 6,
       dpi = 300)

#----------------------------------------------------------
swift <- readRDS("~/ALS/genes/arranged/swift_clean.RDS")
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
swift <- swift %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
swift <- swift %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- swift %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
swift_plt <-ggplot(swift, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = symbol),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    coord_cartesian(
        xlim = c(-10,10)
    ) +
    
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE229095",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )

ggsave("GSE229095_swift.jpeg",
       plot = swift_plt,
       width = 10,
       height = 6,
       dpi = 300)


#----------------------------------------------------------
m6a <- readRDS("~/ALS/genes/arranged/M6a_clean.RDS")
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
m6a <- m6a %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
m6a <- m6a %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- m6a %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
m6a_plt <- ggplot(m6a, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = symbol),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    coord_cartesian(
        xlim = c(-10,10)
    ) +
    
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE242766",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )

ggsave("GSE242766_M6a.jpeg",
       plot = m6a_plt,
       width = 10,
       height = 6,
       dpi = 300)

#----------------------------------------------------------------
fus <- readRDS("~/ALS/genes/arranged/fus_final_clean.RDS")
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
fus <- fus %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
fus <- fus %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- fus %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
fus_plt <- ggplot(fus, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = symbol),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    coord_cartesian(
        xlim = c(-10,10) 
        ) +
    
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE94888",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )

ggsave("GSE94888_fus.jpeg",
       plot = fus_plt,
       width = 10,
       height = 6,
       dpi = 300)
#----------------------------------------------------
multi <- readRDS("~/ALS/genes/arranged/soma_axon.RDS")
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
multi <- multi %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
multi <- multi %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- multi %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
multi_plt <- ggplot(multi, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = symbol),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    # Symmetric X-axis
    coord_cartesian(xlim=c(-10,10)) +
    
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE276214",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )
ggsave("GSE276214_multi.jpeg",
       plot = multi_plt,
       width = 10,
       height = 6,
       dpi = 300)
#------------------------------------------------------------
astro <- readRDS("~/ALS/genes/arranged/astrores_final_clean.RDS")
library(ggplot2)
library(ggrepel)
library(dplyr)
# ===============================
# 1. Define statistical thresholds
# ===============================
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Avoid Inf values
astro <- astro %>%
    mutate(
        padj = ifelse(padj == 0, 1e-300, padj),
        neg_log10_padj = -log10(padj)
    )

# ==================================
# 2. Classify genes biologically
# ==================================
astro<- astro %>%
    mutate(
        significance = case_when(
            padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Upregulated",
            padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
            TRUE ~ "Not Significant"
        )
    )

# ==================================
# 3. Select top genes for labeling
# ==================================
top_genes <- astro %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    slice_head(n = 10)

# ==================================
# 4. Plot
# ==================================
astro_plt <-ggplot(astro, aes(x = log2FoldChange, y = neg_log10_padj)) +
    
    # Points
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.5) +
    
    # Threshold lines
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "black",
               size = 0.4) +
    
    # Labels
    geom_text_repel(
        data = top_genes,
        aes(label = Gene.Name),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        max.overlaps = Inf
    ) +
    
    # Manual clean colors
    scale_color_manual(values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not Significant" = "grey70"
    )) +
    
    # Symmetric X-axis
    coord_cartesian(xlim=c(-10,10)) +
    # Labels
    labs(
        title = "Volcano Plot of Differential Expression : GSE196219",
        subtitle = "Significant genes defined as |log2FC| > 1 and FDR < 0.05",
        x = expression(Log[2]~Fold~Change),
        y = expression(-Log[10]~Adjusted~P~Value),
        color = ""
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "top"
    )


ggsave("GSE196219_astro.jpeg",
       plot = astro_plt,
       width = 10,
       height = 6,
       dpi = 300)