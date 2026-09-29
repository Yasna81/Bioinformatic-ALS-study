#compartive plots
neuron <- readRDS("~/ALS/top_clusters.RDS")
n_2 <- readRDS("~/ALS/pathway_summary.RDS") #cluster_summary_labeled
library(dplyr)
cluster_rank <- n_2 %>%
    group_by(Description,ID) %>%
    summarise(n_datasets = n_distinct(dataset),
              mean_abs_NES = mean(abs(mean_NES)
                                  )
    ) %>%
    arrange(desc(n_datasets),desc(mean_abs_NES))
# 15 top clusters that apeared in 4 or more than 4 datasets (reproducible)
top_clusters <- cluster_rank %>%
    filter(n_datasets >= 4) %>%
    slice_head( n = 15)
top_neuron <- top_clusters$ID
#----------------------------------------
astro_subset <- astro_gse_sig %>%
    filter(ID %in% top_neuron) %>%
    mutate(in_B = ifelse(p.adjust < 0.05, TRUE, FALSE)) %>%
    select(ID,in_B)

astro_full <- top_clusters %>%
    left_join(astro_subset, by = "ID") %>%
    mutate(in_B = ifelse(is.na(in_B),FALSE,in_B))
#plot 
library(ggplot2)
library(dplyr)
library(ggrepel)
#make sure true and false are factors

astro_full <- astro_full %>%
    mutate(
        n_datasets = factor(n_datasets)
    )
# Select top shared pathways to annotate (e.g., top 5 by mean_abs_nes)
top_shared <- astro_full %>%
    filter(in_B == "TRUE") %>%
    arrange(desc(mean_abs_NES))


#---------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr) # For wrapping long labels
library(forcats) # For perfect sorting

# 1. Pre-process the data for perfect sorting and readability
plot_data <- astro_full %>%
    mutate(
        # Wraps long pathway names into multiple lines (30 characters max)
        Description = str_wrap(Description, width = 30),
        # Strictly reorders the Y-axis by the NES value
        Description = fct_reorder(Description, mean_abs_NES),
        # Ensure n_dataset is treated as a discrete factor for coloring
        n_dataset = as.factor(n_datasets)
    )

# 2. Generate the plot
z <- ggplot(plot_data, aes(x = mean_abs_NES, y = Description)) +
    # Lollipop lines (cleaned up)
    geom_segment(
        aes(x = 0, xend = mean_abs_NES, yend = Description),
        color = "gray85",
        linewidth = 0.6
    ) +
    # Points with improved aesthetic mapping
    geom_point(
        aes(color = n_datasets, shape = in_B),
        size = 4,
        stroke = 1.2
    ) +
    # FIX: Use _d for discrete factors to avoid the error you saw
    scale_color_manual(
        values = c("4" = "#56B4E9","5"="#D55E44")
    ) +
    # Clean up shapes: 21 is a circle that allows for distinct TRUE/FALSE visual
    scale_shape_manual(
        values = c("FALSE" = 21, "TRUE" = 16),
        labels = c("FALSE" = "No", "TRUE" = "Yes")
    ) +
    labs(
        x = "Mean Absolute Normalized Enrichment Score (|NES|)",
        y = NULL,
        shape = "Enriched in Astrocytes",
        color = "N-Datasets Sharing Pathway"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), # Cleans the background for the lollipops
        axis.text.y = element_text(color = "black", size = 9, lineheight = 0.8),
        axis.title.x = element_text(margin = margin(t = 10), face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5)
    )


ggsave("pathway_asrovsmotor.jpeg",
       plot = z,
       width = 10,
       height = 6,
       dpi = 300)
#-----------------------------------------------------------------------

