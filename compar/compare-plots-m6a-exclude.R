n_2 <- readRDS("cluster_summary_labeled.RDS")

library(dplyr)
cluster_rank_2 <- n_2 %>%
    group_by(Description,ID) %>%
    summarise(n_datasets = n_distinct(dataset),
              mean_abs_NES = mean(abs(mean_NES)
              )
    ) %>%
    arrange(desc(n_datasets),desc(mean_abs_NES))


library(dplyr)
# all shared pathways 
top_clusters_2 <- cluster_rank_2 %>%
    filter(n_datasets >= 3)
top_neuron <- top_clusters_2$ID
#looking at only top 15
top_clusters_4 <- cluster_rank_2 %>%
    filter(n_datasets >= 3) %>%
    slice_head(n = 15)
#-------------------------------------------------
# now preparing astro part :
astro_gse_sig <- readRDS("~/ALS/down_stream/astro_gse_sig.RDS")

astro_subset <- astro_gse_sig %>%
    filter(ID %in% top_neuron) %>%
    mutate(in_B = ifelse(p.adjust < 0.05, TRUE, FALSE)) %>%
    dplyr::select(ID,in_B)

astro_full <- top_clusters_2 %>%
    left_join(astro_subset, by = "ID") %>%
    mutate(in_B = ifelse(is.na(in_B),FALSE,in_B))
#----------plotting :

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
#--------------------
plot_data <- astro_full %>%
    mutate(
        Description = paste0(
            toupper(substr(Description, 1, 1)),
            substr(Description, 2, nchar(Description))
        ),
        Description = str_wrap(Description, width = 35),
        Description = fct_reorder(Description, mean_abs_NES),
        n_dataset = as.factor(n_datasets)
    )
z <- ggplot(plot_data, aes(x = mean_abs_NES, y = Description)) +
    
    geom_segment(
        aes(x = 0, xend = mean_abs_NES, yend = Description),
        color = "gray88",
        linewidth = 0.5
    ) +
    
    geom_point(
        aes(color = n_datasets, shape = in_B),
        size = 4,
        stroke = 1.2
    ) +
    
    scale_color_manual(
        values = c(
            "3" = "#56B4E9",
            "4" = "#D55E44"
        )
    ) +
    
    scale_shape_manual(
        values = c(
            "FALSE" = 21,
            "TRUE" = 16
        ),
        labels = c(
            "FALSE" = "No",
            "TRUE" = "Yes"
        )
    ) +
    
    labs(
        x = "Mean Absolute Normalized Enrichment Score (|NES|)",
        y = NULL,
        shape = "Also enriched in astrocytes",
        color = "Motor-neuron datasets sharing pathway"
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        
        # Bigger, clearer pathway names
        axis.text.y = element_text(
            color = "black",
            size = 10.5,
            lineheight = 0.9
        ),
        
        axis.text.x = element_text(
            color = "black",
            size = 11
        ),
        
        axis.title.x = element_text(
            size = 13,
            margin = margin(t = 10),
            face = "bold"
        ),
        
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 10.5),
        
        legend.position = "right",
        
        plot.margin = margin(
            t = 10,
            r = 10,
            b = 10,
            l = 10
        )
    )

ggsave(
    "pathway_asrovsmotor.jpeg",
    plot = z,
    width = 11,
    height = 9.5,
    dpi = 300
)
