library(RubsetRankAggreg)
install.packages("RubsetRankAggreg")
library(remotes)
remotes::install_github("raivokolde/RobustRankAggreg") 
library(RobustRankAggreg)

#data sets are : 
nrow(res_fus_clean)#18436
nrow(m6a_res_clean)#25358
nrow(res_df_2_clean)#12052
nrow(multiple_unique)#31224
nrow(polokinas_res_clean)#25829
#polokinas_res_clean$gene_id <- polokinas_res_clean$X
#defining up and down regulated genes fr each dataset
#polokinase : GSE272827
library(dplyr)
up_polo <- polokinas_res_clean %>% filter(stat > 0) %>% pull(gene_id)
down_polo <- polokinas_res_clean %>% arrange(stat) %>%filter(stat< 0) %>% pull(gene_id)
#swift GSE229095
up_swift <- res_df_2_clean %>% filter(stat > 0) %>% pull(gene_id)
down_swift <- res_df_2_clean %>% arrange(stat) %>%filter(stat< 0) %>% pull(gene_id)
#M6A GS242766
up_M6A <- m6a_res_clean %>% filter(stat > 0) %>% pull(gene_id)
down_M6a <- m6a_res_clean %>% arrange(stat) %>%filter(stat< 0) %>% pull(gene_id)
#fus
up_FUS <- res_fus_clean %>% filter(stat > 0) %>% pull(gene_id)
down_FUS <- res_fus_clean %>% arrange(stat) %>%filter(stat< 0) %>% pull(gene_id)
#mltiple
up_multiple <- multiple_unique %>% filter(stat > 0) %>% pull(gene_id)
down_multiple <- multiple_unique %>% arrange(stat) %>%filter(stat< 0) %>% pull(gene_id)
# running RRA
up_list <-list(up_polo,up_swift,up_M6A,up_FUS,up_multiple)
down_list <- list(down_polo,down_swift,down_M6a,down_FUS,down_multiple)
rra_up <- aggregateRanks(up_list)
rra_down <- aggregateRanks(down_list)
#
rra_up$padj <- p.adjust(rra_up$Score, method = "BH")
rra_down$padj <- p.adjust(rra_down$Score, method = "BH")
sig_up <- rra_up %>% filter(padj <0.05)
sig_down <- rra_down %>% filter(padj <0.05)
sig_up$direction <- "UP"
sig_down$direction <- "Down"
meta_genes <- bind_rows(sig_up,sig_down)
library(org.Hs.eg.db)
meta_genes$symbol <- mapIds(org.Hs.eg.db,
                            keys = rownames(meta_genes),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
