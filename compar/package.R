# 1. Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install GOSemSim (this is the stable, supported version)
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("stringr")
install.packages("forcats")
