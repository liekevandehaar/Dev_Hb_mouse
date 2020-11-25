##### Seurat object for developmental hypothalamus dataset #####
# author: Juliska E Boer
# date: 09 Nov 2020

#load packages
setwd("E:/")
library(Seurat)

#developmental mouse hypothalamus dataset of Kim, et al (2020)

#read in published meta data and counts
meta <- read.csv("data/GSE132355_E10-P45_umap_data.csv", header=T, row.names=1)
data <- readRDS("data/GSE132355_E10-P45_log2cpm.rds")

#create Seurat object
hypo <- CreateSeuratObject(data, assay = "RNA", meta.data = meta, names.field = 1, names.delim = "_", project="DevHypothalamus")
hypo@active.ident <- hypo@meta.data$Cluster

#save full object
save(hypo, file="data/output/GWAS/DevHypothalamus_obj.RData")

#show age distribution per cluster for sampling
age_clusters <- data.frame(unclass(table(hypo@meta.data$Cluster, hypo@meta.data$Age)))

#select 14 clusters that resemble the age distribution of the dev Hb dataset
Idents(hypo) <- hypo@meta.data$Cluster
hypo_14 <- subset(hypo, idents = c("Endothelial cells", "Hypothalamic neurons", "Immune cells", "Non-neuronal cells", "Ventral hypothalamic NPC",
                                   "LH", "POA", "PMN", "VMH", "Astrocytes", "Ependymal cells", "Microglia", "Oligodendrocytes", "Tanycytes"))

#sample 5% of cells in each selected cluster
ident <- unique(hypo_14@active.ident)
subcells <- c()
for(item in ident){
  set.seed(42)
  cells <- sample(WhichCells(hypo_14, idents = item), size=round(length(WhichCells(hypo_14, idents=item))*0.05), replace=F)
  subcells <- c(subcells, cells)
}
subhypo <- subset(hypo_14, cells=subcells)

#write average expression per cluster to file
hypo_avg <- AverageExpression(subhypo, assays = "RNA")[[1]]
write.table(hypo_avg, file="data/output/GWAS/Nov2020_GWAS_devhypo-samp.csv")

#save the sampled object
save(subhypo, file="data/output/GWAS/DevHypothalamus-samp_obj.RData")