##### create Seurat object after downstream analysis in SCANPY #####
# author: Juliska E Boer
# date: 03 Nov 2020

#load packages
setwd("E:/")
library(Seurat)

#read raw- and batch corrected expression matrices and metadata
raw_expr <- read.table("data/input/May2020_EmbryoHb_Seurat_rawexpr.csv", header=TRUE, sep=",", dec=".", row.names = 1)
batch_expr <- read.table("data/input/May2020_EmbryoHb_Seurat_batchexpr.csv", header=TRUE, sep=",", dec=".", row.names = 1)
meta <- read.table("data/input/May2020_EmbryoHb_Seurat_meta.csv", header=TRUE, sep=",", dec=".", row.names = 1)

#create Seurat object using raw counts and metadata
final.embryo <- CreateSeuratObject(counts = raw_expr, assay="RNA", meta.data=meta)
#include batch corrected matrix as scaled data matrix
final.embryo@assays$RNA@scale.data <- as.matrix(batch_expr)
#set Louvain clustering as identity
Idents(final.embryo) <- final.embryo[["louvain"]]

#save Seurat object
save(final.embryo, file="data/output/merge_adult/Embryo_Scanpy_Seurat_obj.RData")

#save average expression per cluster for MAGMA analysis
embryo_avg <- AverageExpression(final.embryo, assays = "RNA")[[1]]
write.table(embryo_avg, file="data/output/GWAS/Jun2020_GWAS_embryo.csv")