##### create Seurat object after downstream analysis in SCANPY #####
# author: Juliska E Boer
# edited: L.L. van de Haar
# date: 21 01 2021

#load packages
setwd("/Volumes/HD_LLH/python/STAR_E18_wt/FINAL")
library(Seurat)

#read raw- and batch corrected expression matrices and metadata
raw_expr <- read.table("/Volumes/HD_LLH/Juliska/data/output/embryo_hb/E18_Brn3a-tLacZ_rawexpr.csv", header=TRUE, sep=",", dec=".", row.names = 1)
#batch_expr <- read.table("data/input/May2020_EmbryoHb_Seurat_batchexpr.csv", header=TRUE, sep=",", dec=".", row.names = 1)
meta <- read.table("/Volumes/HD_LLH/Juliska/data/output/embryo_hb/E18_Brn3a-tLacZ_metadata.csv", header=TRUE, sep=",", dec=".", row.names = 1)

raw_expr <- as.data.frame(t(raw_expr))
#create Seurat object using raw counts and metadata
final.E18Brn3a <- CreateSeuratObject(counts = raw_expr, assay="RNA", meta.data=meta) 
#include batch corrected matrix as scaled data matrix
#final.embryo@assays$RNA@scale.data <- as.matrix(batch_expr)
#set Louvain clustering as identity
Idents(final.E18Brn3a) <- final.E18Brn3a[["louvain"]]

#save Seurat object
save(final.E18Brn3a, file="data/output/final.E18Brn3a_Scanpy_Seurat_obj.RData")

#save average expression per cluster for MAGMA analysis
embryo_avg <- AverageExpression(final.embryo, assays = "RNA")[[1]]
write.table(embryo_avg, file="data/output/MAGMA/Jun2020_GWAS_embryo.csv")