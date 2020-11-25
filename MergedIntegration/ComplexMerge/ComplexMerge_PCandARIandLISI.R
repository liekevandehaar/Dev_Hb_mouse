##### batch effect correction for complex merge using Seurat v3 #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)
library(mclust)
library(lisi)
library(dplyr)
library(plyr)

#read the data
#Seurat object of developmental Hb dataset for PC validation
load("data/output/DevelopmentalHb/Embryo_Scanpy_Seurat_obj.RData")
hb_pc <- WhichCells(final.embryo, idents=c("12", "13"))
#Seurat object after Seurat v3 batch effect correction and downstream clustering in Seurat
load("data/output/MergedIntegration/ComplexMerge/Merged_DevHB_H_W.RData")
#meta data and UMAP coordinates after BBKNN batch effect correction and downstream clustering in SCANPY
meta_bbknn <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_BBKNN_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_bbknn <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_BBKNN_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and t-SNE coordinates after ComBat batch effect correction and downstream clustering in SCANPY
meta_combat <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_ComBat_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
tsne_combat <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_ComBat_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and t-SNE coordinates before batch effect correction. downstream clustering is performed in SCANPY
meta_raw <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_raw_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
tsne_raw <- read.csv("data/output/MergedIntegration/ComplexMerge/BatchTest_raw_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
colnames(tsne_raw) <- c("tSNE_1", "tSNE_2")
colnames(umap_bbknn) <- c("UMAP_1", "UMAP_2")
colnames(tsne_combat) <- c("tSNE_1", "tSNE_2")
colnames(tsne_scgen) <- c("tSNE_1", "tSNE_2")

#PC, ARI- and LISI batch and cell type for Seurat v3
pc_merged = data.frame(unclass(table(FetchData(merged, vars="seurat_clusters", cells = hb_pc))))
colnames(pc_merged) = c("n/o_pc")
#highest percentage PC clustered together
(max(pc_merged)/length(hb_pc))*100

cluster = sapply(pc_merged, function(x) head(row.names(pc_merged)[order(x, decreasing = TRUE)], 1))
highest_cluster = WhichCells(merged, idents=cluster)
#percentage of PC in that cluster
(length(intersect(highest_cluster, hb_pc))/length(highest_cluster))*100

bARI_data = select(merged@meta.data, orig.ident, seurat_clusters)
bARI_data$orig.ident <- mapvalues(ARI_data$orig.ident, from = c("Habenula", "10X_LHb", "hab"), to = c(0, 1, 2))
bARI_res = adjustedRandIndex(as.numeric(ARI_data$orig.ident), as.numeric(ARI_data$seurat_clusters))
#batch ARI
(bARI_res-(-1))/(1-(-1))

cARI_data <- select(merged@meta.data, orig.celltype, seurat_clusters)
cARI_data$orig.celltype <- mapvalues(cARI_data$orig.celltype, from = c(levels(cARI_data$orig.celltype)), to = seq(0,34))
cARI_res <- adjustedRandIndex(as.numeric(cARI_data$orig.celltype), as.numeric(cARI_data$seurat_clusters))
#cell type ARI
(cARI_res-(-1))/(1-(-1)) 

LISI_data <- Embeddings(merged[["tsne"]])
bLISI_res <- compute_lisi(LISI_data, object[["orig.ident"]], c("orig.ident"))
#batch LISI
median(bLISI_res$orig.ident)

cLISI_res <- compute_lisi(LISI_data, object[["orig.celltype"]], c("orig.celltype"))
#cell type LISI
median(cLISI_res$orig.celltype)

#ARI- and LISI batch and cell type for BBKNN
bARI_data <- select(meta_bbknn, orig.ident, louvain)
bARI_data$orig.ident <- mapvalues(bARI_data$orig.ident, from = c("Developmental Hb", "Hashikawa", "Wallace"), to = c(0, 1, 2))
bARI_res <- adjustedRandIndex(as.numeric(bARI_data$orig.ident), as.numeric(bARI_data$louvain))
#batch ARI
(bARI_res-(-1))/(1-(-1))

cARI_data <- select(meta_bbknn, orig.celltype, louvain)
cARI_data$orig.celltype <- mapvalues(cARI_data$orig.celltype, from = c(levels(cARI_data$orig.celltype)), to = seq(0,34))
cARI_res <- adjustedRandIndex(as.numeric(cARI_data$orig.celltype), as.numeric(cARI_data$louvain))
#cell type ARI
(cARI_res-(-1))/(1-(-1))

#lisi
bLISI_res <- compute_lisi(umap_bbknn, meta_bbknn, c("orig.ident"))
#batch LISI
median(bLISI_res$orig.ident)

cLISI_res <- compute_lisi(umap_bbknn, meta_bbknn, c("orig.celltype"))
#cell type LISI
median(cLISI_res$orig.celltype)

#ARI- and LISI batch and cell type for ComBat
bARI_data <- select(meta_combat, orig.ident, louvain)
bARI_data$orig.ident <- mapvalues(bARI_data$orig.ident, from = c("Developmental Hb", "Hashikawa", "Wallace"), to = c(0, 1, 2))
bARI_res <- adjustedRandIndex(as.numeric(bARI_data$orig.ident), as.numeric(bARI_data$louvain))
#batch ARI
(bARI_res-(-1))/(1-(-1))

cARI_data <- select(meta_combat, orig.celltype, louvain)
cARI_data$orig.celltype <- mapvalues(cARI_data$orig.celltype, from = c(levels(cARI_data$orig.celltype)), to = seq(0,34))
cARI_res <- adjustedRandIndex(as.numeric(cARI_data$orig.celltype), as.numeric(cARI_data$louvain))
#cell type ARI
(cARI_res-(-1))/(1-(-1))

#lisi
bLISI_res <- compute_lisi(tsne_combat, meta_combat, c("orig.ident"))
#batch LISI
median(bLISI_res$orig.ident)

cLISI_res <- compute_lisi(tsne_combat, meta_combat, c("orig.celltype"))
#cell type LISI
median(cLISI_res$orig.celltype)

#ARI- and LISI batch and cell type for raw
bARI_data <- select(meta_raw, orig.ident, louvain)
bARI_data$orig.ident <- mapvalues(bARI_data$orig.ident, from = c("Developmental Hb", "Hashikawa", "Wallace"), to = c(0, 1, 2))
bARI_res <- adjustedRandIndex(as.numeric(bARI_data$orig.ident), as.numeric(bARI_data$louvain))
#batch ARI
(bARI_res-(-1))/(1-(-1))

cARI_data <- select(meta_raw, orig.celltype, louvain)
cARI_data$orig.celltype <- mapvalues(cARI_data$orig.celltype, from = c(levels(cARI_data$orig.celltype)), to = seq(0,34))
cARI_res <- adjustedRandIndex(as.numeric(cARI_data$orig.celltype), as.numeric(cARI_data$louvain))
#cell type ARI
(cARI_res-(-1))/(1-(-1))

#lisi
bLISI_res <- compute_lisi(tsne_raw, meta_raw, c("orig.ident"))
#batch LISI
median(bLISI_res$orig.ident)

cLISI_res <- compute_lisi(tsne_raw, meta_raw, c("orig.celltype"))
#cell type LISI
median(cLISI_res$orig.celltype)
