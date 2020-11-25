##### perform validation using ARI and LISI for batch effect correction in simple merge #####
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
#Seurat object after Seurat v3 batch effect correction and downstream clustering in SCANPY
load("data/output/MergedIntegration/SimpleMerge/Diencephalon_Merged_Seuratobject.RData")
#meta data and UMAP coordinates after BBKNN batch effect correction and downstream clustering in SCANPY
meta_bbknn <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_BBKNN_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_bbknn <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_BBKNN_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and UMAP coordinates after ComBat batch effect correction and downstream clustering in SCANPY
meta_combat <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_ComBat_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_combat <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_ComBat_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and UMAP coordinates after scGEN batch effect correction and downstream clustering in SCANPY
meta_scgen <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_scgen_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_scgen <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_scgen_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and UMAP coordinates before batch effect correction. downstream clustering is performed in SCANPY
meta_raw <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_raw_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_raw <- read.csv("data/output/MergedIntegration/SimpleMerge/BatchTest_raw_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
colnames(umap_raw) <- c("UMAP_1", "UMAP_2")
colnames(umap_bbknn) <- c("UMAP_1", "UMAP_2")
colnames(umap_combat) <- c("UMAP_1", "UMAP_2")
colnames(umap_scgen) <- c("UMAP_1", "UMAP_2")

#ARI and LISI for Seurat v3
ARI_data = select(merged@meta.data, orig.ident, seurat_clusters)
ARI_data$orig.ident <- mapvalues(ARI_data$orig.ident, from = c("Habenula", "diencephalon"), to = c(0, 1))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$orig.ident), as.numeric(ARI_data$seurat_clusters))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(Embeddings(merged[["umap"]]), merged@meta.data, c("orig.ident"))
median(LISI_res$orig.ident)

#ARI and LISI for BBKNN
ARI_data = select(meta_bbknn, batch, louvain)
ARI_data$batch <- mapvalues(ARI_data$batch, from = c("Developmental Hb", "Diencephalon"), to = c(0, 1))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$batch), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(umap_bbknn, meta_bbknn, c("batch"))
median(LISI_res$batch)

#ARI and LISI for ComBat
ARI_data = select(meta_combat, batch, louvain)
ARI_data$batch <- mapvalues(ARI_data$batch, from = c("Developmental Hb", "Diencephalon"), to = c(0, 1))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$batch), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(umap_combat, meta_combat, c("batch"))
median(LISI_res$batch)

#ARI and LISI for scGEN
ARI_data = select(meta_scgen, batch, louvain)
ARI_data$batch <- mapvalues(ARI_data$batch, from = c("Developmental Hb", "Diencephalon"), to = c(0, 1))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$batch), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(umap_scgen, meta_scgen, c("batch"))
median(LISI_res$batch)

#ARI and LISI for raw
ARI_data = select(meta_raw, batch, louvain)
ARI_data$batch <- mapvalues(ARI_data$batch, from = c("Developmental Hb", "Diencephalon"), to = c(0, 1))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$batch), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(umap_raw, meta_raw, c("batch"))
median(LISI_res$batch)
