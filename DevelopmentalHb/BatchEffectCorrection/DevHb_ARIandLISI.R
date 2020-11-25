##### perform validation using ARI and LISI for sampling time point batch effect correction #####
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
load("data/output/DevelopmentalHb/Embryo_Scanpy_Seurat_obj.RData")
#meta data and UMAP coordinates after BBKNN batch effect correction and downstream clustering in SCANPY
meta_bbknn <- read.csv("data/output/DevelopmentalHb/BatchTest_BBKNN_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
umap_bbknn <- read.csv("data/output/DevelopmentalHb/BatchTest_BBKNN_umap.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and t-SNE coordinates after ComBat batch effect correction and downstream clustering in SCANPY
meta_combat <- read.csv("data/output/DevelopmentalHb/BatchTest_ComBat_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
tsne_combat <- read.csv("data/output/DevelopmentalHb/BatchTest_ComBat_tsne.csv", header=TRUE, row.names = 1, sep=",", dec=".")
#meta data and t-SNE coordinates before batch effect correction. downstream clustering is performed in SCANPY
meta_raw <- read.csv("data/output/DevelopmentalHb/BatchTest_raw_metadata.csv", header=TRUE, row.names = 1, sep=",", dec=".")
tsne_raw <- read.csv("data/output/DevelopmentalHb/BatchTest_raw_tsne.csv", header=TRUE, row.names = 1, sep=",", dec=".")
colnames(tsne_raw) <- c("tSNE_1", "tSNE_2")
colnames(umap_bbknn) <- c("UMAP_1", "UMAP_2")
colnames(tsne_combat) <- c("tSNE_1", "tSNE_2")

#ARI and LISI for Seurat v3
ARI_data = select(final.embryo@meta.data, "stage", "louvain")
ARI_data$stage <- mapvalues(ARI_data$stage, from = c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult"), to = c(0, 1, 2, 3, 4, 5, 6, 7))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$stage), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(Embeddings(final.embryo[["tsne"]]), final.embryo@meta.data, c("stage"))
median(LISI_res$stage)

#ARI and LISI for BBKNN
ARI_data = select(meta_bbknn, "stage", "louvain")
ARI_data$stage <- mapvalues(ARI_data$stage, from = c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult"), to = c(0, 1, 2, 3, 4, 5, 6, 7))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$stage), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(umap_bbknn, meta_bbknn, c("stage"))
median(LISI_res$stage)

#ARI and LISI for ComBat
ARI_data = select(meta_combat, "stage", "louvain")
ARI_data$stage <- mapvalues(ARI_data$stage, from = c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult"), to = c(0, 1, 2, 3, 4, 5, 6, 7))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$stage), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(tsne_combat, meta_combat, c("stage"))
median(LISI_res$stage)

#ARI and LISI for raw
ARI_data = select(meta_raw, "stage", "louvain")
ARI_data$stage <- mapvalues(ARI_data$stage, from = c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult"), to = c(0, 1, 2, 3, 4, 5, 6, 7))
ARI_res = adjustedRandIndex(as.numeric(ARI_data$stage), as.numeric(ARI_data$louvain))
(ARI_res-(-1))/(1-(-1))

LISI_res <- compute_lisi(tsne_raw, meta_raw, c("stage"))
median(LISI_res$stage)
