##### batch effect correction for complex merge using Seurat v3 #####
# author: Juliska E Boer
# date: 03 Nov 2020

#load packages
setwd("E:/")
library(Seurat)
library(ggplot2)
library(dplyr)

#load data
load("data/output/merge_adult/Embryo_Scanpy_Seurat_obj.RData") #developmental Hb
stuber <- readRDS("data/input/Habenula_neuron_Seurat_object.rds") #Hashikawa, et al 2020
load("data/output/merge_adult/Wallace_Seurat_obj.RData") #Wallace, et al 2020

#remove mitochondrial genes from both adult mouse Hb datsets
st_counts <- GetAssayData(stuber, assay = "RNA")
st_counts <- st_counts[rownames(st_counts)[!grepl("^mt", rownames(st_counts))],]
st_nomito <- subset(stuber, features = rownames(st_counts))

w_counts <- GetAssayData(wallace_hb, assay = "RNA")
w_counts <- w_counts[rownames(w_counts)[!grepl("^mt", rownames(w_counts))],]
w_nomito <- subset(wallace_hb, features = rownames(w_counts))

#perform simple merge and export to SCANPY for batch effect correction using BBKNN, ComBat and scGEN
scanpy <- merge(final.embryo, c(st_nomito, w_nomito))
columns.to.remove <- c("nFeature", "nCount", "n_total_counts", "n_counts_norm", "plate", "stage", "ERCC_genes", "percent_mito", "percent_ribo", "stim",
                       "percent.mito", "nCount_integrated", "nFeature_integrated", "integrated_snn_res.0.8", "louvain", "clusters", "celltype")
for(i in columns.to.remove) {
  scanpy[[i]] <- NULL
}
SeuratDisk::SaveH5Seurat(scanpy, filename = "data/output/merge_adult/BatchTest_Merged_Oct20.h5Seurat")
SeuratDisk::Convert("data/output/merge_adult/BatchTest_Merged_Oct20.h5Seurat", dest="h5ad")

#export text file holding all earlier identified clusters and cell types
idents_devhb <- final.embryo@active.ident
idents_hashi <- st_nomito@active.ident
idents_wallace <- w_nomito@active.ident
all_idents <- c(idents_devhb[[1]], idents_hashi[[1]], as.character(idents_wallace[[1]]))
data.table::fwrite(list(all_idents), file="data/output/merge_adult/BatchTest_identlist_26Oct20.txt")

#select highly variable genes for each object
final.embryo <- FindVariableFeatures(final.embryo, assay="RNA", selection.method="mvp")

st_nomito <- ScaleData(st_nomito, assay="RNA")
st_nomito <- FindVariableFeatures(st_nomito, assay="RNA", selection.method="mvp")

w_nomito <- FindVariableFeatures(w_nomito, assay="RNA", selection.method = "mvp")

#search, filter and score neighbors as anchors for integration
merge.anchors <- FindIntegrationAnchors(object.list = list(final.embryo, st_nomito, w_nomito), dims = 1:25, k.anchor=5, k.filter=40, k.score=20)
#perform integration using identified anchors
merged <- IntegrateData(anchorset =  merge.anchors, dims = 1:30, features.to.integrate = c(rownames(st_nomito), rownames(w_nomito), rownames(final.embryo)))
DefaultAssay(merged) <- "integrated"
#set the earlier identified clusters as metadata
merged[["orig.celltype"]] <- merged@active.ident

#scale data and perform PCA
merged <- ScaleData(merged, verbose = TRUE)
merged <- RunPCA(merged, npcs = 30, verbose = TRUE)

#calculate t-SNE coordinates
merged <- RunTSNE(merged)

#calculate nearest neighborhood graph and perform Louvain clustering
merged <- FindNeighbors(merged)
merged <- FindClusters(merged, resolution = 0.8)

#save the merged object for validation using ARI and LISI
save(merged, file="data/output/merge_adult/Merged_DevHB_H_W.RData")

#export the merged object to SCANPY for plotting
SeuratDisk::SaveH5Seurat(merged, filename = "data/output/merge_adult/Merged4TSNE.h5Seurat")
SeuratDisk::Convert("data/output/merge_adult/Merged4TSNE.h5Seurat", dest="h5ad")

#save average expression of the clusters in the merged object for MAGMA analysis
avg_expr <- AverageExpression(merged, assays="RNA", slot="data")[[1]]
write.table(avg_expr, file="data/output/merge_adult/Oct2020_GWAS_merged.csv")