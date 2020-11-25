##### batch effect correction in simple merge using Seurat v3 #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)
library(SeuratDisk)

#load filtered Diencephalon data
Convert("data/output/MergedIntegration/SimpleMerge/Jun2020_merge_diencephalon_filtered_di.h5ad", dest = "h5seurat", overwrite = TRUE)
dienc <- LoadH5Seurat("data/output/MergedIntegration/SimpleMerge/Jun2020_merge_diencephalon_filtered_di.h5seurat")

#load Developmental Hb data
load("data/output/DevelopmentalHb/Embryo_Scanpy_Seurat_obj.RData")

#perform normalizdation on diencephalon data
dienc <- NormalizeData(dienc)
#and determine highly variable genes
dienc <- FindVariableFeatures(dienc, assay="RNA", selection.method="mvp")

#filter the Developmental Hb dataset to only include cells from E12
e12 <- subset(final.embryo, subset = stage=="E12")
#and determine highly variable genes
e12 <- FindVariableFeatures(e12, assay="RNA", selection.method = "mvp")

#find, filter and score neighbors for locating anchors (default parameters)
merge.anchors <- FindIntegrationAnchors(object.list = list(e12, dienc))
#integrate two datasets together
merged <- IntegrateData(anchorset =  merge.anchors, dims = 1:30, features.to.integrate = c(rownames(dienc), rownames(e12)))
DefaultAssay(merged) <- "integrated"
merged[["orig.ident"]] <- replace_na(merged@meta.data$orig.ident, "diencephalon")

#scale data and perform PCA
merged <- ScaleData(merged, verbose = TRUE)
merged <- RunPCA(merged, npcs = 30, verbose = TRUE)

#calculate UMAP coordinates
merged <- RunUMAP(merged, dims=1:5)

#calculate nearest neighborhood graph and perform Louvain clustering
merged <- FindNeighbors(merged)
merged <- FindClusters(merged, resolution = 0.8)

#save the merged object for validation using ARI and LISI
save(merged, file="data/output/MergedIntegration/SimpleMerge/Diencephalon_Merged_Seuratobject.RData")

#export the merged object to SCANPY for plotting
DefaultAssay(merged) <- "RNA" #we need the uncorrected data matrix for gene plotting
SaveH5Seurat(merged, filename = "data/output/MergedIntegration/SimpleMerge/Diencephalon_MergedSeurat.h5Seurat")
Convert("data/output/MergedIntegration/SimpleMerge/Diencephalon_MergedSeurat.h5Seurat", dest="h5ad")