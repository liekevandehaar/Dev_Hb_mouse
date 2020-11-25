##### batch effect correction in multi time point dataset using Seurat v3 #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)
library(SeuratDisk)

#read in raw expression matrix exported from SCANPY
embryo_raw <- read.table("data/output/DevelopmentalHb/Apr2020_embryoHb_raw_expr_matrix.csv", 
                         header=TRUE, sep=",", dec=".", row.names=1)
#create Seurat object using the raw counts
embryo <- CreateSeuratObject(embryo_raw, project="EmbryoHb", assay="RNA", names.field = 5, names.delim = "_")
#set the sampling time point (stage) as meta data
embryo[["stage"]] <- embryo@meta.data$orig.ident

#normalize and Scale data
embryo <- NormalizeData(embryo, assay="RNA")
embryo <- ScaleData(embryo, assay="RNA")

#netermine highly variable genes (HVG) and plot
embryo <- FindVariableFeatures(embryo, assay="RNA", selection.method="mvp")
VariableFeaturePlot(embryo, assay="RNA", selection.method = "mvp")

#show t-SNE embeddings with Louvain clusters before batch effect correction
embryo_no_corr <- RunPCA(embryo, assay="RNA")
embryo_no_corr <- RunTSNE(embryo_no_corr, assay="RNA")
embryo_no_corr <- FindNeighbors(embryo_no_corr)
embryo_no_corr <- FindClusters(embryo_no_corr, resolution = 0.8)
TSNEPlot(embryo_no_corr, label=TRUE)
#export the uncorrected Seurat object to SCANPY for final figure plotting
SaveH5Seurat(embryo_no_corr, filename = "data/output/DevelopmentalHb/DevHb_RawSeurat.h5Seurat")
Convert("data/output/DevelopmentalHb/DevHb_RawSeurat.h5Seurat", dest="h5ad")

#INTEGRATION
#divide the object into sampling time poin subsets and determine the HVG per object
e11 <- subset(embryo, idents = "E11")
e11 <- FindVariableFeatures(e11, assay="RNA", selection.method="mvp")
e12 <- subset(embryo, idents = "E12")
e12 <- FindVariableFeatures(e12, assay="RNA", selection.method="mvp")
e13 <- subset(embryo, idents = "E13")
e13 <- FindVariableFeatures(e13, assay="RNA", selection.method="mvp")
e15 <- subset(embryo, idents = "E15")
e15 <- FindVariableFeatures(e15, assay="RNA", selection.method="mvp")
e18 <- subset(embryo, idents = "E18")
e18 <- FindVariableFeatures(e18, assay="RNA", selection.method="mvp")
p4 <- subset(embryo, idents = "P4")
p4 <- FindVariableFeatures(p4, assay="RNA", selection.method="mvp")
p7 <- subset(embryo, idents = "P7")
p7 <- FindVariableFeatures(p7, assay="RNA", selection.method="mvp")
adult <- subset(embryo, idents = "adult")
adult <- FindVariableFeatures(adult, assay="RNA", selection.method="mvp")

#search, filter and score neighbors as anchors for integration
embryo.anchors <- FindIntegrationAnchors(object.list = list(e11, e12, e13, e15, e18, p4, p7, adult), dims = 1:20, k.filter=30, k.score=15, k.anchor=4, anchor.features = 2000)
#select all the genes to return in the integrated object
all_genes <- rownames(e13@assays$RNA@data) #the genesets are identical for all time point objects
#perform integration using identified anchors
embryo.combined <- IntegrateData(anchorset = embryo.anchors, dims = 1:20, features.to.integrate = all_genes)
DefaultAssay(embryo.combined) <- "integrated"

#export batch corrected matrix for SCANPY analysis
fwrite(x = as.data.frame(as.matrix(embryo.combined@assays$integrated@data)), file = "data/input/Seurat_inputAdataX.txt")