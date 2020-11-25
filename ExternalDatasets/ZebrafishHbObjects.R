##### Seurat object for zebrafish Hb datasets #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)

#read in larval zebrafish Hb object
load("data/input/Larva10x_0805.RObj")

#object is Seurat v1, create Seurat v3 object
pandey <- CreateSeuratObject(counts = Larva10xF@data, meta.data = Larva10xF@data.info, project = "LarvalZebrafish", assay="RNA")
pandey@assays$RNA@scale.data <- Larva10xF@scale.data
Idents(pandey) <- pandey[["res.0.5"]]

#Normalize data and perform PCA
pandey <- NormalizeData(pandey, assay="RNA")
pandey <- RunPCA(pandey)

#select highly variable genes and scale data
pandey <- FindVariableFeatures(pandey)
pandey <- ScaleData(pandey, assay="RNA")

#calculate t-SNE coordinates
pandey <- RunTSNE(pandey)

#change cluster names (numbers) to actual names from publication
pandey.id <- RenameIdents(object = pandey, "1"="La_Hb01", "2"="La_Hb02", "3"="La_Hb03", "4"="La_Hb04","5"="La_Hb05", "6"="La_Hb06", 
                          "7"="La_Hb07", "8"="La_Hb08", "9"="La_Hb09", "10"="La_Hb10", "11"="La_Hb11", "12"="La_Hb12", "13"="La_Hb13", 
                          "14"="La_Hb14", "15"="La_Hb15", "16"="Olf")
#save larval zebrafish Hb Seurat object
save(pandey.id, file="data/output/ExternalDatasets/Pandey_Seurat_obj.RData")

#read in adult zebrafish Hb object
load("data/input/Adult_Hab_NewSeurat.RObj")

#object in Seurat v2, update to Seurat v3
pandey.ad.inter <- UpdateSeuratObject(Adult_CB)

#remove subset of unidentified cells
CellsX <- WhichCells(pandey.ad.inter, idents = "x")
pandey.ad <- subset(pandey.ad.inter, cells = CellsX, invert = TRUE)

#save adult zebrafish Hb Seurat object
save(pandey.ad, file="data/output/ExternalDatasets/PandeyAdult_Seurat_obj.RData")