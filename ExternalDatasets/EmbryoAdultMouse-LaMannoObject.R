##### Seurat object for developmental/adult (dopamine) mouse dataset from La Manno/Linnarsson 2016 #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)

#adult file
adultDA_all <- read.csv("data/input/GSE76381_MouseAdultDAMoleculeCounts.cef.txt", stringsAsFactors=FALSE)

#create adult metadata
adultDA_meta <- adultDA_all[2:3,2:dim(adultDA_all)[2]]
adultDA_meta <- data.frame(t(rbind(adultDA_meta, c("stage", rep("adult", dim(adultDA_all)[2])))), row.names=1)
colnames(adultDA_meta) <- c("celltype", "stage")
adultDA_meta <- adultDA_meta[-1,]

#create adult count matrix
adultDA_count <- adultDA_all[5:dim(adultDA_all)[1],-2]
rownames(adultDA_count) <- adultDA_count[,1]
adultDA_count <- adultDA_count[,-1]
colnames(adultDA_count) <- row.names(adultDA_meta)

#create Seurat object
mouse_ad <- CreateSeuratObject(adultDA_count, assay = "RNA", meta.data = adultDA_meta, project="AdultMouseDA")
mouse_ad@active.ident <- mouse_ad@meta.data$celltype

#embryonic file
embryo_all <- read.csv("data/input/GSE76381_MouseEmbryoMoleculeCounts.cef.txt", stringsAsFactors=FALSE, sep="\t")

#create embryo metadata
embryo_meta <- embryo_all[2:4,2:dim(embryo_all)[2]]
embryo_meta <- data.frame(t(embryo_meta), row.names=1)
colnames(embryo_meta) <- c("celltype", "stage")
embryo_meta <- embryo_meta[-1,]

#create embryo count matrix
embryo_count <- embryo_all[7:dim(embryo_all)[1],-2]
rownames(embryo_count) <- embryo_count[,1]
embryo_count <- embryo_count[,-1]
colnames(embryo_count) <- row.names(embryo_meta)

#create Seurat object
mouse_em <- CreateSeuratObject(embryo_count, assay = "RNA", meta.data = embryo_meta, project="EmbryoMouse")
mouse_em@active.ident <- mouse_em@meta.data$celltype

#merge objects
lamanno <- merge(mouse_ad, mouse_em, add.cell.ids = NULL)

lamanno <- FindVariableFeatures(lamanno)
lamanno <- NormalizeData(lamanno)
lamanno <- ScaleData(lamanno)

#remove cells from unknown ident
Unk <- WhichCells(lamanno, idents = "mUnk")
lamanno <- subset(lamanno, cells = Unk, invert = TRUE)

#export the average expression per cluster for MAGMA analysis
lamanno_avg <- AverageExpression(lamanno, assays = "RNA")[[1]]
write.table(lamanno_avg, file="data/output/MAGMA/Nov2020_GWAS_lamanno2016.csv")

#save the Seurat object
save(lamanno, file="data/output/ExternalDatasets/DevAdMouse-LaManno2016_obj.RData")
