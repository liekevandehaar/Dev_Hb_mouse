##### function for calculating the rank-order correlations between two performed MAGMA protocols #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(Seurat)
library(MAST)

#load downloaded Seurat object
#neuronal celltypes
load("data/input/facs_Brain_Microglia_seurat_tiss.Robj")
tiss_mglia <- UpdateSeuratObject(tiss)

#non-neuronal cell types
load("data/input/facs_Brain_Nonmicroglia_seurat_tiss.Robj")
tiss_non <- UpdateSeuratObject(tiss)
rm(tiss)

#merge the neuronal and non-neuronal objects
brain_tiss <- merge(tiss_mglia, tiss_non)
Idents(brain_tiss) <- brain_tiss[["annotation"]]

#filter out the cells from the "Unknown" cluster
CellstoRemove <- WhichCells(brain_tiss, idents = "unknown")
brain_tiss_filtered <- brain_tiss[,!colnames(brain_tiss) %in% CellstoRemove]

#save Seurat object holding all neuronal and non-neuronal cell types from the Tabula Muris, et al (2018) publication
save(brain_tiss_filtered, file="data/output/ExternalDatasets/TabulaMuris_Seurat_object.RData")

#save average expression for MAGMA analysis
tabulam_avg <- AverageExpression(brain_tiss_filtered, assays = "RNA")[[1]]
write.table(tabulam_avg, file="data/output/MAGMA/TabulaMuris_avgexpr.csv")