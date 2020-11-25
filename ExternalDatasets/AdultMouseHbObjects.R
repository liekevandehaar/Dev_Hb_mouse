##### Seurat object for mouse adult Hb datasets #####
# author: Juliska E Boer
# date: 03 Nov 2020

#load packages
setwd("E:/")
library(Seurat)

#adult mouse Hb dataset of wallace, et al (2020)

#read in published Seurat object
wallace <- readRDS("data/input/hab_batch1.rds")
wallace <- UpdateSeuratObject(wallace)

#read in supplied metadata file holding the original heterogeneous annotations for each cell (e.g. microglia, neurons, pericytes, etc)
load("data/input/W_Habenula_Seurat_meta.RData")
wallace[["orig_clusters"]] <- c(meta$CellClassNames_filtered)
#read in supplied metadata file holding the lateral Hb annotations
load("data/input/lhb_Seurat_meta.RData")
#read in supplied metadata file holding the medial Hb annotations
load("data/input/mhb_Seurat_meta.RData")

#filter Seurat object on medial and lateral Hb cells
meta_lhb <- transform(meta_lhb, tree.ident = sprintf('LHb%d', tree.ident))
meta_mhb <- transform(meta_mhb, tree.ident = sprintf('MHb%d', tree.ident))
wallace_hb <- subset(wallace, cells = c(rownames(meta_mhb),rownames(meta_lhb)))

#set identities of said cells to annotated clusters (format = numbers)
wallace_hb[["clusters"]] <- c(meta_mhb$tree.ident,meta_lhb$tree.ident)
Idents(wallace_hb) <- wallace_hb[["clusters"]]

#rename clusters (identities) to annotated clusters (format = names)
wallace_hb <- RenameIdents(object = wallace_hb, "LHb1"="Oval/Medial(LHb1)", "LHb2"="Marginal(LHb2)", "LHb3"="Lateral(LHb3)", "LHb4"="HBX(LHb4)",
                           "MHb1"="Ventral2/3(MHb1)", "MHb2"="Ventrolateral(MHb2)", "MHb3"="Lateral(MHb3)", "MHb4"="Dorsal(MHb4)", "MHb5"="Superior(MHb5)")
#save Seurat object holding the adult mouse Hb of the Wallace, et al (2020) publication.
save(wallace_hb, file="data/output/merge_adult/Wallace_Seurat_obj.RData")

#save average expression per cluster for MAGMA analysis
wallace_avg <- AverageExpression(wallace_hb, assays = "RNA")[[1]]
write.table(wallace_avg, file="data/output/GWAS/Jun2020_GWAS_wallace.csv")

#adult mouse Hb dataset of Hashikawa, et al (2020)
stuber <- readRDS("data/input/Habenula_neuron_Seurat_object.rds") #only neuronal (MHB and LHB) cell types
stuber_all <- readRDS("data/input/Habenula_Seurat_object.rds") #all cell types

#save average expression per cluster for MAGMA analysis
stuber_avg <- AverageExpression(stuber, assays = "RNA")[[1]]
write.table(stuber_avg, file="data/output/GWAS/Jun2020_GWAS_stuber.csv")
