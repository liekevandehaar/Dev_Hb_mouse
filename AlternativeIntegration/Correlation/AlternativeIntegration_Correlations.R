##### alternative integration: correlation analysis #####
# author: Juliska E Boer
# date: 03 Nov 2020

#load packages
setwd("E:/")
library(MAST)
library(Seurat)
library(tidyverse)
library(heatmaply)

#load functions for identifying differentially expressed genes (DEGs) using MAST and for calculating the correlations
source("Correlation_CalcDEG.R")
source("Correlation_CalcCorr.R")

#read in the data
load("data/output/merge_adult/Embryo_Scanpy_Seurat_obj.RData") #developmental Hb
stuber <- readRDS("data/input/Habenula_neuron_Seurat_object.rds") #adult mouse Hb: Hashikawa, et al (2020)
load("data/output/merge_adult/Wallace_Seurat_obj.RData") #adult mouse Hb: Wallace, et al (2020)
load("data/output/merge_zebrafish/Pandey_Seurat_obj.RData") #larval zebrafish Hb: Pandey, et al (2018)
load("data/output/merge_zebrafish/PandeyAdult_Seurat_obj.RData") #adult zebrafish Hb: Pandey, et al (2018)

#determine DEGs for each dataset, and save information (this function takes some time)
DEG.final.embryo <- DE_Gene_Union(final.embryo, levels(final.embryo@active.ident), data.info.col = 1:13)
save(DEG.final.embryo, file = "data/output/merge_adult/DEG-final-embryo.RData")

DE_Gene_Stuber <- DE_Gene_Union(stuber, levels(stuber@active.ident), data.info.col = 1:9)
save(DE_Gene_Stuber, file = "data/output/merge_adult/DEG_Stuber.RData")

DEG.wallace <- DE_Gene_Union(wallace_hb, levels(wallace_hb@active.ident), data.info.col = 1:5)
save(DEG.wallace, file = "data/output/merge_adult/DEG-wallace.RData")

DEG.pandey <- DE_Gene_Union(pandey.id, levels(pandey.id@active.ident), data.info.col = 1:7)
save(DEG.pandey, file="data/output/merge_zebrafish/DEG_pandey.RData")

DEG.pandey.ad <- DE_Gene_Union(pandey.ad, levels(pandey.ad@active.ident), data.info.col = 1:6)
save(DEG.pandey.ad, file="data/output/merge_zebrafish/DEG_pandey_adult.RData")

#prepare input for correlation function
#average expression matrices filtered on DEG

#load DEG result (vectors behind the objects are the cluster names, these are used for plotting)
load("data/output/merge_zebrafish/DEG_pandey_adult.RData") #c("Ad_Hb01", "Ad_Hb02A", "Ad_Hb02B", "Ad_Hb04", "Ad_Hb05", "Ad_Hb06", "Ad_Hb08", "Ad_Hb09", "Ad_Hb10", "Ad_Hb11", "Ad_Hb13", "Ad_Hb14", "Ad_VHb01", "Ad_VHb02", "Ad_VHb03", "Ad_VHb04", "Ad_Hb16")
load("data/output/merge_zebrafish/DEG_pandey.RData") #c("La_Hb01", "La_Hb02", "La_Hb03", "La_Hb04", "La_Hb05", "La_Hb06", "La_Hb07", "La_Hb08", "La_Hb09", "La_Hb10", "La_Hb11", "La_Hb12", "La_Hb13", "La_Hb14", "La_Hb15", "Olf")
load("data/output/merge_adult/DEG_Stuber.RData") #c("MHb1", "MHb2", "MHb3", "MHb4", "MHb5", "MHb6", "LHb1", "LHb2", "LHb3", "LHb4", "LHb5", "LHb6")
load("data/output/merge_adult/DEG-final-embryo.RData") #c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
load("data/output/merge_adult/DEG-wallace.RData") #c("Ventral 2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")

#set the datasets and DEGs you want to correlate
expr_table1 <- AverageExpression(pandey.ad, assays = "RNA", slot="counts")[[1]]
expr_table2 <- AverageExpression(final.embryo, assays = "RNA", slot="counts")[[1]]
DEgenes1 <- DEG.pandey.ad$`union_fdr<0.05`
DEgenes2 <- DEG.final.embryo$`union_fdr<0.05`

#if dataset #1 is a zebrafish dataset, then we need to convert the gene names to 1:1 mouse orthologue genes
#uncomment this section to do so:
load("data/output/merge_zebrafish/dre_mmu_1-1orts.RData")
expr_table1 <- subset(expr_table1, rownames(expr_table1) %in% mart.ort$Gene.name)
expr_table1 <- expr_table1 %>% rownames_to_column() %>% left_join(mart.ort, by = c("rowname" = "Gene.name"))
rownames(expr_table1) <- expr_table1$Mouse.gene.name
expr_table1 <- expr_table1[, -c(1, length(expr_table1))]
inbetween <- left_join(data.frame(DEgenes1), mart.ort, by=c("DEgenes_sp1" = "Gene.name"))
inbetween <- na.omit(inbetween)
DEgenes1 <- inbetween$Mouse.gene.name
rm(inbetween)

#execute correlation function
comp.intersect <- SpPermute(expr_table1, DEgenes1, expr_table2, DEgenes2, nPermutations=1000, genes.use= "intersect",corr.method="spearman")
#get correlation matrix
comp_table.intersect <- t(comp.intersect[[1]][1:ncol(expr_table1),(ncol(expr_table1)+1):nrow(comp.intersect[[1]])])
#change order of correlation matrix (take the commented vectors behind the DEG results)
comp_table.intersect <- comp_table.intersect[c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
                                             c("Ad_Hb01", "Ad_Hb02A", "Ad_Hb02B", "Ad_Hb04", "Ad_Hb05", "Ad_Hb06", "Ad_Hb08", "Ad_Hb09", "Ad_Hb10", "Ad_Hb11", "Ad_Hb13", "Ad_Hb14", "Ad_VHb01", "Ad_VHb02", "Ad_VHb03", "Ad_VHb04", "Ad_Hb16")]
#if the developmental Hb is correlated, then we add the "Cluster" prefix.
#uncomment line below (COL = expr_table1  ~ ROW = expr_table2) 
rownames(comp_table.intersect) <- paste("Cluster", rownames(comp_table.intersect), sep = "_")

#get signficance matrix
p_table.intersect <- t(comp.intersect[[3]][1:ncol(expr_table1),(ncol(expr_table1)+1):nrow(comp.intersect[[1]])])
#change order of significance matrix (take the commented vectors behind the DEG results)
p_table.intersect <- p_table.intersect[c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
                                       c("Ad_Hb01", "Ad_Hb02A", "Ad_Hb02B", "Ad_Hb04", "Ad_Hb05", "Ad_Hb06", "Ad_Hb08", "Ad_Hb09", "Ad_Hb10", "Ad_Hb11", "Ad_Hb13", "Ad_Hb14", "Ad_VHb01", "Ad_VHb02", "Ad_VHb03", "Ad_VHb04", "Ad_Hb16")]
#if the developmental Hb is correlated, then we add the "Cluster" prefix.
#uncomment line below (COL = expr_table1  ~ ROW = expr_table2) 
rownames(p_table.intersect) <- paste("Cluster", rownames(p_table.intersect), sep = "_")

#change minor details in significance matrix for nice plot results
p_vals <- p_table.intersect
p_table.intersect <- (-p_table.intersect)
p_vals[p_vals>0.05] <- NA
p_vals[p_vals<=0.05] <- "."

#plot correlations in correlation matrix ordered by hierarchical clustering
cols <- colorRampPalette(c("darkblue", "white","darkred"))

heatmaply_cor(comp_table.intersect, cellnote=p_vals, colors = cols(200), node_type="scatter", 
              label_names=c("Embryo", "Pandey", "Correlation"), limits = c(min(comp_table.intersect), max(comp_table.intersect)), 
              cellnote_color="black", cellnote_size=8, cellnote_textposition="middle center", point_size_mat=p_table.intersect,
              point_size_name="p-value")
