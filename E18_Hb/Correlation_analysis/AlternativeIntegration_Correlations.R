##### alternative integration: correlation analysis #####
# author: Juliska E Boer
# edidited: L.L. van de Haar
# date: 21 01 2021
setwd("/Volumes/HD_LLH/python/STAR_E18_wt/FINAL/Correlation")

#load packages
library(MAST)
library(Seurat)
library(tidyverse)
library(heatmaply)

#load functions for identifying differentially expressed genes (DEGs) using MAST and for calculating the correlations
source("Correlation_CalcDEG.R")
source("Correlation_CalcCorr.R")

#read in the data
load("../data/output/E18WT_Scanpy_Seurat_obj.RData") #E18WT
load("../data/output/final.E18Brn3a_Scanpy_Seurat_obj.RData") #E18 Brn3a

#determine DEGs for each dataset, and save information (this function takes some time)
DEG.final.E18WT <- DE_Gene_Union(final.E18WT, levels(final.E18WT@active.ident), data.info.col = 1:9)
save(DEG.final.E18WT, file = "../data/output/DEG-final-E18WT.RData")

DEG.final.E18Brn3a <- DE_Gene_Union(final.E18Brn3a, levels(final.E18Brn3a@active.ident), data.info.col = 1:5)
save(DEG.final.E18Brn3a, file = "../data/output/DEG-final-E18Brn3a.RData")

#prepare input for correlation function
#average expression matrices filtered on DEG

#load DEG result (vectors behind the objects are the cluster names, these are used for plotting)
load("../data/output/DEG-final-E18WT.RData") #c("Ad_Hb01", "Ad_Hb02A", "Ad_Hb02B", "Ad_Hb04", "Ad_Hb05", "Ad_Hb06", "Ad_Hb08", "Ad_Hb09", "Ad_Hb10", "Ad_Hb11", "Ad_Hb13", "Ad_Hb14", "Ad_VHb01", "Ad_VHb02", "Ad_VHb03", "Ad_VHb04", "Ad_Hb16")
load("../data/output/DEG-final-E18Brn3a.RData") #c("La_Hb01", "La_Hb02", "La_Hb03", "La_Hb04", "La_Hb05", "La_Hb06", "La_Hb07", "La_Hb08", "La_Hb09", "La_Hb10", "La_Hb11", "La_Hb12", "La_Hb13", "La_Hb14", "La_Hb15", "Olf")

#set the datasets and DEGs you want to correlate
expr_table1 <- AverageExpression(final.E18WT, assays = "RNA", slot="counts")[[1]]
expr_table2 <- AverageExpression(final.E18Brn3a, assays = "RNA", slot="counts")[[1]]
DEgenes1 <- DEG.final.E18WT$`union_fdr<0.05`
DEgenes2 <- DEG.final.E18Brn3a$`union_fdr<0.05`

DEgenesSpecies1 <- DEgenes1
DEgenesSpecies2 <- DEgenes2
#ExpressionTable1 <- expr_table1
#ExpressionTable2 <- expr_table2

nDESp1 <- 2729
nDESp2 <- 512

#execute correlation function
comp.intersect <- SpPermute(expr_table1, DEgenes1, expr_table2, DEgenes2, nPermutations=1000, genes.use= "intersect",corr.method="spearman")
#get correlation matrix
comp_table.intersect <- t(comp.intersect[[1]][1:ncol(expr_table1),(ncol(expr_table1)+1):nrow(comp.intersect[[1]])])
#change order of correlation matrix (take the commented vectors behind the DEG results)
comp_table.intersect <- comp_table.intersect[c("0", "1", "6", "2", "4", "3", "5", "7", "8"),
                                             c("0", "3", "4", "1", "2")]
#if the developmental Hb is correlated, then we add the "Cluster" prefix.
#uncomment line below (COL = expr_table1  ~ ROW = expr_table2) 
rownames(comp_table.intersect) <- paste("Cluster", rownames(comp_table.intersect), sep = "_")

#get signficance matrix
p_table.intersect <- t(comp.intersect[[3]][1:ncol(expr_table1),(ncol(expr_table1)+1):nrow(comp.intersect[[1]])])
#change order of significance matrix (take the commented vectors behind the DEG results)
p_table.intersect <- p_table.intersect[c("0", "1", "6", "2", "4", "3", "5", "7", "8"),
                                       c("0", "3", "4", "1", "2")]
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
              width = 20, height = 8, point_size_name="p-value")
