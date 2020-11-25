##### classify developmental Hb cells to previously annotated clusters using a RandomForest classifier #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(randomForest)
library(Seurat)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)

#load function for plotting confusion matrix
source("GitHub/AlternativeIntegration/RandomForest/RandomForest_ConfusionMatrix.R")

#load data
load("data/output/DevelopmentalHb/Embryo_Scanpy_Seurat_obj.RData")
stuber <- readRDS("data/input/Habenula_neuron_Seurat_object.rds")
load("data/output/ExternalDatasets/Wallace_Seurat_obj.RData")
load("data/output/ExternalDatasets/Pandey_Seurat_obj.RData")
load("data/output/ExternalDatasets/PandeyAdult_Seurat_obj.RData")

#select for highly variable genes (HVG)
final.embryo <- FindVariableFeatures(final.embryo)
stuber <- FindVariableFeatures(stuber, assay="RNA")
wallace_hb <- FindVariableFeatures(wallace_hb)
pandey.id <- FindVariableFeatures(pandey.id)
pandey.ad <- FindVariableFeatures(pandey.ad)

#load DEG results
load("data/output/DevelopmentalHb/DEG-final-embryo.RData") #c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
load("data/output/ExternalDatasets/DEG_Stuber.RData") #c("MHb1", "MHb2", "MHb3", "MHb4", "MHb5", "MHb6", "LHb1", "LHb2", "LHb3", "LHb4", "LHb5", "LHb6")
load("data/output/ExternalDatasets/DEG-wallace.RData") #c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)") c("Oval.Medial.LHb1.", "Marginal.LHb2.", "Lateral.LHb3.", "HBX.LHb4.", "Ventral2.3.MHb1.", "Ventrolateral.MHb2.", "Lateral.MHb3.", "Dorsal.MHb4.", "Superior.MHb5.")
load("data/output/ExternalDatasets/DEG_pandey.RData") #c("La_Hb01", "La_Hb02", "La_Hb03", "La_Hb04", "La_Hb05", "La_Hb06", "La_Hb07", "La_Hb08", "La_Hb09", "La_Hb10", "La_Hb11", "La_Hb12", "La_Hb13", "La_Hb14", "La_Hb15", "Olf")
load("data/output/ExternalDatasets/DEG_pandey_adult.RData") #c("Ad_Hb01", "Ad_Hb02A", "Ad_Hb02B", "Ad_Hb04", "Ad_Hb05", "Ad_Hb06", "Ad_Hb08", "Ad_Hb09", "Ad_Hb10", "Ad_Hb11", "Ad_Hb13", "Ad_Hb14", "Ad_VHb01", "Ad_VHb02", "Ad_VHb03", "Ad_VHb04", "Ad_Hb16")

#set the datasets for RF classification
#annotated dataset for training the classifier
RFobject = wallace_hb
#dataset used for predicting with the classifier
Predobject = final.embryo

#if the annotated dataset is a zebrafish dataset, then we need to convert the gene names to 1:1 mouse orthologue genes
#uncomment this section to do so:
#load("data/input/dre_mmu_1-1orts.RData")
#pandey_df = data.frame(RFobject@assays$RNA@counts)

#pandey_df <- subset(pandey_df, rownames(pandey_df) %in% mart.ort$Gene.name)
#pandey_df <- pandey_df %>% rownames_to_column() %>% left_join(mart.ort, by = c("rowname" = "Gene.name"))
#rownames(pandey_df) <- pandey_df$Mouse.gene.name
#pandey_df <- pandey_df[, -c(1, length(pandey_df))]
#pandey_matrix <- as.matrix(pandey_df)

#inbetween <- left_join(data.frame(RFobject@assays$RNA@var.features), mart.ort, by=c("RFobject.assays.RNA.var.features" = "Gene.name"))
#inbetween <- na.omit(inbetween)
#var.feat_pandey <- inbetween$Mouse.gene.name
#rm(inbetween)

#uncomment below when using the HVG selected by Seurat with a zebrafish dataset
#genes.use = intersect(var.feat_pandey, Predobject@assays$RNA@var.features)

#uncomment below when using the HVG selected by Seurat
genes.use = intersect(RFobject@assays$RNA@var.features, Predobject@assays$RNA@var.features)

#uncomment below when using the DEGs selected by MAST
#genes.use = intersect(DE_Gene_Stuber$`union_fdr<0.05`, DEG.wallace$`union_fdr<0.05`)

#create training and validation set
training.set = c(); valid.set=c()
training.label = c(); valid.label=c();
for (i in levels(RFobject@active.ident)){
  cells.in.clust = WhichCells(RFobject,idents=i);
  n = min(500, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  valid.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); valid.set=c(valid.set,valid.temp)
  training.label = c(training.label, rep(i,length(train.temp))); valid.label = c(valid.label, rep(i, length(valid.temp)));
  training.set = str_replace(training.set, "-", ".")
  valid.set = str_replace(valid.set, "-", ".")
}

#uncomment line below for zebrafish dataset
#rf_data = pandey_matrix[genes.use,]
rf_data = as.matrix(RFobject@assays$RNA@data[genes.use,])

#uncomment line below if you want to perform extra sampling
#sampsizes = rep(min(as.vector(table(training.label))),length(as.vector(table(training.label))))

#train RandomForest classifier on 70% of the annotated dataset
rf_output=randomForest(x=t(rf_data[,training.set]), y=factor(training.label), ntree=1000, mtry=sqrt(length(genes.use)))
#let the RF predict on the remaining 30% of the annotated dataset
test.predict = predict(rf_output,t(rf_data[,test.set]), )
Conf_test = table(test.label,test.predict)
rf_output

#use the trained RF to classify the developmental Hb cells
pred_data = as.matrix(Predobject@assays$RNA@counts[genes.use,])
probabilities <- predict(rf_output,t(pred_data), type="prob")
#classify cells by the number of votes. NOT majority of votes but only if it has > 25% of votes
#create confusion matrix
cols <- colnames(probabilities)[max.col(probabilities) * NA^!rowSums(probabilities > 0.2) > 0]
cols[is.na(cols)] <- 'unassigned'
pred_idents = factor(Predobject@active.ident)
predictions = data.frame(unclass(table(pred_idents, cols)))

#if cells are not mapped to a classification, it is not included in the confusion matrix. we add these classes to prevent confusion
#missing clusters to Hashikawa (DevHb)
#predictions$LHb4 = rep(0,14)
#predictions$MHb1 = rep(0,14)
#predictions$MHb4 = rep(0,14)

#missing clusters to Adult Pandey (DevHb)
#predictions$Ad_Hb02A = rep(0,14)
#predictions$Ad_Hb05 = rep(0,14)
#predictions$Ad_Hb08 = rep(0,14)
#predictions$Ad_Hb10 = rep(0,14)
#predictions$Ad_Hb11 = rep(0,14)
#predictions$Ad_VHb01 = rep(0,14)
#predictions$Ad_VHb02 = rep(0,14)
#predictions$Ad_VHb03 = rep(0,14)
#predictions$Ad_VHb04 = rep(0,14)

#missing clusters to Larval Pandey (DevHb)
#predictions$La_Hb03 = rep(0,14)
#predictions$La_Hb05 = rep(0,14)
#predictions$La_Hb07 = rep(0,14)
#predictions$La_Hb08 = rep(0,14)
#predictions$La_Hb09 = rep(0,14)
#predictions$La_Hb10 = rep(0,14)
#predictions$La_Hb11 = rep(0,14)
#predictions$La_Hb12 = rep(0,14)
#predictions$La_Hb14 = rep(0,14)
#predictions$La_Hb15 = rep(0,14)
#predictions$Olf = rep(0,14)

#change order of confusion matrix (take the commented vectors behind the DEG results)
predictions <- predictions[c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
                           c("La_Hb01", "La_Hb02", "La_Hb03", "La_Hb04", "La_Hb05", "La_Hb06", "La_Hb07", "La_Hb08", "La_Hb09", "La_Hb10", 
                             "La_Hb11", "La_Hb12", "La_Hb13", "La_Hb14", "La_Hb15", "Olf", "unassigned")]
#if the developmental Hb is classified, then we add the "Cluster" prefix. uncomment line below
rownames(predictions) <- paste("Cluster", rownames(predictions), sep = "_")
predictions <- predictions %>% relocate(which(colnames(predictions)=="unassigned"), .after = last_col())

#plot confusion matrix
plotConfusionMatrix(predictions,row.scale=TRUE, max.size = 12, xlab.use="Larval - Pandey et al 2018", ylab.use="Developmental Hb")
