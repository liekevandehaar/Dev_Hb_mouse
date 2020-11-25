##### preprocess scRNAseq data for MAGMA #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(data.table)
library(dplyr)
library(Seurat)
library(limma)
library(stats)
library(tibble)

source("GitHub/MAGMA/Preprocess/scRNAseq/MAGMA_PreprocessFunct.R")

#read in files holding 1:1 mouse orthologs
load("data/input/mm2hs.RData")

#read in Seurat objects for selection of informative genes using ANOVA
load("data/output/DevelopmentalHb/Embryo_Scanpy_Seurat_obj.RData")
stuber <- readRDS("data/input/Habenula_neuron_Seurat_object.rds")
load("data/output/ExternalDatasets/Wallace_Seurat_obj.RData")
load("data/output/ExternalDatasets/TabulaMuris_Seurat_object.RData")
load("data/output/ExternalDatasets/Merged_DevHB_H_W.RData")
load("data/output/ExternalDatasets/DevHypothalamus-samp_obj.RData")
load("data/output/ExternalDatasets/DevAdMouse-LaManno2016_obj.RData")

#preprocess data into two files, one for each protocol

#read in average expression data
#developmental Hb
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Jun2020_GWAS_embryo.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Jun2020_GWAS_embryo.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- final.embryo
#Hashikawa, et al (2020)
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Jun2020_GWAS_stuber.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Jun2020_GWAS_stuber.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- stuber
#Wallace, et al (2020)
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Jun2020_GWAS_wallace.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Jun2020_GWAS_wallace.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- wallace_hb
#Tabula Muris, et al (2018)
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/TabulaMuris_avgexpr.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/TabulaMuris_avgexpr.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- brain_tiss_filtered
#Merged mouse Hb dataset
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Oct2020_GWAS_merged.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Oct2020_GWAS_merged.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- merged
#Developmental Hypothalamus dataset (2020)
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Nov2020_GWAS_devhypo-samp.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Nov2020_GWAS_devhypo-samp.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- subhypo
#Embryonic and Adult Mouse dataset, La Manno (2016)
mmugenes <- data.frame(mm.symbol=fread("data/output/MAGMA/Nov2020_GWAS_lamanno2016.csv", header=T)[[1]], stringsAsFactors = F)
mean <- read.table("data/output/MAGMA/Nov2020_GWAS_lamanno2016.csv", header=TRUE, sep=" ", dec=".", row.names=1)
object <- lamanno

#conversion of mouse genes to 1:1 human orthologs
hsgenes <- MMUgenes_toHSgenes(mmugenes)

#uncomment two lines below for first protocol by Nathan Skene. only use informative genes
ANOVAres <- ANOVA_infogenes(object)
mean <- mean[intersect(rownames(mean), ANOVAres),]

#uncomment line below for first protocol by Nathan Skene. create gene-cell type specificity matrix
mean <- sweep(mean,MARGIN=1,FUN="/",STATS=rowSums(mean))

#continue with conversion of mouse genes to 1:1 human orthologs
exp <- rownames_to_column(mean, "GENE")
exp$GENE <- hsgenes$hs.ensg[match(exp$GENE, hsgenes$mm.symbol)]
exp <- na.omit(exp)

#subsample to average number of genes used in adult datasets
exp <- exp[sample(nrow(exp), 11304, replace=FALSE),]

#uncomment line below for second protocol by Kyoko Watanabe. add column holding average expression per gene
#exp$AVERAGE <- apply(exp[,-1], 1, mean)

#write input file for MAGMA analysis using first protocol by Nathan Skene
write.table(exp, "data/input/Nov2020_GWAS_TabulaM-samp_ANOVA.txt", quote=F, row.names=F, sep="\t") #CHANGE FILE NAME TO DATASET
#write input file for MAGMA analysis using second protocol by Kyoko Watanabe
#write.table(exp, "data/input/Nov2020_GWAS_DevHb_AVG.txt", quote=F, row.names=F, sep="\t") #CHANGE FILE NAME TO DATASET