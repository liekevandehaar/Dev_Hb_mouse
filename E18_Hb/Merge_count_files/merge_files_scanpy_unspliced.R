# How to merge files...
#install.packages("rio")

##Set working directory
setwd("/Volumes/HD_LLH/python/STAR_E18_wt/FINAL")

##START
# GLASTEGFmixSEZ
#import files
L1 <- read.csv("./data/cout/Habenula_E18_L1_velo_cbc_trimmed_star_unspliced.coutt.tsv",sep="\t",header=TRUE)
L2 <- read.csv("./data/cout/Habenula_E18_L2_velo_cbc_trimmed_star_unspliced.coutt.tsv",sep="\t",header=TRUE)
R1 <- read.csv("./data/cout/Habenula_E18_R1_velo_cbc_trimmed_star_unspliced.coutt.tsv",sep="\t",header=TRUE)
R2 <- read.csv("./data/cout/Habenula_E18_R2_velo_cbc_trimmed_star_unspliced.coutt.tsv",sep="\t",header=TRUE)
LR <- read.csv("./data/cout/Habenula_E18_Pl5_velo_cbc_trimmed_star_unspliced.coutt.tsv",sep="\t",header=TRUE)

#modify row names
rownames(L1) = L1[,1]
L1 = L1[,-1]
rownames(L2) = L2[,1]
L2 = L2[,-1]
rownames(R1) = R1[,1]
R1 = R1[,-1]
rownames(R2) = R2[,1]
R2 = R2[,-1]
rownames(LR) = LR[,1]
LR = LR[,-1]


#modify columns
colnames(L1) <- paste("Habenula_001___L1__left_", seq(1, 384), sep = "")
colnames(L2) <- paste("Habenula_003___L2__left_", seq(1, 384), sep = "")
colnames(R1) <- paste("Habenula_002___R1__right_", seq(1, 384), sep = "")
colnames(R2) <- paste("Habenula_004___R2__right_", seq(1, 384), sep = "")
colnames(LR) <- paste("Habenula_005___LR__mix_", seq(1, 384), sep = "")


#merge unspliced version
Series_Jun2019_merge <- merge(L1,L2, by = 'row.names', all = TRUE)
rownames(Series_Jun2019_merge) = Series_Jun2019_merge[,1]
Series_Jun2019_merge = Series_Jun2019_merge[,-1]

Series_Jun2019_merge <- merge(Series_Jun2019_merge,R1, by = 'row.names', all = TRUE)
rownames(Series_Jun2019_merge) = Series_Jun2019_merge[,1]
Series_Jun2019_merge = Series_Jun2019_merge[,-1]

Series_Jun2019_merge <- merge(Series_Jun2019_merge,R2, by = 'row.names', all = TRUE)
rownames(Series_Jun2019_merge) = Series_Jun2019_merge[,1]
Series_Jun2019_merge = Series_Jun2019_merge[,-1]

Series_Jun2019_merge <- merge(Series_Jun2019_merge,LR, by = 'row.names', all = TRUE)
rownames(Series_Jun2019_merge) = Series_Jun2019_merge[,1]
Series_Jun2019_merge = Series_Jun2019_merge[,-1]

#readjust
Series_Jun2019_merge[is.na(Series_Jun2019_merge)] <- 0

#save the plot
pdf("SeriesJun2019_merge unique_reads combined.pdf")
barplot(colSums(Series_Jun2019_merge), main = 'Series_Jun2019_merge')
dev.off()

#save
library(rio)
export(Series_Jun2019_merge, "merged_E18_Habenula.coutt_unspliced_2.tsv", row.names = T)
