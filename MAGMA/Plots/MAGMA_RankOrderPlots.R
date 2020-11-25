##### rank-order correlations between two performed MAGMA protocols #####
# author: Juliska E Boer
# date: 04 Nov 2020

#load packages
setwd("E:/")
library(dplyr)
library(ggplot2)

#read in function to perform rank-order correlations
source("MAGMA_RankOrderFunct.R")

#calculate and plot rank-order correlations between two performed protocols per GWAS tested

#Major Depressive Disorder (MDD)
#read in results from each protocol for each dataset (excluding merged mouse Hb dataset)
sMDD_kyoko <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sMDD_skene <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_ANOVAstuber.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

wMDD_kyoko <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wMDD_skene <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_ANOVAwallace.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

eMDD_kyoko <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eMDD_skene <- read.table("data/output/GWAS/MDD2-18/Jul2020_GWAS_MDD_ANOVAembryo.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

tmMDD_kyoko <- read.table("data/output/GWAS/MDD2-18/TabulaMuris_MDD_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmMDD_skene <- read.table("data/output/GWAS/MDD2-18/TabulaMuris_MDD_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

#give correct names
eMDD_kyoko$Dataset <- rep("Developmental Hb", times=14)
sMDD_kyoko$Dataset <- rep("Hashikawa", times=12)
wMDD_kyoko$Dataset <- rep("Wallace", times=9)
tmMDD_kyoko$Dataset <- rep("Tabula Muris", times=9)
#rbind() files into one file for the first protocol by Nathan Skene
all_skene <- rbind(sMDD_skene, wMDD_skene, eMDD_skene, tmMDD_skene)
rm(eMDD_skene, sMDD_skene, wMDD_skene, tmMDD_skene)

#rbind() files into one file for the second protocl by Kyoko Watanabe
all_kyoko <- rbind(sMDD_kyoko, wMDD_kyoko, eMDD_kyoko, tmMDD_kyoko)
rm(eMDD_kyoko, sMDD_kyoko, wMDD_kyoko, tmMDD_kyoko)

#perform correlations
corr_allMDD <- calc_rank_cor(all_kyoko, all_skene, all=TRUE)

#save the rs-value and the rankings
rs_allMDD <- corr_allMDD[[1]]
ranks_allMDD <- corr_allMDD[[2]]

#plot result for MDD
ggplot(ranks_allMDD, aes(y=Rkyoko, x=Rskene)) + geom_point(aes(colour = Dataset), size=3) + geom_smooth(method='lm', se=FALSE, color="black") + 
  xlab("Protocol 2 (Skene)") + ylab("Protocol 1 (Kyoko)") + annotate(geom="text", x=32, y=22, label=round(rs_allMDD,digits=4),color="black", size=5) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + labs(colour = "Dataset") +
  scale_color_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/GWAS/MAGMAmethod_corr_MDD.pdf", height=5, width=7)

#Schizophrenia
#read in results from each protocol for each dataset (excluding merged mouse Hb dataset)
sSCZ_kyoko <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sSCZ_skene <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_ANOVAstuber.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

wSCZ_kyoko <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wSCZ_skene <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_ANOVAwallace.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

eSCZ_kyoko <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eSCZ_skene <- read.table("data/output/GWAS/SCZ/Jul2020_GWAS_SCZ_ANOVAembryo.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

tmSCZ_kyoko <- read.table("data/output/GWAS/SCZ/TabulaMuris_SCZ_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmSCZ_skene <- read.table("data/output/GWAS/SCZ/TabulaMuris_SCZ_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

#give correct names
eSCZ_kyoko$Dataset <- rep("Developmental Hb", times=14)
sSCZ_kyoko$Dataset <- rep("Hashikawa", times=12)
wSCZ_kyoko$Dataset <- rep("Wallace", times=9)
tmSCZ_kyoko$Dataset <- rep("Tabula Muris", times=9)
#rbind() files into one file for the first protocol by Nathan Skene
all_skene <- rbind(sSCZ_skene, wSCZ_skene, eSCZ_skene, tmSCZ_skene)
rm(eSCZ_skene, sSCZ_skene, wSCZ_skene, tmSCZ_skene)

#rbind() files into one file for the second protocl by Kyoko Watanabe
all_kyoko <- rbind(sSCZ_kyoko, wSCZ_kyoko, eSCZ_kyoko, tmSCZ_kyoko)
rm(eSCZ_kyoko, sSCZ_kyoko, wSCZ_kyoko, tmSCZ_kyoko)

#perform correlations
corr_allSCZ <- calc_rank_cor(all_kyoko, all_skene, all=TRUE)

#save the rs-value and the rankings
rs_allSCZ <- corr_allSCZ[[1]]
ranks_allSCZ <- corr_allSCZ[[2]]

#plot result for SCZ
ggplot(ranks_allSCZ, aes(y=Rkyoko, x=Rskene)) + geom_point(aes(colour = Dataset), size=3) + geom_smooth(method='lm', se=FALSE, color="black") + 
  xlab("Protocol 2 (Skene)") + ylab("Protocol 1 (Kyoko)") + annotate(geom="text", x=35, y=24, label=round(rs_allSCZ,digits=4),color="black", size=5) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + labs(colour = "Dataset") + 
  scale_color_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/GWAS/MAGMAmethod_corr_SCZ.pdf", height=5, width=7)

#Bipolar disorder
#read in results from each protocol for each dataset (excluding merged mouse Hb dataset)
sBIP_kyoko <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sBIP_skene <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_ANOVAstuber.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

wBIP_kyoko <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wBIP_skene <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_ANOVAwallace.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

eBIP_kyoko <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eBIP_skene <- read.table("data/output/GWAS/BIP/Jun2020_GWAS_BIP_ANOVAembryo.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

tmBIP_kyoko <- read.table("data/output/GWAS/BIP/TabulaMuris_BIP_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmBIP_skene <- read.table("data/output/GWAS/BIP/TabulaMuris_BIP_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

#give correct names
eBIP_kyoko$Dataset <- rep("Developmental Hb", times=14)
sBIP_kyoko$Dataset <- rep("Hashikawa", times=12)
wBIP_kyoko$Dataset <- rep("Wallace", times=9)
tmBIP_kyoko$Dataset <- rep("Tabula Muris", times=9)
#rbind() files into one file for the first protocol by Nathan Skene
all_skene <- rbind(sBIP_skene, wBIP_skene, eBIP_skene, tmBIP_skene)
rm(eBIP_skene, sBIP_skene, wBIP_skene, tmBIP_skene)

#rbind() files into one file for the second protocl by Kyoko Watanabe
all_kyoko <- rbind(sBIP_kyoko, wBIP_kyoko, eBIP_kyoko, tmBIP_kyoko)
rm(eBIP_kyoko, sBIP_kyoko, wBIP_kyoko, tmBIP_kyoko)

#perform correlations
corr_allBIP <- calc_rank_cor(all_kyoko, all_skene, all=TRUE)

#save the rs-value and the rankings
rs_allBIP <- corr_allBIP[[1]]
ranks_allBIP <- corr_allBIP[[2]]

#plot result for BIP
ggplot(ranks_allBIP, aes(y=Rkyoko, x=Rskene)) + geom_point(aes(colour = Dataset), size=3) + geom_smooth(method='lm', se=FALSE, color="black") + 
  xlab("Protocol 2 (Skene)") + ylab("Protocol 1 (Kyoko)") + annotate(geom="text", x=40, y=23, label=round(rs_allBIP,digits=4),color="black", size=5) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + labs(colour = "Dataset") + 
  scale_color_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/GWAS/MAGMAmethod_corr_BIP.pdf", height=5, width=7)

#Body-Mass Index
#read in results from each protocol for each dataset (excluding merged mouse Hb dataset)
sBMI_kyoko <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sBMI_skene <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_ANOVAstuber.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

wBMI_kyoko <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wBMI_skene <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_ANOVAwallace.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

eBMI_kyoko <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eBMI_skene <- read.table("data/output/GWAS/BMI/Jun2020_GWAS_BMI_ANOVAembryo.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

tmBMI_kyoko <- read.table("data/output/GWAS/BMI/TabulaMuris_BMI_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmBMI_skene <- read.table("data/output/GWAS/BMI/TabulaMuris_BMI_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

#give correct names
eBMI_kyoko$Dataset <- rep("Developmental Hb", times=14)
sBMI_kyoko$Dataset <- rep("Hashikawa", times=12)
wBMI_kyoko$Dataset <- rep("Wallace", times=9)
tmBMI_kyoko$Dataset <- rep("Tabula Muris", times=9)
#rbind() files into one file for the first protocol by Nathan Skene
all_skene <- rbind(sBMI_skene, wBMI_skene, eBMI_skene, tmBMI_skene)
rm(eBMI_skene, sBMI_skene, wBMI_skene, tmBMI_skene)

#rbind() files into one file for the second protocl by Kyoko Watanabe
all_kyoko <- rbind(sBMI_kyoko, wBMI_kyoko, eBMI_kyoko, tmBMI_kyoko)
rm(eBMI_kyoko, sBMI_kyoko, wBMI_kyoko, tmBMI_kyoko)

#perform correlations
corr_allBMI <- calc_rank_cor(all_kyoko, all_skene, all=TRUE)

#save the rs-value and the rankings
rs_allBMI <- corr_allBMI[[1]]
ranks_allBMI <- corr_allBMI[[2]]

#plot result for BMI
ggplot(ranks_allBMI, aes(y=Rkyoko, x=Rskene)) + geom_point(aes(colour = Dataset), size=3) + geom_smooth(method='lm', se=FALSE, color="black") + 
  xlab("Protocol 2 (Skene)") + ylab("Protocol 1 (Kyoko)") + annotate(geom="text", x=30, y=23, label=round(rs_allBMI,digits=4),color="black", size=5) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + labs(colour = "Dataset") + 
  scale_color_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/GWAS/MAGMAmethod_corr_BMI.pdf", height=5, width=7)

#Height
#read in results from each protocol for each dataset (excluding merged mouse Hb dataset)
sHeight_kyoko <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sHeight_skene <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_ANOVAstuber.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

wHeight_kyoko <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wHeight_skene <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_ANOVAwallace.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

eHeight_kyoko <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eHeight_skene <- read.table("data/output/GWAS/Height/Jun2020_GWAS_Height_ANOVAembryo.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

tmHeight_kyoko <- read.table("data/output/GWAS/Height/TabulaMuris_Height_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmHeight_skene <- read.table("data/output/GWAS/Height/TabulaMuris_Height_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]

#give correct names
eHeight_kyoko$Dataset <- rep("Developmental Hb", times=14)
sHeight_kyoko$Dataset <- rep("Hashikawa", times=12)
wHeight_kyoko$Dataset <- rep("Wallace", times=9)
tmHeight_kyoko$Dataset <- rep("Tabula Muris", times=9)
#rbind() files into one file for the first protocol by Nathan Skene
all_skene <- rbind(sHeight_skene, wHeight_skene, eHeight_skene, tmHeight_skene)
rm(eHeight_skene, sHeight_skene, wHeight_skene, tmHeight_skene)

#rbind() files into one file for the second protocl by Kyoko Watanabe
all_kyoko <- rbind(sHeight_kyoko, wHeight_kyoko, eHeight_kyoko, tmHeight_kyoko)
rm(eHeight_kyoko, sHeight_kyoko, wHeight_kyoko, tmHeight_kyoko)

#perform correlations
corr_allHeight <- calc_rank_cor(all_kyoko, all_skene, all=TRUE)

#save the rs-value and the rankings
rs_allHeight <- corr_allHeight[[1]]
ranks_allHeight <- corr_allHeight[[2]]

#plot result for Height
ggplot(ranks_allHeight, aes(y=Rkyoko, x=Rskene)) + geom_point(aes(colour = Dataset), size=3) + geom_smooth(method='lm', se=FALSE, color="black") + 
  xlab("Protocol 2 (Skene)") + ylab("Protocol 1 (Kyoko)") + annotate(geom="text", x=41, y=32, label=round(rs_allHeight,digits=4),color="black", size=5) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=14)) + labs(colour = "Dataset") + 
  scale_color_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/GWAS/MAGMAmethod_corr_Height.pdf", height=5, width=7)

