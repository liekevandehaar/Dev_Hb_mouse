##### plotting resulting p-values after MAGMA analysis #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(ggplot2)

#read in the results from each GWAS for each dataset, and rbind() into one dataframe in melted format
#change file names for protocol used

#Developmental Hb
eMDD <- read.table("data/output/MAGMA/MDD2-18/Jul2020_GWAS_MDD_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eMDD$GWAS <- c("MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018",
               "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018")
eMDD$VARIABLE <- c("1", "12", "0", "4", "8", "13", "6", "11", "3", "2", "5", "7", "10", "9")
#eMDD$VARIABLE <- c("early Hb2[1]", "PC2[12]", "early Hb1[0]", "early dMHb[4]", "iHb2[8]", "PC1[13]", "iHb1[6]",
 #                  "LHb[11]", "vMHb2[3]", "vMHb1[2]", "dMHb2[5]", "dMHb3[7]", "dMHb1[10]", "vMHb3[9]")
eBMI <- read.table("data/output/MAGMA/BMI/Jun2020_GWAS_BMI_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eBMI$GWAS <- c("BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015",
               "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015")
eBMI$VARIABLE <- c("1", "12", "0", "4", "8", "13", "6", "11", "3", "2", "5", "7", "10", "9")
#eBMI$VARIABLE <- c("early Hb2[1]", "PC2[12]", "early Hb1[0]", "early dMHb[4]", "iHb2[8]", "PC1[13]", "iHb1[6]",
 #                  "LHb[11]", "vMHb2[3]", "vMHb1[2]", "dMHb2[5]", "dMHb3[7]", "dMHb1[10]", "vMHb3[9]")
eSCZ <- read.table("data/output/MAGMA/SCZ/Jul2020_GWAS_SCZ_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eSCZ$GWAS <- c("SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018",
               "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018")
eSCZ$VARIABLE <- c("1", "12", "0", "4", "8", "13", "6", "11", "3", "2", "5", "7", "10", "9")
#eSCZ$VARIABLE <- c("early Hb2[1]", "PC2[12]", "early Hb1[0]", "early dMHb[4]", "iHb2[8]", "PC1[13]", "iHb1[6]",
 #                  "LHb[11]", "vMHb2[3]", "vMHb1[2]", "dMHb2[5]", "dMHb3[7]", "dMHb1[10]", "vMHb3[9]")
eBIP <- read.table("data/output/MAGMA/BIP/Jun2020_GWAS_BIP_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eBIP$GWAS <- c("BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016",
               "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016")
eBIP$VARIABLE <- c("1", "12", "0", "4", "8", "13", "6", "11", "3", "2", "5", "7", "10", "9")
#eBIP$VARIABLE <- c("early Hb2[1]", "PC2[12]", "early Hb1[0]", "early dMHb[4]", "iHb2[8]", "PC1[13]", "iHb1[6]",
 #                  "LHb[11]", "vMHb2[3]", "vMHb1[2]", "dMHb2[5]", "dMHb3[7]", "dMHb1[10]", "vMHb3[9]")
eHeight <- read.table("data/output/MAGMA/Height/Jun2020_GWAS_Height_embryo_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
eHeight$GWAS <- c("Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015",
                  "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015")
eHeight$VARIABLE <- c("1", "12", "0", "4", "8", "13", "6", "11", "3", "2", "5", "7", "10", "9")
#eHeight$VARIABLE <- c("early Hb2[1]", "PC2[12]", "early Hb1[0]", "early dMHb[4]", "iHb2[8]", "PC1[13]", "iHb1[6]",
 #                     "LHb[11]", "vMHb2[3]", "vMHb1[2]", "dMHb2[5]", "dMHb3[7]", "dMHb1[10]", "vMHb3[9]")
gsa_embryo <- rbind(eMDD, eBMI, eSCZ, eBIP, eHeight)
rm(eMDD, eBMI, eSCZ, eBIP, eHeight)

#Hashikawa, et al (2020)
sMDD <- read.table("data/output/MAGMA/MDD2-18/Jul2020_GWAS_MDD_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sMDD$GWAS <- c("MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018",
                     "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018")
sBMI <- read.table("data/output/MAGMA/BMI/Jun2020_GWAS_BMI_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sBMI$GWAS <- c("BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015",
                     "BMI 2015","BMI 2015","BMI 2015","BMI 2015")
sSCZ <- read.table("data/output/MAGMA/SCZ/Jul2020_GWAS_SCZ_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sSCZ$GWAS <- c("SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018",
               "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018")
sBIP <- read.table("data/output/MAGMA/BIP/Jun2020_GWAS_BIP_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sBIP$GWAS <- c("BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016",
               "BIP 2016","BIP 2016","BIP 2016","BIP 2016")
sHeight <- read.table("data/output/MAGMA/Height/Jun2020_GWAS_Height_stuber_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
sHeight$GWAS <- c("Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015",
                  "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015")
gsa_stuber <- rbind(sMDD, sBMI, sSCZ, sBIP, sHeight)
rm(sMDD, sBMI, sSCZ, sBIP, sHeight)

#Wallace, et al (2020)
wMDD <- read.table("data/output/MAGMA/MDD2-18/Jul2020_GWAS_MDD_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wMDD$GWAS <- c("MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018",
               "MDD2 2018")
wMDD$VARIABLE <- c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", 
                   "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")
wBMI <- read.table("data/output/MAGMA/BMI/Jun2020_GWAS_BMI_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wBMI$GWAS <- c("BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015",
               "BMI 2015")
wBMI$VARIABLE <- c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", 
                   "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")
wSCZ <- read.table("data/output/MAGMA/SCZ/Jul2020_GWAS_SCZ_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wSCZ$GWAS <- c("SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018",
               "SCZ 2018")
wSCZ$VARIABLE <- c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", 
                   "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")
wBIP <- read.table("data/output/MAGMA/BIP/Jun2020_GWAS_BIP_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wBIP$GWAS <- c("BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016",
               "BIP 2016")
wBIP$VARIABLE <- c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", 
                   "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")
wHeight <- read.table("data/output/MAGMA/Height/Jun2020_GWAS_Height_wallace_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
wHeight$GWAS <- c("Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015",
                  "Height 2015","Height 2015")
wHeight$VARIABLE <- c("Ventral2/3(MHb1)", "Ventrolateral(MHb2)", "Lateral(MHb3)", "Dorsal(MHb4)", "Superior(MHb5)", 
                   "Oval/Medial(LHb1)", "Marginal(LHb2)", "Lateral(LHb3)", "HBX(LHb4)")
gsa_wallace <- rbind(wMDD, wBMI, wSCZ, wBIP, wHeight)
rm(wMDD, wBMI, wSCZ, wBIP, wHeight)

#Tabula Muris, et al (2018)
tmMDD <- read.table("data/output/MAGMA/MDD2-18/TabulaMuris_MDD_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmMDD$GWAS <- c("MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018",
               "MDD2 2018")
tmBMI <- read.table("data/output/MAGMA/BMI/TabulaMuris_BMI_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmBMI$GWAS <- c("BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015",
               "BMI 2015")
tmSCZ <- read.table("data/output/MAGMA/SCZ/TabulaMuris_SCZ_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmSCZ$GWAS <- c("SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018",
               "SCZ 2018")
tmBIP <- read.table("data/output/MAGMA/BIP/TabulaMuris_BIP_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmBIP$GWAS <- c("BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016",
               "BIP 2016")
tmHeight <- read.table("data/output/MAGMA/Height/TabulaMuris_Height_AVG.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
tmHeight$GWAS <- c("Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015",
                  "Height 2015","Height 2015")
gsa_tabulam <- rbind(tmMDD, tmBMI, tmSCZ, tmBIP, tmHeight)
rm(tmMDD, tmBMI, tmSCZ, tmBIP, tmHeight)

#transform p-values to -log10 values
#give each file the dataset name
gsa_embryo$logP <- -log10(gsa_embryo$P)
gsa_embryo$Dataset <- rep("Developmental Hb", times=70)
gsa_stuber$logP <- -log10(gsa_stuber$P)
gsa_stuber$Dataset <- rep("Hashikawa", times=60)
gsa_wallace$logP <- -log10(gsa_wallace$P)
gsa_wallace$Dataset <- rep("Wallace", times=45)
gsa_tabulam$logP <- -log10(gsa_tabulam$P)
gsa_tabulam$Dataset <- rep("Tabula Muris", times=45)

#rbind() all files into one big file
gsa_all <- rbind(gsa_embryo, gsa_stuber, gsa_wallace, gsa_tabulam)
level_data <- gsa_all$VARIABLE
gsa_all$VARIABLE <- factor(gsa_all$VARIABLE, levels = unique(level_data))

#make plot
ggplot(gsa_all, aes(fill=Dataset, y=logP, x=VARIABLE)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~GWAS, nrow=5, scales = "free_y") + geom_hline(yintercept=-log10(0.05/(44*5)), linetype="dashed") +
  xlab("") + ylab("-log10(P)") + theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
                                       axis.text.y = element_text(color="black", size=12),
                                       axis.title = element_text(size=14),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       strip.text = element_text(face="bold", size=14),
                                       legend.text = element_text(size=12, color="black"),
                                       legend.title = element_text(size=13, color="black")) + 
  scale_fill_manual(values = c("#66c2a5", "#7690ca", "#e78ac3", "#fc8d62")) +
  ggsave("figures/MAGMA/GWAS_methodKyoko-num.pdf", height=8, width=12)


#Merged mouse Hb dataset

#read in results and rbind() into melted format
mMDD <- read.table("data/output/MAGMA/MDD2-18/Oct2020_GWAS_MDD2-18_merged_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
mMDD$GWAS <- c("MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018",
               "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018", "MDD2 2018","MDD2 2018","MDD2 2018","MDD2 2018", "MDD2 2018")
mMDD$VARIABLE <- c('dorsal/lateral MHb(0)', 'ventro/lateral MHb1(1)','vMHb1(2)', 'dMHb(3)', 'vMHb2(4)', 'early Hb1(5)','LHb subtype1(6)', 
                   'LHb subtype3(7)', 'early Hb2(8)', 'ventro/lateral MHb2(9)', 'lateral MHb(10)', 'superior MHb(11)', 'LHb subtype2(12)', 'HBX(13)',
                   'iHb1 (to MHb)(14)', 'iHb2 (to LHb)(15)', 'PC(16)')
mBMI <- read.table("data/output/MAGMA/BMI/Oct2020_GWAS_BMI_merged_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
mBMI$GWAS <- c("BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015","BMI 2015",
               "BMI 2015","BMI 2015","BMI 2015","BMI 2015", "BMI 2015","BMI 2015","BMI 2015","BMI 2015", "BMI 2015")
mBMI$VARIABLE <- c('dorsal/lateral MHb(0)', 'ventro/lateral MHb1(1)','vMHb1(2)', 'dMHb(3)', 'vMHb2(4)', 'early Hb1(5)','LHb subtype1(6)', 
                   'LHb subtype3(7)', 'early Hb2(8)', 'ventro/lateral MHb2(9)', 'lateral MHb(10)', 'superior MHb(11)', 'LHb subtype2(12)', 'HBX(13)',
                   'iHb1 (to MHb)(14)', 'iHb2 (to LHb)(15)', 'PC(16)')
mSCZ <- read.table("data/output/MAGMA/SCZ/Oct2020_GWAS_SCZ_merged_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
mSCZ$GWAS <- c("SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018",
               "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018", "SCZ 2018","SCZ 2018","SCZ 2018","SCZ 2018", "SCZ 2018")
mSCZ$VARIABLE <- c('dorsal/lateral MHb(0)', 'ventro/lateral MHb1(1)','vMHb1(2)', 'dMHb(3)', 'vMHb2(4)', 'early Hb1(5)','LHb subtype1(6)', 
                   'LHb subtype3(7)', 'early Hb2(8)', 'ventro/lateral MHb2(9)', 'lateral MHb(10)', 'superior MHb(11)', 'LHb subtype2(12)', 'HBX(13)',
                   'iHb1 (to MHb)(14)', 'iHb2 (to LHb)(15)', 'PC(16)')
mBIP <- read.table("data/output/MAGMA/BIP/Oct2020_GWAS_BIP_merged_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
mBIP$GWAS <- c("BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016","BIP 2016",
               "BIP 2016","BIP 2016","BIP 2016","BIP 2016", "BIP 2016","BIP 2016","BIP 2016","BIP 2016", "BIP 2016")
mBIP$VARIABLE <- c('dorsal/lateral MHb(0)', 'ventro/lateral MHb1(1)','vMHb1(2)', 'dMHb(3)', 'vMHb2(4)', 'early Hb1(5)','LHb subtype1(6)', 
                   'LHb subtype3(7)', 'early Hb2(8)', 'ventro/lateral MHb2(9)', 'lateral MHb(10)', 'superior MHb(11)', 'LHb subtype2(12)', 'HBX(13)',
                   'iHb1 (to MHb)(14)', 'iHb2 (to LHb)(15)', 'PC(16)')
mHeight <- read.table("data/output/MAGMA/Height/Oct2020_GWAS_Height_merged_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
mHeight$GWAS <- c("Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015","Height 2015",
                  "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015", "Height 2015","Height 2015","Height 2015","Height 2015","Height 2015")
mHeight$VARIABLE <- c('dorsal/lateral MHb(0)', 'ventro/lateral MHb1(1)','vMHb1(2)', 'dMHb(3)', 'vMHb2(4)', 'early Hb1(5)','LHb subtype1(6)', 
                      'LHb subtype3(7)', 'early Hb2(8)', 'ventro/lateral MHb2(9)', 'lateral MHb(10)', 'superior MHb(11)', 'LHb subtype2(12)', 'HBX(13)',
                      'iHb1 (to MHb)(14)', 'iHb2 (to LHb)(15)', 'PC(16)')
gsa_merge <- rbind(mMDD, mBMI, mSCZ, mBIP, mHeight)
rm(mMDD, mBMI, mSCZ, mBIP, mHeight)

#transform p-values to -log10 values
gsa_merge$logP <- -log10(gsa_merge$P)

#make plot
ggplot(gsa_merge, aes(y=logP, x=VARIABLE)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~GWAS, nrow=5, scales = "free_y") + geom_hline(yintercept=-log10(0.05/(17*5)), linetype="dashed") +
  xlab("") + ylab("-log10(P)") + theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
                                       axis.text.y = element_text(color="black", size=12),
                                       axis.title = element_text(size=14),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       strip.text = element_text(face="bold", size=14),
                                       legend.text = element_text(size=12, color="black"),
                                       legend.title = element_text(size=13, color="black")) + ggsave("figures/MAGMA/GWAS_merged-annotated_ANOVA.pdf")

#Developmental Hypothalamus dataset

#read in results and rbind() into melted format
hMDD <- read.table("data/output/MAGMA/MDD2-18/Nov2020_GWAS_MDD2-18_DevHypo_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hMDD$GWAS <- rep("MDD2 2018", 14)
hBMI <- read.table("data/output/MAGMA/BMI/Nov2020_GWAS_BMI_DevHypo_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hBMI$GWAS <- rep("BMI 2015", 14)
hSCZ <- read.table("data/output/MAGMA/SCZ/Nov2020_GWAS_SCZ_DevHypo_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hSCZ$GWAS <- rep("SCZ 2018", 14)
hBIP <- read.table("data/output/MAGMA/BIP/Nov2020_GWAS_BIP_DevHypo_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hBIP$GWAS <- rep("BIP 2016", 14)
hHeight <- read.table("data/output/MAGMA/Height/Nov2020_GWAS_Height_DevHypo_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hHeight$GWAS <- rep("Height 2015", 14)
gsa_hypo <- rbind(hMDD, hBMI, hSCZ, hBIP, hHeight)
rm(hMDD, hBMI, hSCZ, hBIP, hHeight)

#transform p-values to -log10 values
gsa_hypo$logP <- -log10(gsa_hypo$P)

#make plot
ggplot(gsa_hypo, aes(y=logP, x=VARIABLE)) + 
  geom_bar(position="dodge", stat="identity", fill="#999999") +
  facet_wrap(~GWAS, nrow=5, scales = "free_y") + geom_hline(yintercept=-log10(0.05/(14*5)), linetype="dashed") +
  xlab("") + ylab("-log10(P)") + theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
                                       axis.text.y = element_text(color="black", size=12),
                                       axis.title = element_text(size=14),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       strip.text = element_text(face="bold", size=14),
                                       legend.text = element_text(size=12, color="black"),
                                       legend.title = element_text(size=13, color="black")) + ggsave("figures/MAGMA/GWAS_DevHypo_ANOVA.pdf")

#La Manno 2016 dataset

#read in results and rbind() into melted format
hMDD <- read.table("data/output/MAGMA/MDD2-18/Nov2020_GWAS_MDD2-18_LaManno2016_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hMDD$GWAS <- rep("MDD2 2018", 31)
hBMI <- read.table("data/output/MAGMA/BMI/Nov2020_GWAS_BMI_LaManno2016_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hBMI$GWAS <- rep("BMI 2015", 31)
hSCZ <- read.table("data/output/MAGMA/SCZ/Nov2020_GWAS_SCZ_LaManno2016_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hSCZ$GWAS <- rep("SCZ 2018", 31)
hBIP <- read.table("data/output/MAGMA/BIP/Nov2020_GWAS_BIP_LaManno2016_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hBIP$GWAS <- rep("BIP 2016", 31)
hHeight <- read.table("data/output/MAGMA/Height/Nov2020_GWAS_Height_LaManno2016_ANOVA.gsa.out", header=TRUE, quote="\"")[,c(1,7)]
hHeight$GWAS <- rep("Height 2015", 31)
gsa_hypo <- rbind(hMDD, hBMI, hSCZ, hBIP, hHeight)
rm(hMDD, hBMI, hSCZ, hBIP, hHeight)

#transform p-values to -log10 values
gsa_hypo$logP <- -log10(gsa_hypo$P)

#make plot
ggplot(gsa_hypo, aes(y=logP, x=VARIABLE)) + 
  geom_bar(position="dodge", stat="identity", fill="#e06666") +
  facet_wrap(~GWAS, nrow=5, scales = "free_y") + geom_hline(yintercept=-log10(0.05/(31*5)), linetype="dashed") +
  xlab("") + ylab("-log10(P)") + theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
                                       axis.text.y = element_text(color="black", size=12),
                                       axis.title = element_text(size=14),
                                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                       strip.text = element_text(face="bold", size=14),
                                       legend.text = element_text(size=12, color="black"),
                                       legend.title = element_text(size=13, color="black")) + ggsave("figures/MAGMA/GWAS_LaManno16_ANOVA.pdf")


#create plots showing "Weak MAGMA signal" and "Strong MAGMA signal"
#define data
signal <- data.frame(c("c1", "c2", "c3", "c4", "c5"), c(2, 1.5, 0.5, 3, 1.2), c(7.8, 3.5, 1, 9.5, 2)) 
colnames(signal) <- c("cluster", "weak", "strong")

#plot weak signal
ggplot(signal, aes(y=weak, x=cluster)) + geom_bar(stat="identity", position="dodge") + geom_hline(yintercept=8, linetype="dashed") +
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14)) + xlab("") + ylab("") + scale_y_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10)) +
  ggsave("figures/MAGMA/MAGMAsignal_weak.pdf")

#plot strong signal
ggplot(signal, aes(y=strong, x=cluster)) + geom_bar(stat="identity", position="dodge") + geom_hline(yintercept=8, linetype="dashed") +
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14)) + xlab("") + ylab("") + scale_y_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10)) +
  ggsave("figures/MAGMA/MAGMAsignal_strong.pdf")
