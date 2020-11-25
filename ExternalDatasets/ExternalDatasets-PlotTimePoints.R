##### plotting time point distribution in each cluster for the Kim et al 2020 and LaManno et al 2016 datasets #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(ggplot2)

#create Kim et al 2020 dataframe
hypo_ages <- data.frame(c(rep("Astrocytes",4), rep("Endothelial cells", 12), rep("Ependymal cells", 4), rep("Hypothalamic neurons", 6), rep("Immune cells", 8),
                     rep("LH", 4), rep("Microglia", 4), rep("Non-neuronal cells", 3), rep("Oligodendrocytes", 4), rep("PMN", 5), rep("POA", 2), 
                     rep("Tanycytes", 4), rep("Ventral hypothalamic NPC", 2), rep("VMH", 5)), 
                          c("P4", "P8", "P14", "P45",
                            "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45",
                            "P4", "P8", "P14", "P45",
                            "E16", "E18", "P4", "P8", "P14", "P45",
                            "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18",
                            "E12", "E13", "E14", "E15",
                            "P4", "P8", "P14", "P45",
                            "E16", "E18", "P4",
                            "P4", "P8", "P14", "P45",
                            "E11", "E12", "E13", "E14", "E15",
                            "E13", "E14",
                            "P4", "P8", "P14", "P45",
                            "E11", "E15",
                            "E11", "E12", "E13", "E14", "E15")) 
colnames(hypo_ages) <- c("cluster", "timepoint")
positions <- c("E10", "E11", "E12", "E13", "E14", "E15", "E16", "E18", "P4", "P8", "P14", "P45")

#plot
ggplot(hypo_ages, aes(y=timepoint, x=cluster)) + geom_point(stat="identity", position="dodge", color="#999999", size=3) +
  theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black")) + xlab("") + ylab("") + scale_y_discrete(limits = rev(positions)) + #use rev() when E10 at top
  ggsave("figures/ExternalDatasets/DevHypo_Sampling-E10.pdf")

#create LaManno et al 2016 dataframe
lamanno_ages <- data.frame(c("DA.SNC", "DA.VTA1", "DA.VTA2", "DA.VTA3", "DA.VTA4", rep("mDA0", 5), rep("mDA1", 5), rep("mDA2", 4), rep("mEndo", 6), 
                             rep("mEpen", 3), rep("mGaba1a", 5), rep("mGaba1b", 5), rep("mGaba2", 5), rep("mMgl", 5), rep("mNbDA", 5),
                             rep("mNbL1", 3), rep("mNbL2", 4), rep("mNbM", 3), rep("mNbML1", 4), rep("mNbML2", 5),rep("mNbML3", 3),rep("mNbML4", 4),
                             rep("mNbML5", 6),rep("mProg", 6), rep("mOMTN", 6),rep("mPeric", 6),rep("mRgl1", 6),rep("mRgl2", 5),rep("mRgl3", 5),
                             rep("mRN", 5), rep("mSert", 5)), 
                        c(rep("adult", 5),
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E14.5", "E15.5", "E18.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5",
                          "E11.5", "E12.5", "E13.5",
                          "E11.5", "E12.5", "E13.5", "E14.5",
                          "E11.5", "E12.5", "E13.5",
                          "E11.5", "E12.5", "E13.5", "E14.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5",
                          "E11.5", "E12.5", "E13.5",
                          "E11.5", "E12.5", "E13.5", "E14.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5",
                          "E11.5", "E12.5", "E13.5", "E14.5", "E15.5",
                          "E12.5", "E13.5", "E14.5", "E15.5", "E18.5")) 
colnames(lamanno_ages) <- c("cluster", "timepoint")
positions <- c("E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E18.5", "adult")

#plot
ggplot(lamanno_ages, aes(y=timepoint, x=cluster)) + geom_point(stat="identity", position="dodge", color="#e06666", size=3) +
  theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black")) + xlab("") + ylab("") + scale_y_discrete(limits = positions) + #use rev() when E10 at top
  ggsave("figures/ExternalDatasets/LaManno_Sampling-adult.pdf")
