##### plotting time point distribution in each cluster for the Kim et al 2020 and LaManno et al 2016 datasets #####
# author: Juliska E Boer
# date: 13 Nov 2020

#load packages
setwd("E:/")
library(ggplot2)

#create Dev Hb dataframe
hb_ages <- data.frame(c(rep("early Hb1[0]",8), rep("early Hb2[1]", 5), rep("vMHb1[2]", 6), rep("vMHb2[3]", 6), rep("early dMHb[4]", 5),
                          rep("dMHb2[5]", 5), rep("iHb1[6]", 7), rep("dMHb3[7]", 4), rep("iHb2[8]", 7), rep("vMHb3[9]", 4), rep("dMHb1[10]", 5), 
                          rep("LHb[11]", 8), rep("PC2[12]", 8), rep("PC1[13]", 6)), 
                        c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18",
                          "E13", "E15", "E18", "P4", "P7", "adult",
                          "E13", "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18",
                          "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18", "P7", "adult",
                          "E15", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18", "P4", "P7",
                          "E18", "P4", "P7", "adult",
                          "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult",
                          "E11", "E12", "E13", "E18", "P4", "adult")) 
colnames(hb_ages) <- c("cluster", "timepoint")
positions <- c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult")

#plot
ggplot(hb_ages, aes(y=timepoint, x=cluster)) + geom_point(stat="identity", position="dodge", color="#66c2a5", size=3) +
  theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black")) + xlab("") + ylab("") + scale_y_discrete(limits = positions) + #use rev() when E10 at top
  ggsave("figures/DevHb_Sampling-adult.pdf")
