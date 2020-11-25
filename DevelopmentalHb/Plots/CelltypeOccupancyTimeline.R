##### plot cell type occupancy over sampling timeline #####
# author: Juliska E Boer
# date: 03 Nov 2020

#load packages
library(ggplot2)
library(Seurat)
library(reshape2)
library(dplyr)
setwd("E:/")

#read in Seurat object holding developmental Hb dataset
load("data/output/merge_adult/Embryo_Scanpy_Seurat_obj.RData")

#format data for ggplot (melted format)
data <- data.frame(row.names = row.names(final.embryo@meta.data),
                       stage = final.embryo$stage,
                       cluster = final.embryo$louvain)
data$cluster <- as.character(data$cluster)
data$day <- as.numeric(as.character(plyr::mapvalues(data$stage, from=c("E11", "E12", "E13", "E15", "E18", "P4", "P7", "adult"),
                                                    to=c(11,12,13,15,18,23,26,29))))
#calculate the percentage
data_devhb <- data  %>%
  group_by(day, cluster) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n))

#assign each cell type a color, corresponding to colors given by Louvain clustering (in SCANPY)
col_devhb <- c("0" = "#1f77b4",
                 "1" = "#ff7f0e",
                 "2" = "#279e68",
                 "3" = "#d62728",
                 "4" = "#aa40fc",
                 "5" = "#8c564b",
                 "6" = "#e377c2",
                 "7" = "#b5bd61",
                 "8" = "#17becf",
                 "9" = "#aec7e8",
                 "10" = "#ffbb78",
                 "11" = "#98df8a",
                 "12" = "#ff9896",
                 "13" = "#c5b0d5")

#order cell types
#this order was found with trial and error
data_devhb$cluster <- factor(data_devhb$cluster, levels = rev(c("13", "8", "0", "1", "12", "2", "3", "11","7", "6", "4", "10", "5",  "9")))

#plot cell type occupancy over sampling timeline
ggplot(data_devhb, aes(x=day, y=percentage, fill=cluster)) + 
  geom_area(alpha=1 , size=.1, colour="black") + 
  geom_vline(xintercept=c(11,12,13,15,18,23,26,29), linetype="dashed", size=0.3) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = col_devhb) +
  scale_x_continuous(limits=c(11, 29),breaks=c(seq(11, 29,1))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) + ylab("Percentage") + xlab("") +
  ggsave("figures/embryo_Hb/DevHb_celltype_occupancy_tp.pdf", height=4, width=12)
