##### plotting the number of genes used in each MAGMA analysis #####
# author: Juliska E Boer
# date: 25 Nov 2020

#load packages
setwd("E:/")
library(ggplot2)

#create dataframe with all (sampled) datasets
genes_used <- data.frame(c("Developmental Hb", "Developmental Hb - sampled", "Hashikawa, et al (2020)", "Wallace, et al (2020)", "Tabula Muris, et al (2018)",
                        "Tabula Muris, et al (2018) - sampled", "Kim, et al (2020) - sampled", "LaManno, et al (2016)"), 
                      c(13330, 11304, 11392, 11216, 14033, 11304, 11381, 12214)) 
colnames(genes_used) <- c("dataset", "genes")

#plot
ggplot(genes_used, aes(y=genes, x=dataset, fill=dataset)) + geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black")) + xlab("") + ylab("") + coord_cartesian(ylim=c(10000, 14500)) +
  scale_fill_manual(values = c("#66c2a5", "#66c2a5", "#7690ca", "#999999", "#e06666", "#e78ac3", "#e78ac3", "#fc8d62")) +
  ggsave("figures/MAGMA/GenesInMAGMA.pdf")


#create dataframe with all (sampled datasets) except Kim, et al (2020) and La Manno, et al (2016)
genes_used <- data.frame(c("Developmental Hb", "Developmental Hb - sampled", "Hashikawa, et al (2020)", "Wallace, et al (2020)", "Tabula Muris, et al (2018)",
                           "Tabula Muris, et al (2018) - sampled"), 
                         c(13330, 11304, 11392, 11216, 14033, 11304)) 
colnames(genes_used) <- c("dataset", "genes")

ggplot(genes_used, aes(y=genes, x=dataset, fill=dataset)) + geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90, hjust=0.95,vjust=0.2, color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.text = element_text(face="bold", size=14),
        legend.text = element_text(size=12, color="black"),
        legend.title = element_text(size=13, color="black")) + xlab("") + ylab("") + coord_cartesian(ylim=c(10000, 14500)) +
  scale_fill_manual(values = c("#66c2a5", "#66c2a5", "#7690ca", "#e78ac3", "#e78ac3", "#fc8d62")) +
  ggsave("figures/MAGMA/GenesInMAGMA-noKimManno.pdf")
