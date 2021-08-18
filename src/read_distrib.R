# vim: syntax=r tabstop=4 expandtab

#--------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 1, 2015
#---------------------------
library(data.table)
library( ggplot2 )
library( reshape2 )

args <- commandArgs( trailingOnly = TRUE )

data <- fread( args[1], header=TRUE, check.names=FALSE)
# data <- fread("~/Downloads/read_distrib.matrix.tab", header = TRUE, check.names = FALSE)
data <- as.data.frame(data)
rownames(data) <- data$Feature
sub_data <- data[ c("Introns", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons"), ]

pdf( args[2], width = 8, height = 8)
suppressMessages(m_data <- melt( sub_data ))

ggplot(data=m_data, aes(variable, value, fill=Feature)) + 
  geom_bar(stat="identity") + 
  xlab("Sample") + 
  ylab("Feature") + 
  ggtitle("Reads Distribution") +
  theme_bw() + 
  coord_flip() +
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
        axis.text.y=element_text(size=12,face = "bold",hjust=1),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")
        )

