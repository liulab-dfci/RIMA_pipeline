#!/usr/bin/env Rscript

#dependencies
library(ggplot2)
library(optparse)

option_list = list(
  make_option(c("-t", "--tin_score"), type="character", default=NULL, 
              help="tin summary score", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##read in data
tin <- read.table(opt$tin_score,sep = "\t", header = TRUE)
#tin <- read.table("~/Desktop/plot/qc/tin_score_summary.txt", sep = "\t", header = TRUE)
tin$Bam_file <- gsub("(_downsampling.bam)","",tin$Bam_file)


##set y axis
upper_limit <- max(tin$TIN.median.)
limits <- seq( 0, upper_limit, length.out=10)
cust_labels <- round(limits)


##bar plot showing median TIN score
Pwidth=(60*nrow(tin))+200
png(paste(opt$outdir,"medTIN_score_plot.png", sep = ""), res = 300, width = Pwidth, height = 1500)
#png(paste("~/Desktop/plot/qc/","medTIN_score_plot.png", sep = ""), res = 300, width = Pwidth, height = 1500)
ggplot(tin, aes(x=Bam_file, y=TIN.median.)) + 
  geom_bar(stat = "identity", position="identity", fill = "steelblue") + 
  scale_y_continuous("", limits=c(0,upper_limit), labels=cust_labels, breaks=limits) + 
  labs(title="median TIN score", x = "Sample Names", y="") + 
  guides(fill=guide_legend(title=NULL)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(angle = -90, size=10,face = "bold",hjust=1),
        axis.text.y=element_text(size=10,face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))
dev.off()
