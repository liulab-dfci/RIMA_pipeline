#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Align report matrix", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data <- read.csv(file = opt$input, sep=",", header=TRUE, check.names=FALSE, row.names = 1)


###function of showing aligned reads fraction
map_stat <- function(data){
  ###get total reads and aligned reads
  x <- data.frame( Sample=names(data), Total_Reads=as.numeric(as.matrix(data["Number_of_input_reads",])), Unique_Reads=as.numeric(as.matrix(data["Uniquely_mapped_reads_number",])))
  x1 <- melt(x, id.var="Sample")
  ###modify scale of y axis
  upper_limit <- max(x$Total_Reads)
  limits <- seq( 0, upper_limit, length.out=10)
  cust_labels <- vector("character",length=length(limits))
  if( nchar(upper_limit) < 7 ) {
    cust_labels <- paste(round(limits/1000),"K",sep="") 
    limits <- round(limits/1000) * 1000
  } else {
    cust_labels <- paste(round(limits/1000000),"M",sep="") 
    limits <- round(limits/1000000) * 1000000
  }
  ###plot
  colors <- c(Total_Reads="Grey", Unique_Reads="steelblue")
  gp <- ggplot(x1, aes(x=Sample, y=value, fill=variable)) + 
    geom_bar( stat = "identity", position="identity") + 
    scale_y_continuous("",limits=c(0,upper_limit), labels=cust_labels, breaks=limits) + 
    scale_fill_manual(values=colors) + 
    labs( title="Read Alignment Report", x = "Sample Names", y="") + 
    guides(fill=guide_legend(title=NULL)) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
            axis.text.x=element_text(angle = -90,size=10,face = "bold",hjust=1),
            axis.text.y=element_text(size=10,face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "top")
  return(gp)
}

###main
Pwidth=60*ncol(data)+200
png(file = opt$output, res = 300, width = Pwidth, height = 1500)
map_stat(data = data)
dev.off()

