#!/usr/bin/env Rscript

library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(optparse)

###set parameters
option_list = list(
  make_option(c("-b", "--before"), type="character", default=NULL, 
              help="expression matrix before batch removal", metavar="character"),
  make_option(c("-a", "--after"), type="character", default=NULL, 
              help="expression matrix after batch removal", metavar="character"),
  make_option(c("-m", "--metainfo"), type="character", default=NULL, 
              help="meta file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

expr_before <- read.table(opt$before,sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
expr_after <- read.table(opt$after,sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
annot <- read.table(opt$metainfo,sep = ",",header = TRUE, row.names = 1)
pca_out_dir <- opt$outdir

####function of pca plotting
pca_plot <- function(exprTable, annot,title) {
  batch_n <- length(unique(as.character(annot[colnames(exprTable),"Batch"])))
  print(paste("there are ", batch_n, " batches in your data"))
  df <- cbind.data.frame(t(exprTable),batch = as.character(annot[colnames(exprTable),"Batch"]))
  pca_plot <- autoplot(prcomp(t(exprTable)), data = df, col = 'batch', size = 1, frame = TRUE, frame.type = 'norm')+ 
    labs(title=title)+ 
    scale_color_manual(values = brewer.pal(name = "Set1", n = 9)[1:batch_n])+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"))
  return(pca_plot)
}

# expr_after <- read.table(file = "~/Documents/pipeline_test/tongji_test/rnaseq_pipeline/tpm_convertID.batch", sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
# expr_before <- read.table(file = "~/Documents/rnaseq/data/rsem_tpm_gene_convertID.txt", sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
# annot <- read.table(file = "~/Documents/pipeline_test/tongji_test/rnaseq_pipeline/metasheet.csv", sep = ",",header = TRUE, row.names = 1)
# pca_out_dir <- "~/Documents/rnaseq/data/"
###print pca plot
png(file = paste(pca_out_dir,"pca_plot_before.png",sep = ""), res = 300, height = 1200, width = 1500)
pca_plot(exprTable = expr_before, annot = annot, title = "PCA plot Before Batch Removal")
dev.off()

png(file = paste(pca_out_dir,"pca_plot_after.png",sep = ""), res = 300, height = 1200, width = 1500)
pca_plot(exprTable = expr_after, annot = annot, title = "PCA plot After Batch Removal")
dev.off()

