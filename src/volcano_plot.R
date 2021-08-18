#!/usr/bin/env Rscript

#dependencies
library(ggplot2)
library(dplyr)
library(ggrepel)
library(optparse)

option_list = list(
  make_option(c("-m", "--deseq2_mat"), type="character", default=NULL, 
              help="signature reference", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####read in data
# data <- read.table("~/Documents/rnaseq/data/DESeq2_ConvertID.txt", sep = "\t", header = TRUE, row.names = 1)
data <- read.table(opt$deseq2_mat, sep = "\t", header = TRUE, row.names = 1)
###set fdr cutoff as 0.05, log2FC as 1
xl <- -1
xr <- 1
yp <- 1.3
###label differentially expressed genes
dif.data<- na.omit(data) %>% 
  mutate(logP = -log10(padj)) %>% 
  mutate(color = ifelse(log2FoldChange > xr & logP > yp, 
                        yes = "Treat", no = ifelse(log2FoldChange < xl & logP > yp, yes = "Control",  no = "none")))
###extract top differentially expressed genes
top_n <- 10
up <- arrange(subset(dif.data, color == "Treat"),desc(log2FoldChange)) 
down <- arrange(subset(dif.data, color == "Control"),log2FoldChange)
if(nrow(up) >= 10 && nrow(down) >= 10){
  top_labelled <- rbind.data.frame(up[1:top_n,], down[1:top_n,])
}else if (nrow(up) < 10 && nrow(down) >= 10){
  top_labelled <- rbind.data.frame(up, down[1:top_n,])
}else if (nrow(up) >= 10 && nrow(down) < 10){
  top_labelled <- rbind.data.frame(up[1:top_n,], down)
}else{
  top_labelled <- rbind.data.frame(up, down)
}
###volcano plot
#pdf(paste(opt$outdir,"diff_volcano_plot.pdf", sep = ""),width = 5, height = 5)
#png(paste("~/Documents/rnaseq/data/","diff_volcano_plot.png", sep = ""), res = 300, width = 1800, height = 1800)
png(paste(opt$outdir), res = 300, width = 1900, height = 1000)

ggplot(dif.data, aes(x = log2FoldChange, y = logP)) + 
  geom_point(aes(color = factor(color)), size = 1.55, alpha = 0.8, na.rm = TRUE) + # add gene points
  theme_bw(base_size = 10) + # clean up theme
  theme(legend.position = "none") + # remove legend 
  xlab(expression(log[2]("FC"))) + # x-axis label
  ylab(expression(-log[10]("adj.P.Val"))) + # y-axis label
  theme(axis.title = element_text(size = 8)) + 
  theme(axis.text = element_text(size = 8)) +
  geom_vline(xintercept = xl, colour = "black", linetype = "dashed") + # add line at 0
  geom_vline(xintercept = xr, colour = "black", linetype = "dashed") + # add line at 0
  geom_hline(yintercept = yp, colour = "black", linetype = "dashed") + # p(0.05) = 1.3
  annotate(geom = "text", 
           label = "Control", 
           x = min(down$log2FoldChange)*0.7, y = max(dif.data$logP)*0.9, 
           size = 3, colour = "#3182bd") + # add Down text
  annotate(geom = "text", 
           label = "Treat", 
           x = max(up$log2FoldChange)*0.7, y = max(dif.data$logP)*0.9, 
           size = 3, colour = "#E64B35") + # add Up text
  scale_color_manual(values = c("Treat" = "#E64B35", 
                                "Control" = "#3182bd", 
                                "none" = "#636363")) + # change colors
  scale_y_continuous(trans = "log1p")+  # Scaled Y-axis with log1p function
  geom_label_repel(data = top_labelled, 
                   aes(label = gene_name), fill = "white", size = 2.0,
                   fontface = 'bold',box.padding = 0.2, color = 'black',
                   label.size = 0.10,point.padding = 0.5,segment.color = 'gold')


dev.off()

