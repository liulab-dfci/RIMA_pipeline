suppressMessages(library(ggplot2))
suppressMessages(library(ggfortify))
suppressMessages(library(RColorBrewer))
suppressMessages(library(optparse))
suppressMessages(library(ggpubr))

###set parameters
option_list = list(
  make_option(c("-b", "--before"), type="character", default=NULL, 
              help="expression matrix before batch removal", metavar="character"),
  make_option(c("-a", "--after"), type="character", default=NULL, 
              help="expression matrix after batch removal", metavar="character"),
  make_option(c("-m", "--metainfo"), type="character", default=NULL, 
              help="meta file", metavar="character"),
  make_option(c("-i", "--out_before"), type="character", default=NULL, 
              help="PCA plot before batch removal", metavar="character"),
  make_option(c("-j", "--out_after"), type="character", default=NULL,
              help="PCA plot after batch removal", metavar="character"),
  make_option(c("-c", "--column"), type="character", default=NULL, 
              help="batch column", metavar="character"),
  make_option(c("-g", "--group"), type="character", default=NULL,
              help="group column", metavar="character")  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

expr_before <- read.table(opt$before,sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
expr_after <- read.table(opt$after,sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
annot <- read.table(opt$metainfo,sep = ",",header = TRUE, row.names = 1)


if (opt$column == "False") {
  annot$Batch <- 1
  Batch <- "Batch"
} else {
  Batch <- opt$column
}

print ('Load data done!')

####function of pca plotting
pca_plot <- function(exprTable, annot,title, Batch) {
  batch_n <- length(unique(as.character(annot[colnames(exprTable),Batch])))
  print(paste("there are ", batch_n, " batches in your data"))
  df <- cbind.data.frame(t(exprTable),batch = as.character(annot[colnames(exprTable),Batch]))
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


###print pca plot

if (opt$column != "False") {
	pdf(file = opt$out_before, height = 5, width = 6)
	p_before <- pca_plot(exprTable = expr_before, annot = annot, title = "PCA plot Before Batch Removal",Batch=Batch)
        print(p_before)
	dev.off()

	pdf(file = opt$out_after, height = 5, width = 6)
	p_after <- pca_plot(exprTable = expr_after, annot = annot, title = "PCA plot After Batch Removal",Batch = Batch)
        print(p_after)
	dev.off()
} else
{
	pdf(file = opt$out_before, height = 5, width = 6)
        plot.new()
	dev.off()

	pdf(file = opt$out_after, height = 5, width = 6)
	plot.new()
	dev.off()
}
###
