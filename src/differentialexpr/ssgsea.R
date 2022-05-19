#!/usr/bin/env Rscript

#dependencies
suppressMessages(library(GSEABase))
suppressMessages(library(GSVA))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-e", "--expression"), type="character", default=NULL,
              help="signature reference", metavar="character"),
              
  make_option(c("-f", "--gmt"), type="character", default=NULL,
              help="pathway database gmt file", metavar="character"),
              
  make_option(c("-t", "--treatment"), type="character", default=NULL,
              help="treatment", metavar="character"),
              
  make_option(c("-c", "--control"), type="character", default=NULL,
              help="control", metavar="character"),
              
  make_option(c("-g", "--condition"), type="character",
              help="column name in the meta to compare"),
              
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="meta info", metavar="character"),
              
  make_option(c("-n", "--top"), type="character", default=NULL,
              help="top n terms for displaying", metavar="character"),

  make_option(c("-a", "--hallmark"), type="character", default=NULL,
              help="hallmark gmt file", metavar="character"),
              
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character")
);



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Outdir <- opt$outdir
gmt_file <- opt$gmt
Condition <- opt$condition
Treatment <- opt$treatment
Control <- opt$control
exprsn <- opt$expression
metadata <- opt$meta
top_n <- opt$top

###parameters
geneSets <- getGmt(gmt_file)
hallmark <- getGmt(opt$hallmark)
expr.dat <- read.table(exprsn, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

meta <- read.table(file = metadata, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

#add tmp column if the metasheet only has two columns
if(ncol(meta) == 1) {
meta$tmp <- 1
}

samples <- subset(meta,meta[,Condition] == Treatment | meta[,Condition] == Control)
samples <- samples[order(samples[[Condition]], decreasing = TRUE),]

#####read in data
print("Running the ssgsea using kegg pathways")
ssgsea_kegg <- gsva(as.matrix(expr.dat), geneSets,
               method="ssgsea",
               ssgsea.norm=TRUE,
               verbose=TRUE)
            

print("Running the ssgsea using hallmark pathways")
ssgsea_hallmark <- gsva(as.matrix(expr.dat), hallmark,
               method="ssgsea",
               ssgsea.norm=TRUE,
               verbose=TRUE)

###look at the most different terms
ssgsea_plot <- function(samples, Condition, Treatment, Control, ssgsea) {
	
  ssgsea <- as.data.frame(ssgsea)
  tx = samples[samples[,Condition] == Treatment,][1]
  ssgsea_tx <- ssgsea[as.character(tx[[1]])]

  cn <- samples[samples[,Condition] == Control,][1]
  ssgsea_cn <- ssgsea[as.character(cn[[1]])]

  ta <- NULL
  for (t in rownames(ssgsea_tx)) {
    wilcox <- wilcox.test(as.numeric(ssgsea_tx[t,]), as.numeric(ssgsea_cn[t,]), alternative = "two.sided")
    tmp_ta <- data.frame(terms = t, pvalue = wilcox$p.value)
    ta <- rbind(ta, tmp_ta)
  }

  ta <- ta[order(ta$pvalue, decreasing = FALSE),]
  ssgsea.sort <- ssgsea[ta$terms,]
  ssgsea.sort <- ssgsea.sort[as.character(samples$SampleName)]
  
  ssgsea.sort <- as.matrix(ssgsea.sort)
  samples <- samples[,c("SampleName",Condition)]
  rownames(samples) <- samples[[1]]
  colnames(samples)[2] <- "Condition"

  melted_cormat <- melt(ssgsea.sort[1:30,], na.rm = TRUE)

  melted_cormat <- merge(melted_cormat, samples, by.x = 2, by.y = 1)
  melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = rev(unique(as.character(ta$terms[1:30]))))

  p <- ggplot(data = melted_cormat, aes(Var2, as.factor(Var1), fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(min(melted_cormat$value),max(melted_cormat$value)), space = "Lab",
                         name="ssgsea score") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                     hjust = 1),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 14)) + facet_grid(.~Condition, scales = "free_x", space='free') +
    labs(x = "", y = "")

  return(list(p,ssgsea.sort))
}

ssgsea_kegg <- ssgsea_plot(samples, Condition, Treatment, Control, ssgsea_kegg)

ssgsea_hallmark <- ssgsea_plot(samples, Condition, Treatment, Control, ssgsea_hallmark)

        
write.table(ssgsea_kegg[[2]], paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_ssgsea.txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)  

write.table(ssgsea_hallmark[[2]], paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_hallmark_ssgsea.txt", sep = ""),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)    




