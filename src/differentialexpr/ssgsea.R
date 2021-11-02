#!/usr/bin/env Rscript

#dependencies
suppressMessages(library(GSEABase))
suppressMessages(library(GSVA))
suppressMessages(library(RColorBrewer))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(gplots))
suppressMessages(library(textshape))
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
expr.dat <- read.table(exprsn, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

meta <- read.table(file = metadata, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

#add tmp column if the metasheet only has two columns
if(ncol(meta) == 1) {
meta$tmp <- 1
}

samples <- subset(meta,meta[,Condition] == Treatment | meta[,Condition] == Control)
samples <- samples[order(samples[[Condition]], decreasing = TRUE),]

#####read in data
ssgsea <- gsva(as.matrix(expr.dat), geneSets,
               method="ssgsea",
               ssgsea.norm=TRUE,
               verbose=TRUE)
            

               
###min-max normalization(make scale from -1 to 1)
ssgsea.nor <- t(apply(ssgsea,1,function(x) (2*(x-min(x))/(max(x)-min(x)))-1))


###get the difference between Clinic traits for each term
  
annot.col <- data.frame(Clinic = samples[,Condition],row.names = rownames(samples))
gr.samples <- lapply(unique(annot.col$Clinic), function(x) return(rownames(annot.col)[which(annot.col$Clinic == x)]))
  



###calculate between-group variation for filtering terms with top variance
SelTerm <- lapply(c(1:dim(ssgsea.nor)[1]), function(t) {
    grandmean <- mean(ssgsea.nor[t,])
    bgV <- lapply(c(1:length(gr.samples)), function(g){
      bg <- ((ssgsea.nor[t,gr.samples[[g]]] - grandmean)^2)*length(gr.samples[[g]])
      return (bg)
    })
    bgsum <- sum(unlist(bgV))
    return (bgsum)
  })
  
  
names(SelTerm) <- rownames(ssgsea.nor)
gr.diff.sort <- sort(unlist(SelTerm), decreasing = TRUE)
ssgsea.sort <- ssgsea[names(gr.diff.sort),]


        
write.table(ssgsea.sort, paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_ssgsea.txt", sep = ""), 
            quote = FALSE, sep = "\t", row.names = FALSE)   

png(paste(Outdir,Condition, "_",Treatment,"_vs_",Control,"_ssgsea.png", sep = ""),width = 1000, height = 880)
col_fun = colorRamp2(c(-1,0,1), c("blue","white", "red"))
print(
Heatmap(ssgsea.sort[1:top_n,] ,cluster_rows = TRUE,cluster_columns = FALSE,show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),
        col = col_fun, name='ssgsea')
        )
dev.off()
    




