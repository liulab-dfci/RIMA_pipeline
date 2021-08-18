#!/usr/bin/env Rscript
###library packages
rm(list=ls())
library(DESeq2)
library(optparse)
library(tximport)
library(dplyr)

###define options
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list input files", metavar="character"),
  make_option(c("-b", "--batch"), type="character", default=NULL, 
              help="control batch effect or not", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="data source", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL, 
              help="comparison design matrix", metavar="character"),
  make_option(c("-c", "--comparison"), type="character", default=NULL, 
              help="comparison type", metavar="character"),
  make_option(c("-x", "--tx2gene"), type="character", default=NULL, 
              help="transcript annotation", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$tx2gene) || is.null(opt$type)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied ", call.=FALSE)
}

######################################--------------programme-------------##################################
###----If you have transcript quantification files, as produced by Salmon, Sailfish, or kallisto, you would use DESeqDataSetFromTximport.
Transcript <- function(files,samples,tx2gene,Type,batch){
  # filelist <- dir(path = "~/Documents/rnaseq/data", pattern = "genes.results$", full.names = TRUE)
  # tx2gene <- read.table(file = "~/Documents/rnaseq/data/Ens2gene.csv",sep = ",",header = TRUE)#~/Documents/rnaseq/data/Ens2gene.csv
  # Type <- "rsem"
  # samples <- read.table(file = "~/Documents/rnaseq/snakemake_files/metasheet_test.csv", sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  # samples$Batch <- factor(samples$Batch)
  # samples$Tissue <- factor(samples$Tissue)
  filelist <- strsplit(files, "\\,")[[1]]
  filelist.samples <- sapply(rownames(samples), function(x) grep(x, filelist, value = TRUE))
  filelist.samples <- filelist.samples[lapply(filelist.samples,length)>0]
  print(filelist.samples)
  txi <- tximport(filelist.samples, type=Type, tx2gene=tx2gene)
  txi$length[txi$length == 0] <- 1
  if(batch == "yes"){
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Batch + Condition)
  }
  else{
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = samples,
                                       design = ~ Condition)
  }
  dds <- DESeq(ddsTxi)
  return (dds)
}

###----If you have count matrix files, the first line would use DESeqDataSetFromMatrix.
CountMatrix <- function(files,samples,batch){
  cts <- read.table(files, sep = ",", header=TRUE, row.names = 1)
  cts <- cts[, rownames(samples)]
  if(batch == "yes"){
    ddsMat <- DESeqDataSetFromMatrix(countData = cts,
                                     colData = samples,
                                     design = ~ Condition + Batch)
  }
  else{
    ddsMat <- DESeqDataSetFromMatrix(countData = cts,
                                     colData = samples,
                                     design = ~ Condition)
  }
  dds <- DESeq(ddsMat)
  return (dds)
}

###----ID conversion
IDConversion <- function(Type, res, tx2gene){
  if(Type == "salmon"){
    res.convertID <- cbind.data.frame(gene_name =rownames(res),res)
  }
  if(Type != "salmon") {
    sub.ref <- tx2gene[,c(1:2)] ##only extract gene symbol and ensembl
    sub.ref <- sub.ref[!duplicated(sub.ref),]
    res.convertID <- cbind.data.frame(gene_name = sub.ref[match(rownames(res), sub.ref$ID),"GENENAME"],res)
  }
  return (res.convertID)
}


###differential gene expression
###input gene mapping file
tx2gene <- read.table(file = opt$tx2gene,sep = "\t",header = TRUE)
###comparison type
comp.type <- opt$comparison
###input comparion design file
comp <- opt$design
# name <- strsplit(comp,"\\/")[[1]][length(strsplit(comp,"\\/")[[1]])]
samples <- read.table(file = comp, sep='\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
samples.all <- subset(samples, Batch != "" & Condition != "")
samples.all$Batch <- factor(samples.all$Batch)
print(paste("There are ", dim(samples.all)[1], " samples to be compared", sep = ""))
###run deseq2
####compare between any two groups
if(comp.type == "between"){
  samples.all$Condition <- factor(samples.all$Condition)
  ss <- t(as.matrix(combn(levels(samples.all$Condition),2)))
  if(opt$type == "salmon" || opt$type == "rsem"){  
    dds <- Transcript(files = opt$input, samples = samples.all, tx2gene = tx2gene, Type = opt$type, batch = opt$batch)  
    print(class(dds))
  }
  if(opt$type == "count.matrix"){
    dds <- CountMatrix(files = opt$input, samples = samples.all, batch = opt$batch)
  }
  for(i in 1:dim(ss)[1]){
    print (paste("Compare ", ss[i,1], " and ", ss[i,2], sep = ""))
    res <- results(dds, contrast = c("Condition",c(ss[i,1],ss[i,2])))
    res.convertID <- IDConversion(opt$type, res, tx2gene)
    write.table(res,file = paste(opt$outpath,"_",ss[i,1],"_VS_",ss[i,2],'_DESeq2_raw.txt',sep = ""), quote = FALSE,sep = "\t")
    write.table(res.convertID,file = paste(opt$outpath,"_",ss[i,1],"_VS_",ss[i,2],'_DESeq2_ConvertID.txt',sep = ""), quote = FALSE,sep = "\t")
  } 
}
####compare one group with all the others
if(comp.type == "loop"){
  for(gr in unique(as.character(samples.all$Condition))){
    samples.loop <- samples.all %>% mutate(Condition_loop = factor(ifelse(Condition == gr, gr, "others"))) %>% mutate(Condition = Condition_loop)
    print(samples.loop$Condition)
    if(opt$type == "salmon" || opt$type == "rsem"){  
      dds <- Transcript(files = opt$input, samples = samples.loop, tx2gene = tx2gene, Type = opt$type, batch = opt$batch)  
    }
    if(opt$type == "count.matrix"){
      dds <- CountMatrix(files = opt$input, samples = samples.loop, batch = opt$batch)
    }
    print (paste("Compare ", gr, " and others", sep = ""))
    res <- results(dds, contrast = c("Condition",c(gr,"others")))
    res.convertID <- IDConversion(opt$type, res, tx2gene)
    write.table(res,file = paste(opt$outpath,"_",gr,"_VS_others_DESeq2_raw.txt",sep = ""), quote = FALSE,sep = "\t")
    write.table(res.convertID,file = paste(opt$outpath,"_",gr,"_VS_others_DESeq2_ConvertID.txt",sep = ""), quote = FALSE,sep = "\t")
  }
}





