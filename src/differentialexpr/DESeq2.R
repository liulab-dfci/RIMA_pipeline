rm(list=ls())
suppressMessages(library(DESeq2))
suppressMessages(library(optparse))
suppressMessages(library(tximport))
suppressMessages(library(dplyr))
suppressMessages(library(sva))
suppressMessages(library(limma))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="list input files", metavar="character"),
  make_option(c("-b", "--batch"), type="character", default=NULL,
              help="control batch effect or not", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="salmon", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=NULL,
              help="metasheet info", metavar="character"),
  make_option(c("-x", "--tx2gene"), type="character", default=NULL,
              help="transcript annotation", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./",
              help="output file path and prefix ", metavar="character"),
  make_option(c("-g", "--condition"), type="character", default="./",
              help="Condition to do comparison", metavar="character"),
  make_option(c("-r", "--treatment"), type="character", default="./",
              help="Treatment", metavar="character"),
  make_option(c("-c", "--control"), type="character", default="./",
              help="Control", metavar="character"),
  make_option(c("-p", "--pcoding"), type="character", default="./",
              help="proding coding gene list", metavar="character")
)
	      

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

batch <- opt$batch 
Condition <- opt$condition
metadata <- opt$meta
gene <- opt$tx2gene
Treatment <- opt$treatment
print(Treatment)
Control <- opt$control
print(Control)
Type <- opt$type
input <- opt$input

print("Reading meta file ...")
meta <- read.table(file = metadata, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
#print(head(meta))
#samples <- subset(meta, meta[,Condition] == Treatment | meta[,Condition] == Control)
samples <- subset(meta, meta[,Condition] != 'NA')
#print(samples)


print ("Reading tx2gene file ...")
tx2gene <- read.table(file = gene,sep = ",",header = TRUE)


if (is.null(opt$input) || is.null(opt$tx2gene) || is.null(opt$type)){
  print_help(opt_parser)
  stop("At least 3 arguments must be supplied ", call.=FALSE)
}

######################################--------------programme-------------##################################
###----If you have transcript quantification files, as produced by Salmon, Sailfish, or kallisto, you would use DESeqDataSetFromTximport.

Transcript <- function(files,samples,tx2gene,Type,batch){

  filelist <- strsplit(files, "\\,")[[1]]
  print(filelist)
  #print(rownames(meta))
  filelist.samples <- sapply(rownames(samples), function(x) grep(paste0("\\b",x,"\\b"), filelist, value = TRUE))
  
  #exactly match may output a list, need to convert list to character
  filelist.samples <- filelist.samples[lapply(filelist.samples,length)>0]
  print(paste("There are ",length(filelist.samples), "samples to be compared ...", sep = ""))

  tmp_chr <- as.character(filelist.samples)
  names(tmp_chr) <- names(filelist.samples)

  filelist.samples <- tmp_chr
  print(filelist.samples)

  txi <- tximport(filelist.samples, type=Type, tx2gene=tx2gene)
  txi$length[txi$length == 0] <- 1
  
  print(head(txi$counts))
  exprsn <- txi$counts
  exprsn_log <- log2(exprsn + 1)

  
  if(batch != "False"){
    print (paste("Running DESeq2 with batch effects ",opt$batch," on ",opt$condition,sep=""))  
     
    colData <- samples[,c(batch ,Condition)]
    colnames(colData) <- c("Batch","Condition")
       
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = colData,
                                       design = ~ Batch + Condition)
                                       
    #print (paste("Generating log transformed TPMS after batch correction of ",batch," on ",Condition,sep=""))
    
    #expr.limma = tryCatch(
    #                 removeBatchEffect(exprsn_log,colData$Batch),
    #                 error = function(e){
    #                 print(e)
    #                 })
    
    print ("Generating TPM matrix ...")
    write.table(exprsn,paste(opt$outpath,opt$condition,'_',opt$treatment,'_vs_',opt$control,'_TPMs.txt',sep = ""),quote = FALSE,sep = "\t")

  }else{
    print(paste("Running DESeq2 on ",opt$condition,sep=""))
    
    colData <- cbind(rownames(samples),samples[,Condition])
    colnames(colData) <- c('sample','Condition')
    print(colData)
    
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = colData,
                                       design = ~ Condition)
                                       
    print ("Generating TPM matrix ...")
    write.table(exprsn,paste(opt$outpath,opt$condition,'_',opt$treatment,'_vs_',opt$control,'_TPMs.txt',sep = ""),quote = FALSE,sep = "\t")
  }
  dds <- DESeq(ddsTxi)
  
  return (dds)
}


print ("Star running DESeq2 ...")
dds <- Transcript(files = input,samples, tx2gene = tx2gene, Type = Type, batch = batch)
print(class(dds))

save.image("git_version.RData")

print (paste("Comparing ",opt$treatment , " VS ", opt$control, sep = ""))
res <- results(dds, contrast = c("Condition",c(opt$treatment,opt$control)))
			     
source("src/differentialexpr/clusteringheatmap.R")
clustering_heatmap(dds, res, opt$pcoding, opt$outpath, opt$treatment, opt$control)

res_final <- as.data.frame(res)
res_final$Gene_name <- rownames(res_final)
res_final <- res_final[c(7, 1:6)]
res_final$`-log10(padj)` <- -log10(res_final$padj)
write.table(res_final,file = paste(opt$outpath,opt$condition,'_',opt$treatment,'_vs_',opt$control,'_DESeq2.txt',sep = ""), quote = FALSE,sep = "\t")
