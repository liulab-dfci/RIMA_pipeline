##Script to subset metasheet into separate "Pre" and "Post" category for running TIDE.
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-m", "--meta"), type="character", default=NULL, 
              help="meta information", metavar="character"),
  make_option(c("-e", "--expression"), type="character",
              help="batch tpm expression data"),
  make_option(c("-o1", "--outdir1"), type="character", default=NULL, 
              help="output directory for pre samples", metavar="character"),
  make_option(c("-o2", "--outdir2"), type="character", default=NULL,
              help="output directory for post samples", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


meta <- read.table(opt$meta,sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
tpm.batch <- as.data.frame(t(read.table(opt$expression,sep=",",header=T,row.names=1, stringsAsFactors = FALSE)))
outdir1 <- opt$outdir1
outdir2<- opt$outdir2

###Subset the metasheet and match with combat results based on Timing
pre.metasheet <- subset(meta, Timing =="Pre")
post.metasheet <- subset(meta, Timing == "Post")
input.tide.pre <- as.data.frame(t(tpm.batch[meta$Timing=="Pre", ]))
input.tide.post <- as.data.frame(t(tpm.batch[meta$Timing=="Post",]))


######Write the output files to corresponding folders
write.table(data.frame("SampleName"=rownames(pre.metasheet),pre.metasheet),paste(opt$outdir1,"pre.metasheet.csv",sep=""), quote=F,sep=",",row.names=FALSE)
write.table(data.frame("gene_id"=rownames(input.tide.pre),input.tide.pre),paste(opt$outdir1,"tide.pre.txt",sep=""),quote=F,sep=",",row.names=FALSE)
write.table(data.frame("SampleName"=rownames(post.metasheet),post.metasheet),paste(opt$outdir2,"post.metasheet.csv",sep=""), quote=F,sep=",",row.names=FALSE)
write.table(data.frame("gene_id"=rownames(input.tide.post),input.tide.post),paste(opt$outdir2,"tide.post.txt",sep=""),quote=F,sep=",",row.names=FALSE)





