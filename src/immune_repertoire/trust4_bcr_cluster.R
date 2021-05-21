suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(optparse))

option_list <- list( 
  make_option(c("-i", "--cdr3"), type="character",
              help="input cdr3 file "),
  make_option(c("-s", "--sampleid"), type="character",
              help="sample id[Required]"),
  make_option(c("-o","--outdir"),type="character",
              help="Output path [Required]")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


##set options
outdir <- opt$outdir
file <- opt$cdr3
ss <- opt$sampleid 

##set work directory
source("src/immune_repertoire/trust4_metric_functions.R")

###function of processing cdr3 for each separated sample
cdr3_process <- function(file){
  print(file)
  cdr3 <- read.table(file = file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (dim(cdr3)[1] != 0){
  
  cdr3$sample <- ss
  print(head(cdr3))
  cdr3 <- subset(cdr3, count > 0) %>% 
    mutate(V = as.character(V), J = as.character(J), C = as.character(C), CDR3aa = as.character(CDR3aa)) 

  cdr3$is_complete <- sapply(cdr3$CDR3aa, function(x) ifelse(x != "partial" && x != "out_of_frame" && !grepl("^_",x) && !grepl("^\\?", x),"Y","N"))
  cdr3.bcr <- subset(cdr3, grepl("^IG",V) | grepl("^IG",J) | grepl("^IG",C))
  } else {
  
  cdr3.bcr <- cdr3
  }
  
  ##save bcr 
  return(cdr3.bcr)
}

##main function of merging processed cdr3 data
cdr3.bcr<- cdr3_process(file)


if( dim(cdr3.bcr)[1] !=  0 ) {
##add lib size and clinic traits
cdr.bcr.new <- cdr3.bcr %>% mutate(lib.size = sum(count)) 
cdr3.bcr.heavy <- subset(cdr.bcr.new, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
print(ss)
sample_bcr_cluster <- BuildBCRlineage(sampleID = ss, Bdata = cdr3.bcr.heavy, start=3, end=10)
print ("Saving results ...")
save(sample_bcr_cluster,file = paste(outdir, "_TRUST4_BCR_heavy_cluster.Rdata", sep = ""))

} else {
sample_bcr_cluster  <- NULL
save(sample_bcr_cluster,file = paste(outdir, "_TRUST4_BCR_heavy_cluster.Rdata", sep = ""))
}

