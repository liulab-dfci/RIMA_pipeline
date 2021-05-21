suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(optparse))

option_list <- list( 
  make_option(c("-i", "--cdr3"), type="character",
              help="input cdr3 file "),
  make_option(c("-s", "--sampleid"), type="character",
              help="sample id[Required]"),
  make_option(c("-r", "--stat"), type="character",
              help="star stat file for number of reads"),
  make_option(c("-o","--outdir"),type="character",
              help="Output path [Required]")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


##set options
outdir <- opt$outdir
file <- opt$cdr3
ss <- opt$sampleid 
stat_f <- opt$stat

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

  return(cdr3.bcr)
}


##main function of merging processed cdr3 data
cdr3.bcr<- cdr3_process(file)

print(ss)
cdr3.bcr.heavy <- subset(cdr3.bcr, grepl("^IGH",V) | grepl("^IGH",J) | grepl("^IGH",C))
cdr3.bcr.light <- subset(cdr3.bcr, grepl("^IG[K|L]",V) | grepl("^IG[K|L]",J) | grepl("^IG[K|L]",C))

if( dim(cdr3.bcr.light)[1] !=  0 ) {
print ("Saving BCR light result...")
cdr3.bcr.light <- cdr3.bcr.light %>% 
                  mutate(lib.size = sum(count))                            
save(cdr3.bcr.light, file = paste(outdir, "_TRUST4_BCR_light.Rdata",sep = ""))

} else {
cdr3.bcr.light <- NULL
save(cdr3.bcr.light, file = paste(outdir, "_TRUST4_BCR_light.Rdata",sep = ""))

}




if(  dim(cdr3.bcr.heavy)[1] !=  0 ) {

print ("Saving bcr heavy results ...")
cdr3.bcr.heavy <- cdr3.bcr.heavy %>% 
                  mutate(lib.size = sum(count))  
save(cdr3.bcr.heavy, file = paste(outdir, "_TRUST4_BCR_heavy.Rdata",sep = ""))


print ("Saving cluster results ...")
sample_bcr_cluster <- BuildBCRlineage(sampleID = ss, Bdata = cdr3.bcr.heavy, start=3, end=10)
save(sample_bcr_cluster,file = paste(outdir, "_TRUST4_BCR_heavy_cluster.Rdata", sep = ""))



print ("Saving clonality results ...")
sample_all_clonality <- getClonality(sampleID = ss, cdr3.bcr.heavy, start=3, end=10)
save(sample_all_clonality,file = paste(outdir, "_TRUST4_BCR_heavy_clonality.Rdata", sep = ""))


 
print ("Saving SHM results ...")
if( is.null(sample_bcr_cluster)){
	ss.ratio <- NULL
  } else if (is.na(sample_bcr_cluster)) {
    ss.ratio <- NULL
  } else { ss.ratio <- getSHMratio(sample_bcr_cluster) } 
save(ss.ratio,file = paste(outdir, "_TRUST4_BCR_heavy_SHMRatio.Rdata", sep = ""))  



print ("Saving library reads results ...")
stat<- read.table(stat_f,sep = "\t",row.names = 1)
map.reads <- stat["reads mapped:","V2"]
lib.size <- mean(cdr3.bcr.heavy$lib.size)
Infil <- signif(as.numeric(lib.size)/as.numeric(map.reads),4)
bcr.lib.reads <- cbind(map.reads,lib.size,Infil) 
save(bcr.lib.reads,file = paste(outdir, "_TRUST4_BCR_heavy_lib_reads_Infil.Rdata", sep = "")) 

} else {
sample_all_clonality <- NULL
sample_bcr_cluster  <- NULL
ss.ratio <- NULL
bcr.lib.reads <- NULL
cdr3.bcr.heavy <- NULL


save(sample_bcr_cluster,file = paste(outdir, "_TRUST4_BCR_heavy_cluster.Rdata", sep = ""))
save(sample_all_clonality,file = paste(outdir, "_TRUST4_BCR_heavy_clonality.Rdata", sep = ""))
save(ss.ratio,file = paste(outdir, "_TRUST4_BCR_heavy_SHMRatio.Rdata", sep = "")) 
save(bcr.lib.reads,file = paste(outdir, "_TRUST4_BCR_heavy_lib_reads_Infil.Rdata", sep = "")) 
save(cdr3.bcr.heavy, file = paste(outdir, "_TRUST4_BCR_heavy.Rdata",sep = ""))

}

