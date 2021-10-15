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
  cdr3.tcr <- subset(cdr3, grepl("^TR",V) | grepl("^TR",J) | grepl("^TR",C))
  } else {  
  cdr3.tcr <- cdr3
  }
  return(cdr3.tcr)
}



##main function of merging processed cdr3 data
cdr3.tcr<- cdr3_process(file)


if( dim(cdr3.tcr)[1] !=  0 ) {
print ("Saving TCR result...")
cdr3.tcr<- cdr3.tcr %>% 
                  mutate(lib.size = sum(count))                            
save(cdr3.tcr, file = paste(outdir, "_TRUST4_TCR.Rdata",sep = ""))
write.table(cdr3.tcr, paste(outdir, "_TRUST4_TCR.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)


print ("Saving clonality results ...")
tcr_clonality <- getClonalityTCR(sampleID = ss,cdr3.tcr)
save(tcr_clonality,file = paste(outdir, "_TRUST4_TCR_clonality.Rdata", sep = ""))
write.table(t(as.matrix(tcr_clonality)),paste(outdir, "_TRUST4_TCR_clonality.txt", sep = ""),sep = '\t',quote=FALSE,row.names = FALSE)




print ("Saving library reads results ...")
stat<- read.table(stat_f,sep = "\t",row.names = 1)
map.reads <- stat["reads mapped:","V2"]
lib.size <- mean(cdr3.tcr$lib.size)
Infil <- signif(as.numeric(lib.size)/as.numeric(map.reads),4)
tcr.lib.reads <- cbind(sample= ss, map.reads,lib.size,Infil) 
save(tcr.lib.reads,file = paste(outdir, "_TRUST4_TCR_lib_reads_Infil.Rdata", sep = "")) 
write.table(tcr.lib.reads,paste(outdir, "_TRUST4_TCR_lib_reads_Infil.txt", sep = ""),sep = '\t',quote=FALSE,row.names = FALSE)


} else {

cdr3.tcr <- data.frame(count=integer(),
                 frequency = double(),
                 CDR3nt = character(),
                 CDR3aa = character(),
                 V = character(),
                 D = character(),
                 J = character(),
                 cid = character(),
                 sample = character(),
                 is_complete = character(),
                 lib.size = integer() )
tcr_clonality  <- data.frame(
                 sample = character(),
                 clonality = character()
                 )
tcr.lib.reads <- data.frame(
                 sample = character(),
                 map.reads = character(),
                 lib.size = character(),
                 Infil = character()
                 )


save(cdr3.tcr, file = paste(outdir, "_TRUST4_TCR.Rdata",sep = ""))
write.table(cdr3.tcr, paste(outdir, "_TRUST4_TCR.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)
save(tcr_clonality,file = paste(outdir, "_TRUST4_TCR_clonality.Rdata", sep = ""))
write.table(tcr_clonality, paste(outdir, "_TRUST4_TCR_clonality.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)
save(tcr.lib.reads,file = paste(outdir, "_TRUST4_TCR_lib_reads_Infil.Rdata", sep = "")) 
write.table(tcr.lib.reads, paste(outdir, "_TRUST4_TCR_lib_reads_Infil.txt",sep=''),sep = '\t',quote=FALSE,row.names = FALSE)

}

