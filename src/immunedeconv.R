library(rlang)
library(R.utils)
library(dplyr)
library(immunedeconv)
library(optparse)



# make option list and parse command line
library(optparse)
option_list <- list(  
  make_option(c("-e", "--expression_dat"), type="character", 
              help="Input path of expression file. [Required]"),
  make_option(c("-o", "--output_dat"), type="character",
              help="Ouput path of corresponding tables. [Required]"),
  make_option(c("-p", "--permutation"), type="numeric", default=NULL, 
              help="permutation for cibersort", metavar="character"),
  make_option(c("-q", "--quantile_normalization"), type="character", default=NULL, 
              help="whether enable quantile normalization for cibersort", metavar="character"),
  make_option(c("-a", "--absolute"), type="character", default=NULL, 
              help="absolute mode for cibersort", metavar="character"),
  make_option(c("-m", "--abs_method"), type="character", default=NULL, 
              help="absolute method for cibersort", metavar="character"),
  make_option(c("-t", "--cancer_type"), type="character",
              help="cancer type used in timer.  [Required]"),
  make_option(c("-d", "--data_type"), type="character", default=NULL, 
              help="data_type for immnuedecov", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="mouse to human ref", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);



# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')

input <- opts$expression_dat
output <- opts$output_dat
cancertype <- opts$cancer_type
cibersort_perm <- opts$permutation
cibersort_qn <- opts$quantile_normalization
cibersort_abs <- opts$absolute
cibersort_abs_method <- opts$abs_method
data_type <- opts$data_type

###read in data
df <- read.table(input, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)

if (data_type == "mm10") {

  print("converting mouse gene symbol to human")
  #save the mouse gene expression matrix
  exp <- df

  #call the references
  genes <- read.csv(opts$reference)
  
  df <- exp[rownames(exp) %in% genes$Gene.name,]
  df <- merge(genes, df, by.x = 1, by.y = 0, all = FALSE)
  rownames(df) <- df[[2]]

  write.table(df[-1], file = paste0(output, "/mm_to_human.tmp", sep = ""), sep = ",", row.names = FALSE, quote = FALSE)
  input <- paste0(output, "/mm_to_human.tmp")
  df <- df[-c(1,2)]
}
# df <- read.table("~/Documents/rnaseq/data/tpm_convertID.txt", sep = ",", row.names = 1,header = TRUE, check.names = FALSE)

#quantiseq
res_quant = as.data.frame(deconvolute(df, "quantiseq"))
write.table(res_quant,paste0(output,'/quantiseq.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run quantiseq")


#xcell
res_xcell = as.data.frame(deconvolute(df, "xcell"))
write.table(res_xcell,paste0(output,'/xcell.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run xcell")


#epic
res_epic = as.data.frame(deconvolute(df, "epic"))
write.table(res_epic,paste0(output,'/epic.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run epic")

if (data_type == "mm10") {
  #mMCP
  library(mMCPcounter)
  res_mMcp <- mMCPcounter.estimate(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1])
  res_mMcp <- cbind(rownames(res_mMcp), res_mMcp)
  colnames(res_mMcp)[1] <- "cell_type"
  write.table(res_mMcp,paste0(output,'/mcp_counter.txt'),sep='\t',row.names=FALSE, quote = FALSE)
  
} else {
  #mcp_counter
  res_mcp = as.data.frame(deconvolute(df, "mcp_counter"))
  write.table(res_mcp,paste0(output,'/mcp_counter.txt'),sep='\t',row.names=FALSE, quote = FALSE)
  print('run mcp_counter')
}

#cibersort
source('src/CIBERSORT.R')
res_ciber <- CIBERSORT("static/immunedeconv/cibersort/LM22.txt",input,perm = cibersort_perm, QN = cibersort_qn, absolute = cibersort_abs, abs_method = cibersort_abs_method)
write.table(res_ciber, file = paste0(output,'/cibersort_abs.txt'), sep="\t", quote = FALSE)
print("run cibersort absolution mode")

##timer
## colnames(df)<-gsub('[.]','-',colnames(df))
#cancer_list <- c("BLCA", "BRCA", "CESC", "COADREAD", "ESCA", "GBM", "HNSC", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", 
#                 "LUSC", "OV", "PAAD", "SARC", "SKCM", "STAD", "UCEC", "UVM")
#if (cancertype %in% cancer_list) {
#  res_timer = as.data.frame(deconvolute(df, "timer",indications=rep(tolower(cancertype),ncol(df))))
#  write.table(res_timer,paste0(output,'/timer.txt'),sep='\t', row.names=FALSE, quote = FALSE)
#} else {
#  write.table(NULL,paste0(output,'/timer.txt'),sep='\t', row.names=FALSE, quote = FALSE)
#  print(paste(cancertype, "is not available for TIMER, please check the document of TIMER"))
#}





