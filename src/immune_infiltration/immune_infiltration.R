
#############################Script to process immune infiltration #####################
suppressMessages(library(rlang))
suppressMessages(library(R.utils))
suppressMessages(library(dplyr))
suppressMessages(library(immunedeconv))
suppressMessages(library(optparse))


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
df <- read.table(input, sep = ",", row.names = 1,header = TRUE, check.names = FALSE)
df <- 2^df -1
print(head(df))

if (data_type == "mm10") {

  print("converting mouse gene symbol to human")
  #save the mouse gene expression matrix
  exp <- df

  #call the references
  genes <- read.csv(opts$reference)
  
  df <- exp[rownames(exp) %in% genes$Gene.name,]
  df <- merge(genes, df, by.x = 1, by.y = 0, all = FALSE)
  rownames(df) <- df[[2]]

  write.table(df[-1], file = "analysis/immune_infiltration/mm_to_human.tmp", sep = ",", row.names = FALSE, quote = FALSE)
  input <- "analysis/immune_infiltration/mm_to_human.tmp"
  df <- df[-c(1,2)]
}

#quantiseq
res_quant = as.data.frame(deconvolute(df, "quantiseq"))
write.table(res_quant,paste0(output,'quantiseq.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run quantiseq")


#xcell
res_xcell = as.data.frame(deconvolute(df, "xcell"))
write.table(res_xcell,paste0(output,'xcell.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run xcell")


#epic
res_epic = as.data.frame(deconvolute(df, "epic"))
write.table(res_epic,paste0(output,'epic.txt'),sep='\t',row.names=FALSE, quote = FALSE)
print("run epic")

#mcp_counter
if (data_type == "mm10") {
  #mMCP
  library(mMCPcounter)
  print('run mouse version of mcp_counter')
  res_mMcp <- mMCPcounter.estimate(exp, features = c("Gene.Symbol","ENSEMBL.ID","Probes")[1])
  res_mMcp <- cbind(rownames(res_mMcp), res_mMcp)
  colnames(res_mMcp)[1] <- "cell_type"
  write.table(res_mMcp,paste0(output,'/mcp_counter.txt'),sep='\t',row.names=FALSE, quote = FALSE)
  
} else {
  #mcp_counter
  print('run mcp_counter')
  res_mcp = as.data.frame(deconvolute(df, "mcp_counter"))
  write.table(res_mcp,paste0(output,'/mcp_counter.txt'),sep='\t',row.names=FALSE, quote = FALSE)
}

source('src/immune_infiltration/CIBERSORT.R')
res_ciber <- CIBERSORT("static/cibersort/LM22.txt",input,perm = cibersort_perm, QN = cibersort_qn, absolute = cibersort_abs, abs_method = cibersort_abs_method)
write.table(res_ciber, file = paste0(output,'cibersort_abs.txt'), sep="\t", quote = FALSE)
print("run cibersort absolution mode")

res_timer = as.data.frame(deconvolute(df, "timer",indications=rep(tolower(cancertype),ncol(df))))
write.table(res_timer,paste0(output,'timer.txt'),sep='\t', row.names=FALSE, quote = FALSE)
