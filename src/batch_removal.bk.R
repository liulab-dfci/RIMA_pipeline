#!/usr/bin/env Rscript

#################################################################
# Function: Remove batch effect
# Call: Rscript combat.r -e expression_dat -b batch_dat 
# R packages used: sva,optparse
#################################################################

## install necessary libraries
library(sva)
library(limma)
library(optparse)

# make option list and parse command line
option_list <- list(  
  make_option(c("-e", "--expression_dat"), type="character", 
              help="Input path of expression file. [Required]"),
  make_option(c("-b", "--batch_dat"), type="character",
              help="Input path of corresponding batch file[Required]"),
  make_option(c("-m", "--method"), type="character",
              help="batch removal method"),
  make_option(c("-c", "--covariates"), type="character",  
              help="covariates needs to be adjusted for"),
  make_option(c("-n", "--covars_n"), type="numeric",
              help="number of covariates needs to be preserved when correded batches effect"),
  make_option(c("-o","--output",type="character", help="Output files [Required]"))
)
opt_parser <- OptionParser(option_list=option_list);
opts <- parse_args(opt_parser);

# paramenter checking
if(is.null(opts$expression_dat)) stop('Expression file required.')  ###if not provide batch file output log expression matrix

###functions for inputing and outputing
ssgsvaFormat <- function(dat){
  dat <- cbind(gene_id=rownames(dat),dat)
  return(dat)
}

writeDF <- function(dat,path){
  write.table(dat,path,quote = FALSE, sep=',', row.names = FALSE)
}

# load data
#if(file_test('-f',opts$expression_dat)){
  expr.dat <- read.table(opts$expression_dat,sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
  expr.dat <- log2(expr.dat + 1)
#} else {
#  expr.dat <- getExpr(opts$expression_dat) 
#}

print('loading data done!')
CVFILTER <- 0
##filtering out genes with low variance among samples
# datatype <- opts$data_type
## Calculate CVs for all genes (rows) - ComBat needs a minimal variance, so selecting ones that have a min var
# if (datatype == "cufflinks") {CVFILTER = 0.01}
# if (datatype == "rsem || salmon") {CVFILTER = 1}
# if (datatype == "star") {CVFILTER = 2}
mean_nolym <- apply(expr.dat,1,mean)
var_nolym <- apply(expr.dat,1,var)
cv_nolym <- abs(var_nolym/mean_nolym)
filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)

## Select those genes that pass variance filtering
exprZero <- expr.dat
expr.dat <- expr.dat[rownames(expr.dat) %in% names(filt_genes),]
exprZero <- subset(exprZero, !(rownames(exprZero) %in% names(filt_genes)))

print('start removal')

## do batch or not
#if(file_test('-f',opts$batch_dat)){
  batch.dat = read.table(opts$batch_dat, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
  batch.dat[,'Gender']
print('metasheet load done!')
 # Get overlap samples
  overlap_sample = intersect(colnames(expr.dat),rownames(batch.dat))
  
  # stop at low sample numbers
  if(length(overlap_sample) < dim(expr.dat)[2]/2) stop("too few samples")
  
  print('Processing samples:')
  print(overlap_sample)
  
  
  # remove batch effect
  errflag = F
#  if(opts$method == "combat"){
#    expr.combat = tryCatch(ComBat(as.matrix(expr.dat[,overlap_sample]),
#                                  batch.dat[overlap_sample,"Batch"]),
#                           error = function(e){
#                             print(e)
#                             errflag <<- T
#                           }
#    )
#  }
#  if(opts$method == "limma"){
    covars <- strsplit(opts$covariates, "\\,")[[1]]
    print(covars)
    covariates_df <- batch.dat[overlap_sample,covars]
    
    fomul <- as.formula(paste("", paste0(covars, collapse = " + "), sep = " "))
    print(fomul)
    
    design <- model.matrix(fomul, data=covariates_df)
    
    if(opts$covars_n == 0){  ##no optional design matrix to be preserved
      batch.design <- design[,-c(1:1)]
      expr.combat = tryCatch(removeBatchEffect(as.matrix(expr.dat[,overlap_sample]), 
                                               covariates = batch.design),
                             error = function(e){
                               print(e)
                               errflag <<- T
                             })
    } else {
      treatment.design <- design[,c(1:(opts$covars_n+1))]
      batch.design <- design[,-c(1:(opts$covars_n+1))]
      expr.combat = tryCatch(removeBatchEffect(as.matrix(expr.dat[,overlap_sample]), 
                                               design =  treatment.design,
                                               covariates = batch.design),
                             error = function(e){
                               print(e)
                               errflag <<- T
                             })
    }
 
#  }
#  print(errflag)
#  if(errflag){
#    print('Cannot do Combat with batch removal error!')
#    expr.dat = rbind(expr.dat,exprZero)
#    writeDF(ssgsvaFormat(expr.dat),paste0(opts$output,'.batch'))
#  } else {
    expr.combat = rbind(expr.combat,exprZero)
    writeDF(ssgsvaFormat(expr.combat),paste0(opts$output,'.batch'))
#  }
#} else{
#  print('Did not do Combat without bacth file!')
#  writeDF(ssgsvaFormat(expr.dat),paste0(opts$output,'.batch'))
#}


