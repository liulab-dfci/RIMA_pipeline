#!/usr/bin/env Rscript

#dependencies
library(dplyr)
library(ggrepel)
library(optparse)
library(ggplot2)
library(tidyverse)
library(data.table)

option_list = list(
  make_option(c("-f", "--fusion"), type="character", default=NULL, 
              help="merged fusion prediction file", metavar="character"),
  make_option(c("-out", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
  
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


###read in input data
fusion <- read.table(file = opt$fusion, header = TRUE)
outdir <- opt$outdir

##########
fusion=setDT(fusion)[, paste0("FusionName", 1:2) := tstrsplit(FusionName, "--")]
fusion=subset(fusion,select=c("FusionName1","FusionName2"))
fusion=fusion[!duplicated(fusion), ]
write.table(fusion, paste(outdir, "pyprada_fusion_table.txt", sep = ""), sep = "\t", quote = FALSE,row.names=FALSE,col.names=FALSE)

