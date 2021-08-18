#!/usr/bin/env Rscript
#dependencies
library(maftools)
library(dplyr)
library(stringr)
library(optparse)
library(ggplot2)

#make option list and parse command line
option_list = list(
  make_option(c("-m", "--maf"), type="character", default=NULL, 
              help="maf file", metavar="character"),
  make_option(c("-t", "--top_n"), type="character", default=NULL, 
              help="top n genes for showing in oncoplot", metavar="character"),
  make_option(c("-n", "--n_pattern"), type="character", default=NULL, 
              help="top n genes for showing mutation pattern", metavar="character"),
  make_option(c("-c", "--cohort_outdir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-s", "--metasheet"), type="character", default=NULL, 
              help="output directory", metavar="character")  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###option setting
top_n <- as.numeric(opt$top_n)
top_n_pattern <- as.numeric(opt$n_pattern)
outdir <- opt$cohort_outdir
meta <- opt$metasheet
maf_file <- read.maf(maf=opt$maf,useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
metasheet <- read.csv(meta)

###test
# maf_file <- read.maf(maf="~/Desktop/merged_processed_dbSNP.maf",useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
# metasheet <- read.csv("~/Desktop/metasheet.csv")

###set color series
  all.samples <- as.character(maf_file@variants.per.sample$Tumor_Sample_Barcode)
  color_list <- c()
  for (ss in all.samples){
    maf <- subsetMaf(maf_file, tsb = c(ss))
    gs = getGeneSummary(maf)
    gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE][,-1]
    color_list <- c(colnames(gs),color_list)
  }
  color_list <- unique(color_list)
  col = RColorBrewer::brewer.pal(n = length(color_list), name = 'Paired')
  names(col) = color_list
  
  ##write somatic mutation maf table
  merge_data <- merge(maf_file@data,metasheet,by.x = "Tumor_Sample_Barcode", by.y = "SampleName")
  maf.tab <- subset(merge_data, select = c('PatName', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode','Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 
                                              'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 
                                              'Tumor_Seq_Allele2', 'HGVSc', 'HGVSp', 'HGVSp_Short'))
  
  write.table(maf.tab, file = paste(outdir, "/", "maftools_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  
  png(file = paste(outdir, "/", "maftools_summary_plot.png", sep = ""),res = 300, width = 2000, height = 1500)
  plotmafSummary(maf = maf_file, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, color = col, textSize = 5, titleSize = c(1.2,1))
  dev.off()
  ###transition and tranversion rate
  png(file = paste(outdir, "/", "maftools_titv_plot.png", sep = ""), res = 300, width = 1200, height = 1200)
  titv = titv(maf = maf_file, plot = FALSE, useSyn = TRUE)
  plotTiTv(res = titv)
  dev.off()
  ###oncoplot
  png(file = paste(outdir, "/", "maftools_onco_plot.png", sep = ""), res = 300, width = 2500, height = 3000)
  oncoplot(maf = maf_file, top = top_n, removeNonMutated = FALSE,fontSize = 1.1, legendFontSize = 1.5,colors = col,gene_mar = 10)
  dev.off()
  ###lollipop plot
  onco.mat <- subset(maf_file@data, select = c('Hugo_Symbol','Tumor_Sample_Barcode'))
  onco.mat <- onco.mat[!duplicated(onco.mat),]
  mut.freq <- onco.mat %>% group_by(Hugo_Symbol) %>% dplyr::summarise(freq = n()) %>% arrange(desc(freq))
  if (length(mut.freq$Hugo_Symbol)<=top_n_pattern){
    top.genes <- mut.freq$Hugo_Symbol
  } else {top.genes <- mut.freq$Hugo_Symbol[1:top_n_pattern]}
  pdf(file = paste(outdir, "maftools_lollipop_plot.pdf", sep = ""), width = 6, height = 3)
  gp <- lapply(top.genes, function(g){
    lollipopPlot(maf_file, AACol = "HGVSp_Short", gene = g, showMutationRate = TRUE,showDomainLabel=FALSE)
  })
  dev.off()
  ###compare with TCGA data
  png(file = paste(outdir, "/", "maftools_tcga_plot.png", sep = ""), res = 200, width = 1500, height = 800)
  tcgaCompare(maf = maf_file, cohortName = 'USER', col = c("gray70", "red"), medianCol = "black")
  dev.off()
  ###Drug-Gene Interactions
  png(file = paste(outdir, "/", "maftools_drug_plot.png", sep = ""), res = 200, width = 1500, height = 800)
  drugInteractions(maf = maf_file, fontSize = 0.75)
  dev.off()
