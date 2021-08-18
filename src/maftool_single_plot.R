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
  make_option(c("-s", "--ss_outdir"), type="character", default=NULL, 
              help="single sample plot output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###read maf file
#maf_file <- read.maf(maf=opt$maf,useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
tryCatch({
  maf_file <- read.maf(maf=opt$maf,useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
  return(maf_file)}, 
  error = function(e){cat("ERROR:", conditionMessage(e),"\n")})

if (exists("maf_file")){
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

###single plot function
dashboard_single_sample = function(maf, res, color = NULL, rmOutlier = TRUE, log_conv = FALSE, titv.color = NULL, 
                                   sfs = statFontSize, fontSize = fs, n = 10, donut = pie, rawcount = TRUE, stat = NULL, 
                                   titleSize = NULL, barcodes = NULL, barcodeSize = NULL, 
                                   plotType = "both", sampleOrder = NULL,  showBarcodes = FALSE, 
                                   textSize = 0.8, baseFontSize = 1, axisTextSize = c(1, 1), plotNotch = FALSE){
  
  if(is.null(color)){
    #hard coded color scheme if user doesnt provide any
    col = get_vcColors()
  }else{
    col = color
  }
  
  vcs = getSampleSummary(maf)
  vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]
  vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event
  vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
  colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')
  
  data.table::setDF(vcs)
  rownames(x = vcs) = vcs$Tumor_Sample_Barcode
  if (length(colnames(vcs)==2)){
    
    colname <- rownames(vcs)
    rowname <- colnames(vcs)
  } 
  vcs = vcs[,-1]
  vcs = t(vcs)
  if (length(colnames(vcs))==0){
    colnames(vcs) <- colname
    rownames(vcs) <- rowname[-1]
  }
  
  lo = matrix(data = 1:4, nrow = 1, byrow = TRUE)
  graphics::layout(mat = lo, heights = c(1), widths = c(3, 2, 2, 2))
  par(cex.axis = fontSize, font = 3, cex.main = titleSize[1], lwd = 1.2)
  
  #--------------------------- variant classification plot -----------------
  vc.plot.dat = rev(rowSums(vcs))
  if(log_conv){
    vc.plot.dat = log10(vc.plot.dat)
  }
  
  xt = pretty(c(0, vc.plot.dat))
  
  par(mar = c(3, 9, 3, 1))
  b = barplot(vc.plot.dat, axes = FALSE, horiz = TRUE, col = col[names(vc.plot.dat)], border = NA,
              xlim = c(0, max(xt)), names.arg = rep("", length(vc.plot.dat)))
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  axis(side = 2, at = b, labels = names(vc.plot.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.9, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Classification", adj = 0, cex.main = titleSize[1], font = 3)
  if(log_conv){
    axis(side = 2, at = 0, labels = "(log10)", lwd = 1.2, font = 3,
         las = 1, cex.axis = fontSize*0.9, hadj = 0.5, padj = 2, line = 0.75, tick = FALSE, outer = FALSE)
  }
  
  
  
  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.frame(value = colSums(titv.counts[,2:7]), stringsAsFactors = FALSE)
  titv.sums$class = rownames(titv.sums)
  if(!rawcount){
    titv.sums$raw_value = titv.sums$value
    titv.sums$value = titv.sums$value/sum(titv.sums$value)
    xt = seq(0, 1, 0.25)
  }else{
    xt = as.integer(seq(0, max(titv.sums$value, na.rm = TRUE), length.out = 4))
  }
  
  if(is.null(titv.color)){
    titv.color = get_titvCol()
  }
  
  par(mar = c(3, 3, 3, 1))
  b = barplot(titv.sums$value, axes = FALSE, horiz = TRUE, col = titv.color[rownames(titv.sums)],
              border = NA, xlim = c(0, xt[length(xt)]))
  axis(side = 2, at = b, labels = rownames(titv.sums), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "SNV Class", adj = 0, cex.main = titleSize[1], font = 3)
  if(!rawcount){
    text(x = titv.sums$value+0.03, y = b, labels = titv.sums$raw_value,
         font = 4, col = "black", cex = fontSize, adj = 0)
  }
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  
  
  
  #--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
  nsamps = as.numeric(maf@summary[ID %in% "Samples", summary])
  gs.load = gs[,.(Hugo_Symbol, AlteredSamples)]
  gs.load[,AlteredSamples := round(AlteredSamples/nsamps, digits = 2) * 100]
  data.table::setDF(x = gs.load, rownames = gs.load$Hugo_Symbol)
  gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]
  
  if(nrow(gs) < n){
    gs.dat = gs
  }else{
    gs.dat = gs[1:n]
  }
  
  data.table::setDF(gs.dat)
  rownames(gs.dat) = gs.dat$Hugo_Symbol
  if (length(colnames(gs.dat)==2)){
    
    colName <- rownames(gs.dat)
    rowName <- colnames(gs.dat)
  } 
  gs.dat = gs.dat[,-1]
  gs.dat = t(gs.dat)
  
  if (length(colnames(gs.dat))==0){
    colnames(gs.dat) <- colName
    rownames(gs.dat) <- rowName[-1]
  }
  gs.dat = gs.dat[names(sort(rowSums(gs.dat), decreasing = TRUE)),, drop = FALSE]
  gs.dat = gs.dat[,names(sort(colSums(gs.dat))), drop = FALSE]
  
  xt = as.integer(seq(0, max(colSums(gs.dat))+2, length.out = 4))
  
  par(mar = c(3, 4, 3, 0))
  gs.load = gs.load[rev(colnames(gs.dat)),,]
  b = barplot(gs.dat, axes = FALSE, horiz = TRUE, col = col[rownames(gs.dat)], border = NA,
              xlim = c(0, max(xt)+(max(xt)*0.15)), names.arg = rep("", ncol(gs.dat)))
  axis(side = 2, at = b, labels = colnames(gs.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = paste0('Top ',  n, '\nmutated genes'), adj = 0, cex.main = titleSize[1], font = 3)
  text(x = colSums(gs.dat)+1, y = b, labels = rev(paste0(gs.load$AlteredSamples, "%")),
       font = 4, col = "black", cex = fontSize*0.9, adj = 0)
  abline(h = b, v = xt,lty = 2, lwd = 0.3,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  
  #--------------------------- Ti-Tv plot -----------------
  col = get_titvCol(alpha = 0.8)
  titv.frac = res$fraction.contribution
  titv.frac.melt = data.table::melt(data = titv.frac, id = "Tumor_Sample_Barcode")
  conv.class = c("Ti", "Ti", "Tv", "Tv", "Tv", "Tv")
  names(conv.class) = c("T>C", "C>T", "T>A", "T>G", "C>A", 
                        "C>G")
  titv.frac.melt$TiTv = conv.class[as.character(titv.frac.melt$variable)]
  data.table::setDT(x = res$TiTv.fractions)
  titv.contrib = suppressMessages(data.table::melt(res$TiTv.fractions, 
                                                   id = "Tumor_Sample_Barcode"))
  titv.frac.melt$variable = factor(x = titv.frac.melt$variable, 
                                   levels = c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G"))
  titv.order = titv.frac.melt[, mean(value), by = .(variable)]
  titv.order = titv.order[order(V1, decreasing = TRUE)]
  orderlvl = as.character(titv.order$variable)
  titv.frac.melt$variable = factor(x = titv.frac.melt$variable, 
                                   levels = orderlvl)
  tf = res$TiTv.fractions
  data.table::setDF(x = tf)
  rownames(tf) = tf$Tumor_Sample_Barcode
  tf = tf[, -1]
  par(mar = c(3, 4, 3, 1))
  b = barplot(c(tf[,1], tf[,2]), axes = FALSE, xlab = "", ylab = "", col = c("gold","lightblue"), border = NA)
  axis(side = 1, at = c(1,2), lwd = 1.2, font = 3, las = 2, labels = names(tf),
       line = -0.3, hadj = 0.6, cex.axis = fontSize, tick = FALSE)
  axis(side = 2, at = as.integer(seq(0, max(tf), length.out = 4)), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3)
  
  title(main = "Ti and Tv", adj = 0.5, cex.main = titleSize[1], font = 3)
  abline(h = as.integer(seq(0, max(tf), length.out = 6)), v = seq(0,2,length.out = 3),lty = 2, lwd = 0.3,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  
}

get_titvCol = function(alpha = 1){
  col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}

###single sample figures
#for (ss in all.samples){}

lapply(all.samples, function(ss){
  tryCatch({
    ss.maf <- subsetMaf(maf_file, tsb = c(ss))
    ss.titv <- titv(maf = ss.maf, plot = FALSE, useSyn = TRUE)
    # maf = ss.maf
    rmOutlier = TRUE
    addStat = 'median'
    titvRaw = FALSE
    color = col
    textSize = 5
    titleSize = c(1.2,1)
    fs = 1
    pie = FALSE
    log_scale = FALSE
    titvColor = NULL
    top = 10
    
    ###write a mutation table for a single sample
    #     ss.maf.sub <- subset(ss.maf@data, select = c('Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2','Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'HGVSc', 'HGVSp', 'HGVSp_Short'))
    #     write.table(ss.maf.sub, file = paste("~/Documents/rnaseq/bootstrap/images/mutation/",ss, "_maf_table.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
    
    ###overview of a single sample
    png(file = paste(opt$ss_outdir,"_maftools_summary_plot.png", sep = ""),res = 300, width = 2000, height = 1000)
    p <- dashboard_single_sample(maf = ss.maf, res = ss.titv, color = col, rmOutlier = TRUE, 
                            log_conv = log_scale, titv.color = titvColor, fontSize = fs, 
                            titleSize = titleSize, sfs = statFontSize, n = top, 
                            donut = pie, rawcount = titvRaw, stat = addStat, 
                            barcodes = showBarcodes, barcodeSize = textSize,plotType = "both", sampleOrder = NULL,  
                            showBarcodes = FALSE, textSize = 0.8, baseFontSize = 1, axisTextSize = c(1, 1), plotNotch = FALSE)
    print(p)
    dev.off()
  }, error = function(e) {cat("ERROR:",conditionMessage(e),"\n")})
})
} else {
  png(file = paste(opt$ss_outdir,"_maftools_summary_plot.png", sep = ""),res = 300, width = 2000, height = 1000)
  print(ggplot())
  dev.off()
}
  

