library(data.table)
library(dplyr)
library(ggpubr)
library(data.table)
library(pheatmap)
library(optparse)
# make option list and parse command line
option_list <- list(  
  make_option(c("-i","--input"), type="character", 
              help="Path of merged microbiota abundance data. [Required]"),
  make_option(c("-o","--outdir"),type="character", 
              help="Output files [Required]"),
  make_option(c("-c","--clinic_col"),type = "character",
              help = "Phenotypes to be compared"),
  make_option(c("-m","--meta"),type = "character",
              help = "Path of meta data")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###test
# outdir<-"~/Desktop/"
# clinic.col <- "Response"
# clinic_col <- "Response"
# input<-"~/Desktop/merged_microbiota_abundance.txt"
# meta<-"~/Desktop/metasheet.csv"


####set parameters
outdir<-opt$outdir
clinic.col <- unlist(strsplit(opt$clinic_col,","))
input <- opt$input
meta <- opt$meta

####read in merged microbiota abundance data
abundance.input <- fread(input, sep = "\t", header = TRUE)

####read in meta data
meta <- read.table(meta,sep = ",", header = TRUE, row.names = 1)

if ("abundance" %in% colnames(abundance.input)){
  ####process and filter data(centrifuge)#####
  abundance <- subset(abundance.input, select = c("sample","name","abundance"),!(name %in% c("synthetic construct","Homo sapiens")))   ###, )
  abundance.nor <- abundance %>% group_by(sample) %>% mutate(abundance.rate = abundance/sum(abundance))
  abundance.new <- reshape2::dcast(abundance.nor,sample~name, value.var = "abundance.rate",fun.aggregate=mean)
  rownames(abundance.new) <- abundance.new$sample
  abundance.new <- abundance.new[,-1]
  abundance.new[is.na(abundance.new)] <- 0
  ##diversity of abundance
  aver.abundance <- apply(abundance.new,2,mean)
  ###filter out microbes without abundance in all samples
  micros <- names(aver.abundance)[which(aver.abundance > 0)]
  m.mad <- apply(abundance.new[,micros],2,mad)
  micros.sel <- micros[which(m.mad > 0)]
  ###output selected microbes with abundance 
  micros.out <- abundance.input %>% dplyr::filter(name %in% micros.sel) %>% dplyr::select(name, taxID) %>% distinct()
  write.table(micros.out, file = paste(outdir, "selected_microbes_taxonID.txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
} else {
  ####process and filter data(pathseq)#####
  abundance <- subset(abundance.input, select = c("sample","name","score_normalized"),!(name %in% c("synthetic construct","Homo sapiens")))  
  colnames(abundance)[3] <- "abundance"
  abundance.nor <- abundance %>% group_by(sample) %>% mutate(abundance.rate = abundance/sum(abundance))
  abundance.new <- reshape2::dcast(abundance.nor,sample~name, value.var = "abundance.rate",fun.aggregate=mean)
  rownames(abundance.new) <- abundance.new$sample
  abundance.new <- abundance.new[,-1]
  abundance.new[is.na(abundance.new)] <- 0
  ##diversity of abundance
  aver.abundance <- apply(abundance.new,2,mean)
  ###filter out microbes without abundance in all samples
  micros <- names(aver.abundance)[which(aver.abundance > 0)]
  m.mad <- apply(abundance.new[,micros],2,mad)
  micros.sel <- micros[which(m.mad > 0)]
  ###output selected microbes with abundance 
  micros.out <- abundance.input %>% dplyr::filter(name %in% micros.sel) %>% dplyr::select(name, tax_id) %>% distinct()
  write.table(micros.out, file = paste(outdir, "selected_microbes_taxonID.txt", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}



###########---------------------- show abundance of selected microbes in all samples --------------------###########
###prepare data
abun.ss <- abundance.new[,micros.sel]
###set colors and breaks for heatmap
bk1 <- c(seq(0,0.01,by=0.0001))
bk2 <- c(seq(0.011,max(abun.ss),by = 0.01))
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing
my_palette1 <- c(colorRampPalette(colors = c("white", "gold"))(n = length(bk1)))   #c(colorRampPalette(colors = c("darkred", "tomato1"))(n = length(bk2)-1))
my_palette2 <- c(colorRampPalette(colors = c("gold", "red"))(n = length(bk2)))   #c(colorRampPalette(colors = c("darkred", "tomato1"))(n = length(bk2)-1))
my_palette <- c(my_palette1,my_palette2)
###set annotation for samples
for (clinic_col in clinic.col){
  
  annot_col1 <- data.frame(Triat = as.character(meta[rownames(abun.ss),clinic_col]),
                           row.names = rownames(abun.ss))
  ###sample group
  phenos <-sort(unique(as.character(meta[,clinic_col])),decreasing = TRUE)
  ###if only one phenotype for a given clinical data
  Pwidth <- nrow(abun.ss)*5+1200
  Pheight <- 30*ncol(abun.ss)
  if(length(phenos) > 1){
    png(paste(outdir,clinic_col, "_microbes_abundance.png", sep = ""), width = Pwidth, height = Pheight)
    # if(length(phenos) != 2){
    #   pheatmap(t(abun.ss), color = my_palette, breaks = bk,border_color = NA, #show_rownames = FALSE,
    #            annotation_col = annot_col1,
    #            fontsize_row = 8, fontsize_col = 8)
    # }
    ####if there two phenotypes for a given clinical data calculating significance level 
    if(length(phenos) == 2){
      tmp <- cbind.data.frame(sample = rownames(abun.ss), abun.ss)
      abun.ss.melt <- reshape2::melt(tmp, id.var = "sample") %>% mutate(Group = as.character(meta[sample, clinic_col]))
      df.tmp <- abun.ss.melt %>% group_by(variable,Group) %>% dplyr::summarise(Mean = mean(value))
      df.tmp.comp <- unlist(lapply(as.character(unique(df.tmp$variable)), function(c){
        gap <- subset(df.tmp, variable == c & Group == phenos[1], select = "Mean") - subset(df.tmp, variable == c & Group == phenos[2], select = "Mean")
        gap <- gap[1,1]
        if(gap >0) {return(paste("Higher_in_",phenos[1], sep = ""))} else if(gap < 0){return (paste("Higher_in_",phenos[2], sep = ""))} else{return ("equal")}
      }))
      names(df.tmp.comp) <- as.character(unique(df.tmp$variable))
      comp.sig <- data.frame(compare_means(value~Group, data = abun.ss.melt, group.by = "variable", paired = FALSE)) %>%
        mutate(Means = df.tmp.comp[as.character(variable)]) %>%
        mutate(Means_sig = ifelse(p.signif != "ns", paste(Means, p.signif, sep = ""), p.signif))
      rownames(comp.sig) <- as.character(comp.sig$variable)
      ###add annotation for significance level
      annot_col2 <- data.frame(Feature = as.character(comp.sig[colnames(abun.ss),"Means_sig"]),
                               row.names = colnames(abun.ss))
      ###set annotation color
      sig_col <- c("#FED976","#FD8D3C","#E31A1C","#800026","lightgray","#C7E9B4","#41B6C4","#225EA8","#081D58")
      names(sig_col) <- c(paste(paste("Higher_in_",phenos[1], sep = ""), c("*","**","***","****"), sep = ""), "ns",
                          paste(paste("Higher_in_",phenos[2], sep = ""), c("*","**","***","****"), sep = ""))
      ann_colors = list(Feature = sig_col)
    } else {
      tmp <- cbind.data.frame(sample = rownames(abun.ss), abun.ss)
      abun.ss.melt <- reshape2::melt(tmp, id.var = "sample") %>% mutate(Group = as.character(meta[sample, clinic_col]))
      comp.sig <- data.frame(compare_means(value~Group, data = abun.ss.melt, group.by = "variable", paired = FALSE, method = "kruskal.test")) 
      rownames(comp.sig) <- as.character(comp.sig$variable)
      sig_col <- c("#FED976","#FD8D3C","#E31A1C","#800026","lightgray")
      names(sig_col) <- c("*","**","***","****","ns")
      ann_colors = list(Feature = sig_col)
      annot_col2 <- data.frame(Feature = as.character(comp.sig[colnames(abun.ss),"p.signif"]),
                               row.names = colnames(abun.ss))
    }
    pheatmap(t(abun.ss), color = my_palette, breaks = bk,border_color = NA, show_colnames = FALSE,
             annotation_row = annot_col2,
             annotation_col = annot_col1, annotation_colors = ann_colors,
             fontsize_row = 20, fontsize = 20)
    dev.off()
  }
   write.csv(t(abun.ss), file = "/mnt/zhao_trial/test_lin/new_heatmap.csv")
  
  ###########---------------------- clustering for phenotypes --------------------###########
  for(clust_phe in phenos){
    ###prepare data
    ss <- rownames(meta)[which(meta[,clinic_col] == clust_phe)]
    ss.abun <- abun.ss[ss,]
    ###set colors and breaks for heatmap
    bk1 <- c(seq(0,0.01,by=0.0001))
    bk2 <- c(seq(0.011,max(ss.abun),by = 0.01))
    bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing
    my_palette1 <- c(colorRampPalette(colors = c("white", "gold"))(n = length(bk1)))   #c(colorRampPalette(colors = c("darkred", "tomato1"))(n = length(bk2)-1))
    my_palette2 <- c(colorRampPalette(colors = c("gold", "red"))(n = length(bk2)))   #c(colorRampPalette(colors = c("darkred", "tomato1"))(n = length(bk2)-1))
    my_palette <- c(my_palette1,my_palette2)
    ###heatmap
    Pheight <- length(micros.sel)*30
    Pwidth <- nrow(ss.abun)*5+1200
    png(paste(outdir,clinic_col, "_", clust_phe, "_microbes_clustering.png", sep = ""), width = Pwidth, height = Pheight)
    pheatmap(t(ss.abun), color = my_palette, breaks = bk,border_color = NA, show_colnames = FALSE, 
             fontsize_row = 20, fontsize = 25, main = paste("Microbiota Abundace in ",clust_phe, sep = "") )
    dev.off()
  }}
