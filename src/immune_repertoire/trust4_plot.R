#!/usr/bin/env Rscript

## load packages
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(reshape))
suppressMessages(library(ggpubr))
suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(ggcorrplot))
suppressMessages(library(optparse))
suppressMessages(library(sna))

## make option list and parse command line
option_list <- list(  
  make_option(c("-i", "--infil_bcr"), type="character", 
              help="BCR fraction [Required]"),
  make_option(c("-j", "--heavy_bcr"), type="character", 
              help="BCR heavy chain[Required]"),
  make_option(c("-s", "--shm"), type="character", 
              help="somatic hypermutation [Required]"),
  make_option(c("-k", "--clone_bcr"), type="character", 
              help="BCR clonality [Required]"),
  
  make_option(c("-l", "--infil_tcr"), type="character", 
              help="TCR fraction [Required]"),
  make_option(c("-m", "--tcr"), type="character", 
              help="TCR  [Required]"),
  make_option(c("-p", "--clone_tcr"), type="character", 
              help="TCR clonality [Required]"),
  
  make_option(c("-o","--outdir",type="character", 
              help="Output files [Required]")),
  make_option(c("-g", "--condition"), type="character", 
              help="Treatment", metavar="character"),
  make_option(c("-t", "--treatment"), type="character", 
              help="Treatment", metavar="character"),
  make_option(c("-c", "--control"), type="character", 
              help="Control", metavar="character")

)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###set parameters
infil_bcr_input<- opt$infil_bcr
heavy_bcr_input <- opt$heavy_bcr
shm_input <- opt$shm
clone_bcr_input <- opt$clone_bcr


infil_tcr_input <- opt$infil_tcr
tcr_input <- opt$tcr
clone_tcr_input <- opt$clone_tcr
outdir<-opt$outdir
Treatment <- opt$treatment
Control <- opt$control
Condition <- opt$condition


###function of comparing BCR metrics in the multiple groups
CompareGroups <- function(dat,col,metric){
 
 p.value <- wilcox.test(dat[,col] ~ clinic, dat)$p.value
 p.value <- sprintf(p.value, fmt = '%#.3f')
 par(oma=c(2,2,2,2))
 y_pos <- max(dat[,col])*1.2
 gp <- ggplot(dat, aes_string(x = factor(dat$clinic, levels=c(Treatment,Control)), y = col, fill = 'clinic'))+
    geom_boxplot(lwd=0.3,size=0.3,outlier.size=-1,width=0.3,alpha=0.3,position = position_dodge2(preserve = "single"))+
    geom_jitter(shape=16, position=position_dodge(0.5), size = 0.6, alpha = 0.6) +
    theme_bw() +
    ylab(metric)+
    ggtitle(paste('wilcoxon p-value = ',p.value,sep=''))+
    theme(plot.title = element_text(size=5,hjust=0.5,vjust = 0.5, margin = margin(l=10,r=5,t=4,b=4),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=8,face = "bold",hjust=1),
          axis.text.y=element_text(size=8,hjust=1,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8,face="bold"),
          legend.position = "bottom",
	  legend.title = element_text(size=6, face = "bold"),
	  legend.text = element_text(size=6, face = "bold")) + guides(fill=FALSE)
    theme(plot.margin = unit(c(0,0.05,0,0.05), "cm"))
 return(gp)
}

##############-------------Metric: fraction of BCR reads------------------------###############
################################################################################################

infil_bcr <- read.csv(file = infil_bcr_input,header = TRUE,sep= '\t')
if (length(unique(infil_bcr$clinic)) < 2 ) {
  gr <- NULL
} else {
gr <- CompareGroups(infil_bcr,"Infil","Fraction of BCR reads")
}

##############-------------Metric: unique BCR cdr3s size and Ig frequency ------------------------###############
################################################################################################
cdr3.bcr.heavy <- read.csv(file = heavy_bcr_input, header = TRUE,sep= '\t')
if (length(unique(cdr3.bcr.heavy$clinic)) < 2) {
  gr <- NULL
  st.Ig <- data.frame(
    sample =character(),
    clinic = character(),
    C = character(),
    Num.Ig = charater()
  )
  write.table(st.Ig,file=paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4_Ig.txt",sep = ""),quote=FALSE, row.names=FALSE,sep="\t")
  
} else {
  dat <- aggregate(CDR3aa ~ sample + clinic, cdr3.bcr.heavy, function(x) length(unique(x))) 
  gc <- CompareGroups(dat,"CDR3aa","Unique CDR3")
  
  st.Ig <- cdr3.bcr.heavy%>% 
    group_by(sample,clinic) %>%
    mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%  
    dplyr::filter(C != ".") %>%
    group_by(sample,clinic, C) %>% 
    dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))
  write.table(st.Ig,file=paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4_Ig.txt",sep = ""),quote=FALSE, row.names=FALSE,sep="\t")
  
  ## Ig colors
  Igs.order <- c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2")
  Ig.color <- c("gold","seagreen","hotpink","brown1","skyblue","firebrick","deeppink","purple","dodgerblue")
  names(Ig.color) <- Igs.order
  st.Ig$C <- factor(st.Ig$C, levels = Igs.order)

  png(file = paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4_BCR_Ig_frequency.png",sep = ""),res = 300, width = 1000, height = 800, pointsize = 4)
  ggplot(st.Ig,aes(x = clinic, y = Num.Ig, fill = C))+
    geom_bar(stat = "identity",position="fill",alpha = 0.8)+
    theme_bw()+
    labs(y = "Normalized Ig Abundance",fill = "Ig")+
    scale_fill_manual(values = Ig.color)+
    theme(legend.position = "right",
          axis.text.x=element_text(angle=70,size=4,face = "bold",hjust=1),
          axis.text.y = element_text(size = 4,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 4,face = "bold")
          )
  dev.off()
  
}

##############-------------Metric: somatic hypermutation rate ------------------------###########
################################################################################################
shm <- read.csv(file = shm_input, header = TRUE,sep= '\t')
if (length(unique(shm$clinic)) < 2) {
  gs <- NULL
} else {
  gs <- CompareGroups(shm,"SHMRatio","SHM Ratio")
}


##############-------------Metric: BCR clonality ------------------------###########################
################################################################################################

clonality.bcr <- read.csv(file = clone_bcr_input, header = TRUE,sep= '\t')
if (length(unique(clonality.bcr$clinic)) < 2) {
  gn <- NULL
} else {
  gn <- CompareGroups(clonality.bcr,"clonality","Clonality")
}

###combine BCR Metrics Plots
png(file = paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4-BCR_mqc.png",sep = ""),res = 300, width = 2200, height = 550)
ggarrange(gr, gc, gs, gn, ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()




##############-------------Metric: fraction of TCR reads ------------------------###########################
################################################################################################
infil_tcr <- read.csv(file = infil_tcr_input,header = TRUE,sep= '\t')
if (length(unique(infil_tcr$clinic)) < 2 ) {
  gr <- NULL
} else {
  gr <- CompareGroups(infil_tcr,"Infil","Fraction of TCR reads")
}


##############-------------Metric: TCR unique cdr3 size and CPK ------------------------###########################
################################################################################################

cdr3.tcr <-  read.csv(file = tcr_input,header = TRUE,sep= '\t')
if (length(unique(cdr3.tcr$clinic)) < 2 ) {
  gc <- NULL
  gk <- NULL
} else {
  dat <- aggregate(CDR3aa ~ sample + clinic+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
    mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
  gc <- CompareGroups(dat,"CDR3aa","Unique CDR3")
  gk <- CompareGroups(dat,"CPK","Clonetypes per kilo reads")
}

##############-------------Metric: TCR clonality ------------------------###########################
################################################################################################
clonality.tcr <- read.csv(file = clone_tcr_input, header = TRUE,sep= '\t')
if (length(unique(clonality.tcr$clinic)) < 2) {
  gn <- NULL
} else {
  gn <- CompareGroups(clonality.tcr,"clonality","Clonality")
}

###combine TCR Metrics Plots
png(file = paste(outdir,Condition,"_",Treatment,"_vs_",Control,"_TRUST4-TCR_mqc.png",sep = ""),res = 300, width = 2200, height = 550)
par(oma=c(1,1,1,1))
ggarrange(gr, gc, gk, gn,ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()


