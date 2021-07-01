suppressMessages(library(optparse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

###make option list and parse command line
option_list <- list(
  make_option(c("-i","--msiscore"), type = "character",
              help = "msi score"),
  make_option(c("-j","--tidescore"), type = "character",
              help = "tide score"),
  make_option(c("-m", "--meta"), type = "character",
              help = "meta data path"),
  make_option(c("-o","--outdir",type="character", 
              help="Output files [Required]")),
  make_option(c("-g", "--condition"), type="character", 
              help="Treatment", metavar="character"),
  make_option(c("-t", "--treatment"), type="character", 
              help="Treatment", metavar="character"),
  make_option(c("-c", "--control"), type="character", 
              help="Control", metavar="character")
)


opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

outdir<-opt$outdir
Treatment <- opt$treatment
Control <- opt$control
Condition <- opt$condition
meta <- opt$meta


###set parameters
meta <- read.table(file = meta, sep=',', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
msi <- read.table(opt$msiscore,stringsAsFactors = FALSE, header = TRUE,sep ='\t',row.names = 1)
tide <- read.table(opt$tidescore,stringsAsFactors = FALSE, header = TRUE,sep ='\t',row.names = 1)
tmp <- merge(tide,msi,by=0)
dat <- cbind(meta[tmp$Row.names,Condition],tmp)
dat$clinic <- dat[,1]

CompareGroups <- function(dat,col,metric){
 
 p.value <- wilcox.test(dat[,col] ~ clinic, dat)$p.value
 p.value <- sprintf(p.value, fmt = '%#.3f')
 par(oma=c(2,2,2,2))
 y_pos <- max(dat[,col])*1.2
 gp <- ggplot(dat, aes_string(x = factor(dat$clinic, levels=c(Treatment,Control)), y = col, fill = 'clinic'))+
    geom_boxplot(lwd=0.3,size=0.3,outlier.size=-1,width=0.3,alpha=0.3,position = position_dodge2(preserve = "single"))+
    geom_jitter(shape=16, position=position_dodge(0.5)) +
    theme_bw() +
    ylab(metric)+
    ggtitle(paste('wilcoxon p-value = ',p.value,sep=''))+
    theme(plot.title = element_text(size=8,hjust=0.5,vjust = 0.5, margin = margin(l=10,r=5,t=4,b=4),face = "bold", colour = "black"),
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

g_msi <- CompareGroups(dat,"MSI_score","MSI_score")
g_tide <- CompareGroups(dat,"TIDE","TIDE")
g_Dys <- CompareGroups(dat,"Dysfunction","Dysfunction")
g_Exl <- CompareGroups(dat,"Exclusion","Exclusion")
g_mdsc <- CompareGroups(dat,"MDSC","MDSC")
g_TAM <- CompareGroups(dat,"TAM.M2","TAM.M2")
g_CAF <- CompareGroups(dat,"CAF","CAF")
g_IFNG <- CompareGroups(dat,"IFNG","IFNG")

png(file = paste(outdir,Condition,"_comparison.png",sep = ""),res = 300, width = 2200, height = 1150)
ggarrange(g_msi, g_tide, g_Dys, g_Exl,g_mdsc,g_TAM,g_CAF,g_IFNG,ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()


