suppressMessages(library(optparse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

###make option list and parse command line
option_list <- list(
  make_option(c("-i","--msiscore"), type = "character",
              help = "msi score path"),
  make_option(c("-m", "--meta"), type = "character",
              help = "meta data path"),
  make_option(c("-p", "--phenotype"), type = "character",
              help = "phenotype information used to compare"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "output path")
)
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

###set parameters
meta<-read.csv(opt$meta,header = TRUE,stringsAsFactors = FALSE)
msi<-read.table(opt$msiscore,stringsAsFactors = FALSE)
phenotype<-unlist(strsplit(opt$phenotype,","))
outdir<-opt$outdir


###preprocess of the input files
colnames(msi)<-c("SampleName","msi_score")
library(ggpubr)
library(RColorBrewer)

###box plot for group comparison
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
msi_process<- function(phenotype){
  msi_phenotype<-merge(meta[,c(phenotype,"SampleName")],msi)
  colnames(msi_phenotype)<-c("SampleName","Phenotype","MSI_Score")
  #write.table(msi_phenotype, paste0(outdir,"msi_boxplot.txt"),sep="\t",quote=F,row.names=F)
  write.table(msi_phenotype, paste(opt$outdir,"msi_boxplot.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
  colourCount = length(unique(msi_phenotype[,"Phenotype"]))
  getPalette = colorRampPalette(brewer.pal(8, "Set3"))
  compare_list<-list(unique(msi_phenotype[,"Phenotype"]))

  p<-ggboxplot(msi_phenotype, x="Phenotype", y = "MSI_Score",color="Phenotype",add="jitter", size = 0.5) +
     scale_fill_manual(values = getPalette(colourCount))+
     theme_bw()+
     theme(axis.text.x=element_text(size=6,hjust=1),
        axis.text.y=element_text(size=6),
        axis.title.x = element_text(size = 8),
	legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.title.y = element_text(size = 8))+
     guides(fill=FALSE)+
     stat_compare_means(comparisons = compare_list,aes(label =paste0("p = ", ..p.format..)), size = 2.5)
     return(p)
}
plot_list<-lapply(phenotype, msi_process)
pwdith <- length(phenotype)*400
png(filename = paste0(outdir,"msi_score_comparison.png"),res = 300 ,width= 1300, height = 1050, pointsize = 1)
ggarrange(plotlist = plot_list)
dev.off()

####density plot for msi score of all samples
colnames(msi)<-c("sample","msi_score")
median<-median(msi[,2])
png(file = paste(outdir,"MSISensor.png",sep = ""), res = 300, width = 1500, height = 1500)
ggplot(msi, aes(x = msi_score))+
  theme_bw()+geom_density(color = "#6ecdb7")+
  geom_vline(xintercept = median,linetype="dashed", size=0.5)+
  geom_text(x=median,y=0.05,label = paste0("median=",median),size=7)+
  theme_bw()+
  labs(title = paste0("msi density plot"))+
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(size=12,face = "bold",hjust=1),
        axis.text.y=element_text(size=12,face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))
dev.off()

##testing

t <- ggplot(msi, aes(x = msi_score))+
  theme_bw()+geom_density(color = "#6ecdb7")+
  geom_vline(xintercept = median,linetype="dashed", size=0.5)+
  geom_text(x=median,y=0.05,label = paste0("median=",median),size=4)+
  theme_bw()+
  labs(title = paste0("MSI Density Plot"), x = "MSI_Score", y = "Density")+
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, size =10, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(size=6,hjust=1),
        axis.text.y=element_text(size=6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8))

p <- ggarrange(plotlist = plot_list)

png(file = paste(outdir,"Testing.png",sep = ""),res = 300, width = 2000, height = 900)
ggarrange(p, t, ncol = 2, nrow = 1)
dev.off()

save.image(paste(outdir,"Testing.Rdata",sep = ""))
