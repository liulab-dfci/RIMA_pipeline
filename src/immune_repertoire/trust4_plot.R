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
  make_option(c("-i", "--input_path"), type="character", 
              help="Input path of processed cdr3 data. [Required]"),
  make_option(c("-c", "--clinic_col"), type="character",
              help="column number of clinic traits in meta file[Required]"),
  make_option(c("-m", "--meta"), type="character",
              help="meta info[Required]"),
  make_option(c("-o","--outdir",type="character", help="Output files [Required]")),
  make_option(c("-s", "--order"), type="character", default=NULL,
              help="sample order for multiqc report", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###set parameters
meta<-opt$meta
inputdir<-opt$input_path
outdir<-opt$outdir
clinic.col<-opt$clinic_col
col <- opt$order

# prepare for plot 
meta <- read.csv(file = meta,header = TRUE)
cols <- length(unique(meta[,clinic.col]))
phenotype<- meta[,clinic.col]
source("src/immune_repertoire/trust4_metric_functions.R")

if (clinic.col == "Timing"){
  meta$Timing <- ifelse(meta$Timing == "pre"|meta$Timing == "Pre", "Apre",
                        ifelse(meta$Timing == "post"|meta$Timing == "Post", "Post", "other"))
}



##############-------------load data for BCR heavy chain  ------------------------###############
load(paste(inputdir,"TRUST4_BCR_heavy_cluster.Rdata", sep = ""))
load(paste(inputdir,"TRUST4_BCR_heavy_clonality.Rdata", sep = ""))
load(paste(inputdir,"TRUST4_BCR_heavy_SHMRatio.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_BCR_heavy_lib_reads_Infil.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_BCR_Ig_CS.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_BCR_heavy_lib_reads_Infil.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_BCR_heavy.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_TCR_lib_reads_Infil.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_TCR.Rdata",sep = ""))
load(paste(inputdir,"TRUST4_TCR_clonality.Rdata",sep = ""))

###function of comparing BCR metrics in the multiple groups
CompareGroups <- function(dat,F1,F2,cols,y_pos,metric){
 par(oma=c(2,2,2,2))
 gp <- ggplot(dat, aes_string(x = F1, y = F2, fill = F1))+
    #geom_violin(trim = TRUE,alpha = 0.3,scale = "width",lwd=0.2,width=0.3)+
    geom_boxplot(lwd=0.3,size=0.3,outlier.size=-1,width=0.3,alpha=0.3,position = position_dodge2(preserve = "single"))+
    theme_bw() +
    ylab(metric)+
    scale_fill_manual(values = brewer.pal(n = 8, name = "Set2")[1:cols])+
    stat_compare_means( aes(label = ..p.format..),
                        label.x = 1.4, label.y = y_pos,
                        tip.length = 0.01,size = 2)+
    scale_x_discrete(name ="", labels=c("Apre" = "Pre")) +

    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=10,r=5,t=4,b=4),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=4,face = "bold",hjust=1),
          axis.text.y=element_text(size=6,hjust=1,face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 6,face="bold"),
          legend.position = "bottom",
	  legend.title = element_text(size=6, face = "bold"),
	  legend.text = element_text(size=6, face = "bold")) + guides(fill=FALSE)
    theme(plot.margin = unit(c(0,0.05,0,0.05), "cm"))
 return(gp)
}


##############-------------Metric: fraction of BCR reads------------------------###############
################################################################################################
cdr3.lib <- merge(bcr.lib.reads,meta,by.x='sample',by.y='SampleName')
y_pos <- max(cdr3.lib$Infil)*1.2
#y_pos<- 0.0001
gr <- CompareGroups(cdr3.lib,clinic.col,"Infil",cols,y_pos,"Fraction of BCR reads")
if (nrow(bcr.lib.reads) == 0) {
  gr <- NULL
}

if(clinic.col=="Timing") {

gr <- gr + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) + 
	geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}
gr<-gr + ylim(0,y_pos)


##############-------------Metric: unique BCR cdr3s size ------------------------###############
################################################################################################

tmp <- aggregate(CDR3aa ~ sample, cdr3.bcr.heavy, function(x) length(unique(x))) 
cdr3.size <- merge(tmp,meta,by.x='sample',by.y='SampleName')
#y_pos=100
y_pos <- max(cdr3.size$CDR3aa)*1.2
gc <- CompareGroups(cdr3.size,clinic.col,"CDR3aa",cols,y_pos,"Unique CDR3")

if(clinic.col=="Timing") {
gc <- gc + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}
gc<- gc + ylim(0,y_pos)


##############-------------Metric: somatic hypermutation rate ------------------------###########
################################################################################################
tmp <- na.omit(cbind.data.frame(sample = names(SHMRatio), SHM.rate = unlist(SHMRatio)))
tmp['name'] <-   str_replace(tmp$sample, "_TRUST4_BCR_heavy_cluster.Rdata", "")
  
all.SHM <- merge(tmp,meta,by.x='name',by.y='SampleName')

y_pos <- max(all.SHM$SHM.rate)*1.0
gs <- CompareGroups(all.SHM,clinic.col,"SHM.rate",cols,y_pos,"SHM rate")

if(clinic.col=="Timing") {
gs <- gs + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}
gs<- gs+ ylim(0,y_pos)



##############-------------Metric: BCR clonality ------------------------###########################
################################################################################################
tmp <- data.frame(do.call("rbind",bcr_clonality)) %>% 
  mutate(clonality = as.numeric(as.character(clonality)))

clonality.all <- merge(tmp,meta,by.x='sample',by.y='SampleName')

if ("NaN"%in%clonality.all$clonality) {
	clonality.all <- subset(clonality.all, clonality.all[[2]] != "NaN")

}


y_pos <- max(clonality.all$clonality)*1.0
gn <- CompareGroups(clonality.all,clinic.col,"clonality",cols,y_pos,"Clonality")
#gn

if(clinic.col=="Timing") {
gn <- gn + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size = 0.3)
}
gn<- gn + ylim(0,y_pos)
###combine BCR Metrics Plots
png(file = paste(outdir,"TRUST4-BCR_mqc.png",sep = ""),res = 300, width = 2200, height = 550)
ggarrange(gr, gc, gs, gn, ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()



#############-------------Metric: Ig fraction ------------------------###########################
################################################################################################
Igs.order <- c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2")
Ig.color <- c("gold","seagreen","hotpink","brown1","skyblue","firebrick","deeppink","purple","dodgerblue")
names(Ig.color) <- Igs.order

tmp <- merge(cdr3.bcr.heavy,meta,by.x='sample',by.y='SampleName') 
tmp['clinic'] <- tmp[,clinic.col]
cdr3.bcr.heavy_tmp <- tmp

st.Ig <- cdr3.bcr.heavy_tmp %>% 
  group_by(clinic,sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*") %>%
  group_by(clinic, C) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))

write.table(st.Ig,file=paste(outdir,"TRUST_Ig.txt",sep = ""),quote=FALSE, row.names=FALSE,sep="\t")

st.Ig$C <- factor(st.Ig$C, levels = Igs.order)
png(file = paste(outdir,"TRUST4_BCR_Ig_frequency.png",sep = ""),res = 300, width = 1000, height = 800, pointsize = 4)
ggplot(st.Ig,aes(x = clinic, y = Num.Ig, fill = C))+
  geom_bar(stat = "identity",position="fill",alpha = 0.8)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values = Ig.color)+
  theme(legend.position = "right",
        axis.text.x=element_text(angle=70,size=6,face = "bold",hjust=1),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6,face = "bold"))
dev.off()





##############-------------Metric: fraction of TCR reads ------------------------###########################
################################################################################################
cdr3.lib <- tcr.lib.reads %>% mutate(clinic = meta[sample,clinic.col])
y_pos <- max(cdr3.lib$Infil)+0.0000001

if (clinic.col == "Timing") {
  cdr3.lib <- merge(cdr3.lib, meta[c("SampleName", "Responder", "PatName")], by = 1, all = FALSE)
}
gr <- CompareGroups(cdr3.lib,"clinic","Infil",cols,y_pos,"Fraction of TCR reads")
if (nrow(tcr.lib.reads) == 0) {
  gr <- NULL
}

if(clinic.col=="Timing") {
gr <- gr + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}

gr<- gr + ylim(0, y_pos)


##############-------------Metric: TCR unique cdr3 size ------------------------###########################
################################################################################################
tmp <- aggregate(CDR3aa ~ sample, cdr3.tcr, function(x) length(unique(x))) 
cdr3.size <- merge(tmp,meta,by.x='sample',by.y='SampleName') 

y_pos <- max(cdr3.size$CDR3aa)-0.4
gc <- CompareGroups(cdr3.size,clinic.col,"CDR3aa",cols,y_pos,"Unique CDR3")

if(clinic.col=="Timing") {
gc <- gc + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}

gc <- gc + ylim(0,y_pos)


##############-------------Metric: TCR CPK ------------------------###########################
################################################################################################
tmp <- aggregate(CDR3aa ~ sample+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
  
cpk <-  merge(tmp,meta,by.x='sample',by.y='SampleName') 
y_pos <- max(cpk$CPK)-0.4

gk <- CompareGroups(cpk,clinic.col,"CPK",cols,y_pos,"Clonetypes per kilo reads")

if(clinic.col=="Timing") {
gk <- gk + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}

gk<- gk + ylim(0,y_pos)


##############-------------Metric: TCR clonality ------------------------###########################
################################################################################################
tmp <- data.frame(do.call("rbind",tcr_clonality)) %>% 
  mutate(clonality = as.numeric(as.character(clonality)))

clonality.all <- merge(tmp,meta,by.x='sample',by.y='SampleName')
if ("NaN"%in%clonality.all$clonality) {
  clonality.all <- subset(clonality.all, clonality.all[[2]] != "NaN")
}
y_pos <- max(clonality.all$clonality)-0.4

gn <- CompareGroups(clonality.all,clinic.col,"clonality",cols,y_pos,"Clonality")

if(clinic.col=="Timing") {
gn <- gn + geom_point(aes(shape = Responder, color = Responder),alpha=0.7, size = 1) +
        geom_line(aes(group = PatName, color = Responder), alpha = 0.5, size =0.3)
}

gn <- gn + ylim(0,y_pos)

###combine TCR Metrics Plots
png(file = paste(outdir,"TRUST4-TCR_mqc.png",sep = ""),res = 300, width = 2200, height = 550)
par(oma=c(1,1,1,1))
ggarrange(gr, gc, gk, gn,ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
dev.off()

#multiqc report

#extract the ig info by samples
Ig_sample <- cdr3.bcr.heavy_tmp %>%
  group_by(sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*") %>%
  group_by(sample, C) %>%
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))

#create multiqc matrix 
ind_matrix = matrix(rep(FALSE, length(meta$SampleName) * length(Igs.order)), nrow = length(meta$SampleName))
rownames(ind_matrix) <- meta$SampleName
colnames(ind_matrix) <- Igs.order
for (i in Igs.order) {
  ig_tmp <- subset(Ig_sample, C == i)
  test <- merge(ig_tmp, as.data.frame(meta$SampleName), by = 1, all = TRUE)
  rownames(test) <- test$sample
  test <- test[rownames(ind_matrix),]
  ind_matrix[,i] <- test$Num.Ig
}

ind_matrix <- cbind(rownames(ind_matrix), ind_matrix)
ind_matrix[is.na(ind_matrix)] <- 0
colnames(ind_matrix)[1] <- "samples"

write.table(ind_matrix, file = paste(outdir,"TRUST_Ig.txt",sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
