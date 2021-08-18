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
  make_option(c("-o","--outdir",type="character", help="Output files [Required]"))
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###set parameters
meta<-opt$meta
inputdir<-opt$input_path
outdir<-opt$outdir
clinic.col<-opt$clinic_col


# prepare for plot 
meta <- read.csv(file = meta,header = TRUE,  row.names = 1)
cols <- length(unique(meta[,clinic.col]))
source("src/trust4_metric_functions.R")



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
  gp <- ggplot(dat, aes_string(x = F1, y = F2, fill = F1))+
    geom_violin(trim = TRUE,alpha = 0.6)+
    geom_boxplot(fill = "white")+
    theme_bw() +
    ylab(metric)+
    scale_fill_manual(values = brewer.pal(n = 8, name = "Set2")[1:cols])+
    stat_compare_means( aes(label = ..p.signif..), 
                        label.x = 1.5, label.y = y_pos,
                        tip.length = 0.01,size = 5)+
    theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
          axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
          axis.text.y=element_text(size=12,face = "bold",hjust=1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none")
  return(gp)
}

##############-------------Metric: fraction of BCR reads------------------------###############
################################################################################################
cdr3.lib <- bcr.lib.reads %>% mutate(clinic = as.character(meta[sample,clinic.col]))
y_pos <- max(cdr3.lib$Infil)*1.2
gr <- CompareGroups(cdr3.lib,"clinic","Infil",cols,y_pos,"Fraction of BCR reads")
gr

##############-------------Metric: unique BCR cdr3s size ------------------------###############
################################################################################################
cdr3.size <- aggregate(CDR3aa ~ sample+clinic, cdr3.bcr.heavy, function(x) length(unique(x))) 
y_pos <- max(cdr3.size$CDR3aa)*1.1
gc <- CompareGroups(cdr3.size,"clinic","CDR3aa",cols,y_pos,"Unique CDR3")
gc

##############-------------Metric: somatic hypermutation rate ------------------------###########
################################################################################################
all.SHM <- na.omit(cbind.data.frame(sample = names(SHMRatio), SHM.rate = unlist(SHMRatio))) %>% 
  mutate(clinic = meta[sample,clinic.col])
y_pos <- max(all.SHM$SHM.rate)*1.1
gs <- CompareGroups(all.SHM,"clinic","SHM.rate",cols,y_pos,"SHM rate")
gs


##############-------------Metric: BCR clonality ------------------------###########################
################################################################################################
clonality.all <- data.frame(do.call("rbind",bcr_clonality)) %>% 
  mutate(clonality = as.numeric(as.character(clonality))) %>%
  mutate(clinic = meta[sample,clinic.col])
y_pos <- max(clonality.all$clonality)*1.5
gn <- CompareGroups(clonality.all,"clinic","clonality",cols,y_pos,"Clonality")
gn

###combine BCR Metrics Plots
png(file = paste(outdir,"TRUST4_BCR_heavy_metrics_plot.png",sep = ""),res = 300, width = 3000, height = 900)
ggarrange(gr, gc, gs, gn,ncol = 4, nrow = 1)
dev.off()

##############-------------Metric: Ig fraction ------------------------###########################
################################################################################################
Igs.order <- c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2")
Ig.color <- c("gold","seagreen","hotpink","brown1","skyblue","firebrick","deeppink","purple","dodgerblue")
names(Ig.color) <- Igs.order
st.Ig <- cdr3.bcr.heavy %>% 
  group_by(clinic,sample) %>%
  mutate(est_clonal_exp_norm = frequency/sum(frequency)) %>%   #as.numeric(sample.clones[filename,2])
  dplyr::filter(C != "*") %>%
  group_by(clinic, C) %>% 
  dplyr::summarise(Num.Ig = sum(est_clonal_exp_norm))
st.Ig$C <- factor(st.Ig$C, levels = Igs.order)
png(file = paste(outdir,"TRUST4_BCR_Ig_frequency.png",sep = ""),res = 300, width = 1500, height = 1500)
ggplot(st.Ig,aes(x = clinic, y = Num.Ig, fill = C))+
  geom_bar(stat = "identity",position="fill",alpha = 0.8)+
  theme_bw()+
  labs(y = "Normalized Ig Abundance",fill = "Ig")+
  scale_fill_manual(values = Ig.color)+
  theme(legend.position = "right",
        axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"))
dev.off()


##############-------------Metric: Ig class switch ------------------------###########################
################################################################################################
bcr_clusters.count <- sum(unlist(lapply(bcr_clusters, length)))
network.matrix <- matrix(0, 9, 9)
rownames(network.matrix) = colnames(network.matrix) = c('IGHA1','IGHA2','IGHE','IGHD','IGHM','IGHG1','IGHG2','IGHG3','IGHG4')
for(i in 1:nrow(bcr.cluster.cs)){
  id=colnames(bcr.cluster.cs)[which(as.numeric(bcr.cluster.cs[i, 3:11])>0)+2]
  network.matrix[id,id]=network.matrix[id,id]+1
}
for(i in 1:9){
  network.matrix[i,i]=0
}
network.edge=melt(network.matrix)
network.edge = network.edge[-which(network.edge$value==0), ]

network.data <- network(network.matrix, directed = TRUE, matrix.type='adjacency')
set.edge.attribute(network.data, attrname = 'as', value = network.edge[, 'value']*50/bcr_clusters.count)
node.attr=rowSums(network.matrix)/sum(network.matrix)
node.attr[which(node.attr==0)]=0.001
set.vertex.attribute(network.data, attrname = 'node', value = node.attr)
network.data %v% "Ig" = Igs.order
png(file = paste(outdir,"Ig_CS_network.png",sep = ""), res = 300, width = 1800, height = 1500)
ggnet2(network.data, mode = 'circle',  color = "Ig", node.alpha = 1,
       label = FALSE, label.color = "black", palette = "Paired", 
       label.size = 4, edge.color = 'grey80', edge.alpha = 0.5) 
#edge.size = 'as'
dev.off()

##############-------------Metric: BCR jacard similarity ------------------------###########################
################################################################################################
cdr3.bcr.complete <- subset(cdr3.bcr.heavy, is_complete == "Y")
share.mat <- getCDR3Jacard(meta, cdr3.bcr.complete)
share.mat$s1 <- factor(share.mat$s1, levels = rownames(meta))
share.mat$s2 <- factor(share.mat$s2, levels = rownames(meta))
sampleN <- length(unique(share.mat[,1]))
Pwidth <- sampleN*54+1080
Pheight <- sampleN*48+960
png(file = paste(outdir,"TRUST4_BCR_Jacard_plot.png", sep = ""), res = 300, width = Pwidth, height = Pheight)
ggplot(data = share.mat, aes(x = s1, y = s2, fill = jacard))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = max(share.mat$jacard)/2, limit = c(0, max(share.mat$jacard)), space = "Lab", 
                       name = "Jacard")+
  theme_minimal()+
  ggtitle("BCR Jacard Similarity")+
  theme(legend.position = "right",
        plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x = element_text(size = 12,face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12,face = "bold", hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()



##############-------------Metric: fraction of TCR reads ------------------------###########################
################################################################################################
cdr3.lib <- tcr.lib.reads %>% mutate(clinic = meta[sample,clinic.col])
y_pos <- max(cdr3.lib$Infil)*1.2
gr <- CompareGroups(cdr3.lib,"clinic","Infil",cols,y_pos,"Fraction of TCR reads")
gr

##############-------------Metric: TCR unique cdr3 size ------------------------###########################
################################################################################################
cdr3.size <- aggregate(CDR3aa ~ sample+clinic, cdr3.tcr, function(x) length(unique(x))) 
y_pos <- max(cdr3.size$CDR3aa)*1.1
gc <- CompareGroups(cdr3.size,"clinic","CDR3aa",cols,y_pos,"Unique CDR3")
gc

##############-------------Metric: TCR CPK ------------------------###########################
################################################################################################
cpk <- aggregate(CDR3aa ~ sample+clinic+lib.size, cdr3.tcr, function(x) length(unique(x))) %>%
  mutate(CPK = signif(CDR3aa/(lib.size/1000),4))
y_pos <- max(cpk$CPK)*1.2
gk <- CompareGroups(cpk,"clinic","CPK",cols,y_pos,"Clonetypes per kilo reads")
gk

##############-------------Metric: TCR clonality ------------------------###########################
################################################################################################
clonality.all <- data.frame(do.call("rbind",tcr_clonality)) %>% 
  mutate(clonality = as.numeric(as.character(clonality))) %>%
  mutate(clinic = meta[sample,clinic.col])
y_pos <- max(clonality.all$clonality)*1.1
gn <- CompareGroups(clonality.all,"clinic","clonality",cols,y_pos,"Clonality")
gn

###combine TCR Metrics Plots
png(file = paste(outdir,"TRUST4_TCR_metrics_plot.png",sep = ""),res = 300, width = 3000, height = 900)
ggarrange(gr, gc, gk, gn,ncol = 4, nrow = 1)
dev.off()




##############-------------Metric: TCR jacard similarity ------------------------###########################
################################################################################################
cdr3.tcr.complete <- subset(cdr3.tcr, is_complete == "Y")
share.mat <- getCDR3Jacard(meta, cdr3.tcr.complete)
share.mat$s1 <- factor(share.mat$s1, levels = rownames(meta))
share.mat$s2 <- factor(share.mat$s2, levels = rownames(meta))
sampleN <- length(unique(share.mat[,1]))
Pwidth <- sampleN*54+1080
Pheight <- sampleN*48+960
png(file = paste(outdir,"TRUST4_TCR_Jacard_plot.png", sep = ""), res = 300, width = Pwidth, height = Pheight)
ggplot(data = share.mat, aes(x = s1, y = s2, fill = jacard))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = max(share.mat$jacard)/2, limit = c(0, max(share.mat$jacard)), space = "Lab", 
                       name = "Jacard")+
  theme_minimal()+
  ggtitle("TCR Jacard Similarity")+
  theme(legend.position = "right",
        plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x = element_text(size = 12,face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12,face = "bold", hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()




