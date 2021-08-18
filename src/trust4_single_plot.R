library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(optparse)
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="trust4 output", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


file <- opt$input
outdir <- opt$outdir
# outdir <- "~/Desktop/"
# file <- "~/Desktop/temp/plot/test/analysis/trust4/SRR3184301/SRR3184301_cdr3.out.processed.txt"
file_in<-read.table(file,header = TRUE,stringsAsFactors = FALSE)
chains<-c("TRA","TRB","IGL","IGH","IGK")
sample_name <- unlist(strsplit(file,'/'))[3]
for (chain in chains){
  print(chain)
  # chain="IGK"
  ### V gene frequency
  chain_V<-file_in[which(str_detect(file_in[,"V"],chain)),]
  chain_V<-str_extract(chain_V[,"V"], "[A-Z]+[0-9]{0,}")
  if (length(chain_V)!=0){
    summary<-as.data.frame(prop.table(table(chain_V)))
    colnames(summary)<-c("Vgene","Frequency")
    summary<-summary[order(summary[,2],decreasing = TRUE),]
    summary <- within(summary, Vgene <- factor(Vgene , levels=Vgene))
    Pwidth=(60*nrow(summary))+140
    if (Pwidth<500){
      Pwidth=500
    }
    png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_Vgene.png"),res = 200, width = Pwidth, height = 1000)
    print(ggplot(data=summary,aes(x=Vgene,y=Frequency))+geom_bar(stat="identity",fill="#6ecdb7",colour="black",position=position_dodge(0.5),width=0.5)+theme_bw()+theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),axis.text.x=element_text(angle = -90,size = 12,hjust = 1,face = "bold"),axis.text.y = element_text(size = 12,face = "bold"),axis.title.y = element_text(size = 12,face = "bold"))+labs(title=paste0(chain," Vgene usage"),x=""))
    dev.off()
  }
  
  ### J gene frequency
  chain_J<-file_in[which(str_detect(file_in[,"J"],chain)),]
  chain_J<-str_extract(chain_J[,"J"], "[A-Z]+[0-9]{0,}")
  if (length(chain_J)!=0){
    summary<-as.data.frame(prop.table(table(chain_J)))
    colnames(summary)<-c("Jgene","Frequency")
    summary<-summary[order(summary[,2],decreasing = TRUE),]
    summary <- within(summary, Jgene <- factor(Jgene , levels=Jgene))
    Pwidth=(60*nrow(summary))+140
    if (Pwidth<500){
      Pwidth=500
    }
    png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_Jgene.png"),res = 200, width = Pwidth, height = 1000)
    print(ggplot(data=summary,aes(x=Jgene,y=Frequency))+geom_bar(stat="identity",fill="#feb24c",colour="black",position=position_dodge(0.5),width=0.5)+theme_bw()+theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),axis.text.x=element_text(angle = -90,size = 12,hjust = 1,face = "bold"),axis.text.y = element_text(size = 12,face = "bold"),axis.title.y = element_text(size = 12,face = "bold"))+labs(title=paste0(chain," Jgene usage"),x=""))
    dev.off()
  }
  
  ### barplot showing VJ pair count
  VJ<-file_in[which(str_detect(file_in[,"V"],"V")&str_detect(file_in[,"J"],"J")),]
  VJ[,"V"]<-str_extract(VJ[,"V"], "[A-Z]+[0-9]{0,}")
  VJ[,"J"]<-str_extract(VJ[,"J"], "[A-Z]+[0-9]{0,}")
   VJ<-VJ[which(str_detect(VJ[,"V"],chain)&str_detect(VJ[,"J"],chain)),]
  if (nrow(VJ)!=0){
    VJ<-VJ[,c("frequency","V","J")]
    colnames(VJ)<-c("Frequency","Vgene","Jgene")
    VJ <- VJ[order(VJ$Vgene,VJ$Jgene),]
    VJ[,1] <- as.numeric(VJ[,1]) 
    if (sum(VJ[,1])!=0){
      VJ[,1]=VJ[,1]/sum(VJ[,1])
      colourCount = length(unique(VJ$Jgene))
      getPalette = colorRampPalette(brewer.pal(9, "Set3"))
      Pwidth=(75*length(unique(VJ[,2])))+140
      if (Pwidth<500){
        Pwidth=500
      }
      png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_VJpair.png"),res = 200, width = Pwidth, height = 1000)
      print(ggplot(VJ,aes(x=Vgene,y=Frequency,fill=Jgene))+geom_bar(stat = 'identity', width = 0.5, position = 'stack')+theme_bw()+theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin =margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),axis.text.x =  element_text(angle = -90,size = 12,hjust = 1,face = "bold"),axis.text.y = element_text(size = 12,face = "bold"),axis.title.y = element_text(size = 12,face = "bold"),legend.key.width=unit(0.08,"inches"),legend.key.height=unit(0.08,"inches"))+scale_fill_manual(values = getPalette(colourCount))+labs(title = paste0(chain," VJgene usage"),x=""))
      dev.off()
      ### heatmap showing VJ pair frequency
      Vgene<-unique(VJ[,2])
      Jgene<-unique(VJ[,3])
      htMatrix<-matrix(nrow = length(Vgene),ncol = length(Jgene))
      for (i in 1:length(Vgene)){
        for (j in 1:length(Jgene)){
          temp<-VJ[which(str_detect(VJ[,2],Vgene[i])&str_detect(VJ[,3],Jgene[j])),]
          htMatrix[i,j]=sum(temp[,1])
        }
      }
      htMatrix<-as.data.frame(htMatrix)
      colnames(htMatrix)<-Jgene
      rownames(htMatrix)<-Vgene
      Pwidth<-length(Jgene)*30+500
      Pheight<-length(Vgene)*30+500
      #treeheight_col = length(Vgene)
      #treeheight_row = length(Jgene)
      if(dim(htMatrix)[1] >= 2){cl_row <- TRUE} else{cl_row <- FALSE}
      if(dim(htMatrix)[2] >= 2){ cl_col <- TRUE} else{cl_col <- FALSE}
      png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_VJpair_ht.png"),res = 200, width = Pwidth, height = Pheight)
      print(pheatmap(htMatrix,main=paste0(chain," VJ pairing usage"), cluster_rows = cl_row, cluster_cols = cl_col))
      dev.off()
    } else {
      png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_VJpair.png"),res = 200, width = 1000, height = 1000)
      print(ggplot())
      dev.off()
      png(paste0(outdir,sample_name,"/",sample_name,"_",chain,"_VJpair_ht.png"),res = 200, width =1000, height = 1000)
      print(ggplot())
      dev.off()
    }
  }
} 
AA<-file_in[which(str_detect(file_in[,"CDR3aa"],"[A-Z]")),]
if (nrow(AA)!=0){
  lenAA_dis<-as.data.frame(prop.table(table(str_length(AA[,"CDR3aa"]))))
  colnames(lenAA_dis)<-c("cdr3Length","Frequency")
  Pwidth=(60*nrow(lenAA_dis))+140
  if (Pwidth<500){
    Pwidth=500
  }
  png(paste0(outdir,sample_name,"/",sample_name,"_cdr3_length.png"),res = 200, width = Pwidth, height = 1000)
  print(ggplot(data=lenAA_dis,aes(x=cdr3Length,y=Frequency))+geom_bar(stat = "identity",fill="#87ceeb",colour="black",position=position_dodge(0.5),width=0.5)+theme_bw()+theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),axis.text.x=element_text(size = 12,hjust = 1,face = "bold"),axis.text.y = element_text(size = 12,face = "bold"),axis.title.y = element_text(size = 12,face = "bold"),axis.title.x = element_text(size = 12,face = "bold"))+labs(title = "cdr3 length"))
  dev.off()
}