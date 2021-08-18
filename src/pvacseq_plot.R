############################################################################
#pre-preparation
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
library(optparse)

###options
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input expression file", metavar="character"),
  make_option(c("-out", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###parameters setting
input <- opt$input
outdir <- opt$outdir

# ### testing
# input = "~/Documents/pipeline_test/kraken_test/Merged.filtered.condensed.ranked.addSample.tsv"
# outdir <- "~/Documents/rnaseq/data/"

##############################################################################
### Reading matrix
input_table = read.table(input, header=TRUE,sep = '\t')
input_table$is_novel = ifelse(input_table$Median.Fold.Change>1 & input_table$Median.MT.Score <500, "novel","no-novel")
COLORS<- brewer.pal(n = 3, name = "Set1")

#########---------------------count the number of epitopes in each patient----------------###########
pat.epit <- input_table %>% group_by(Sample,is_novel) %>% dplyr::summarise(count = n())
Pwidth <- 800+80*length(unique(pat.epit[,1]))
png(file = paste(outdir, "Patient_count_epitopes_plot.png",sep = ""), res = 300, width = Pwidth, height = 1000)
ggplot(pat.epit, aes(x = Sample,y=count,fill = is_novel)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  theme_minimal()+
  scale_fill_manual(values = brewer.pal(n=8, name = "Set1")[2:1])+
  ylab("#Neoepitopes")+xlab("Patient")+
  theme(legend.position = "top",
        plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(size=12,face = "bold",hjust=1),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"))
dev.off()
 
##########--------------------predicted epitope Affinity scatter plot------------------##############
png(file = paste(outdir, "epitopes_affinity_plot.png",sep = ""), res = 300, width = 1500, height = 1200)
ggplot(input_table, aes(y = Median.MT.Score, x = Median.WT.Score,colour = is_novel)) +
  geom_point()+
  xlab("Normal epitope affinity(nM)")+ylab("Tumor epitope affinity(nM)")+
  theme_minimal()+
  scale_x_log10() +
  scale_y_log10()+
  geom_hline(aes(yintercept=500)) +
  geom_vline(aes(xintercept=500))+
  geom_abline(intercept=0,slope=1 ,linetype='dashed')+ 
  annotate(geom="text", x=250, y=5, label="x=500")+
  annotate(geom="text", x=50, y=450, label="y=500")+
  theme(plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(size=12,face = "bold",hjust=1),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"))
dev.off()

##########--------------------the HLA allele fraction for patients(how many patients have this HLA?)------------------##############
hla.fraction <- input_table %>% subset(select=c(Sample,HLA.Allele,is_novel)) %>% distinct() %>%
  group_by(HLA.Allele, is_novel) %>% dplyr::summarise(count = n()) %>% 
  group_by(HLA.Allele, is_novel) %>% mutate(fraction = count/sum(count)) %>%
  mutate(HLA_type = substring(HLA.Allele, 5,5))
hla.fraction$HLA.Allele <- factor(hla.fraction$HLA.Allele,levels = as.character(unique(hla.fraction$HLA.Allele)))
g1 <- ggplot(hla.fraction, aes(x = HLA.Allele, y = fraction, fill = HLA_type))+
  geom_bar(aes(alpha = is_novel),stat="identity", position=position_dodge())+
  scale_fill_manual(values = COLORS)+
  scale_alpha_manual(values = c('novel' = 0.9,'no-novel' = 0.4), name = "Type")+
  scale_y_reverse()+
  theme_minimal()+
  ylab("Patient Fraction")+xlab("")+
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.3,0.3,0.4),"cm"),
        plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x=element_text(angle=70,size=12,face = "bold",hjust=1),
        axis.text.y=element_text(size=12,face = "bold",hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"))

#########-------------the number of epitopes for each HLA(each HLA can be binded with how many epitopes)------------------##############
hla.epit <- input_table %>% 
  group_by(HLA.Allele,is_novel) %>% 
  dplyr::summarise(count = n()) %>% 
  mutate(HLA_type = substring(HLA.Allele, 5,5))
hla.epit$HLA.Allele <- factor(hla.epit$HLA.Allele, levels = as.character(unique(hla.fraction$HLA.Allele)))
g2 <- ggplot(hla.epit, aes(x = HLA.Allele, y = count, fill = HLA_type))+
  geom_bar(aes(alpha = is_novel),stat="identity", position=position_dodge())+
  scale_fill_manual(values = COLORS)+
  scale_alpha_manual(values = c('novel' = 0.9,'no-novel' = 0.4), name = "Type")+
  theme_minimal()+
  ylab("#Neoepitopes")+xlab("")+
  theme(legend.position = "top")+
  theme(plot.margin = unit(c(0.5,0.3,0.3,0.65),"cm"),
        plot.title = element_text(hjust=0.5,vjust = 0.5, margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"),
        axis.text.x= element_blank(),
        axis.text.y=element_text(size=12,face = "bold",hjust=1),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"))

###merge g1 and g2
Pwidth <- length(unique(hla.fraction[[1]]))*60+1440
png(file = paste(outdir, "HLA_epitopes_fraction_plot.png",sep = ""), res = 300, width = Pwidth, height = 2000)
ggarrange(g2,g1,ncol=1,nrow=2)
dev.off()



 




