library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(reshape)
library(optparse)

# make option list and parse command line
option_list <- list(  
  make_option(c("-i", "--input"), type="character", 
              help="output from TIDE. [Required]"),
  make_option(c("-c", "--cc"), type="character",
              help="which cancer type to compare with in TCGA[Required]"),
  make_option(c("-e", "--expression"), type="character",
              help="normalized expression data"),
  make_option(c("-o","--outdir"),type="character", 
              help="Output files [Required]")
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


###set parameters
input <- opt$input
outdir <- opt$outdir
expr <- opt$expression
cc <- opt$cc


###read in user's data
tide.mat <- read.table(file = input, header = TRUE, sep = "\t")
colnames(tide.mat)[1] <- "sample"
tide.mat$group <- ifelse(tide.mat$TIDE >= 0, "non-responder","responder")

### read in normalized expression data and calculate CTL
expr.norm <- read.table(expr, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
CTL.entrez <- c("925","926","3001","3002","5551")  #CD8A:925;CD8B:926;GZMA:3001;GZMB:3002;PRF1;5551
expr.norm.ctl <- colMeans(expr.norm[CTL.entrez,])

####barplot of showing tide score
Pwidth <- nrow(tide.mat)*40+500
png(file=paste(outdir,"tide_score_mqc.png",sep = ""), res = 300, width = Pwidth, height = 1400)
ggplot(tide.mat,aes(x = reorder(sample, -TIDE), y = TIDE, fill = group))+
  geom_bar(stat = "identity", position="identity",alpha = 0.8)+
  theme_minimal()+
  labs(y = "TIDE Score",x = "")+
  scale_fill_manual(values = c("non-responder" = "blue", "responder" = "red"), name = "")+
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle=-90,size=6,face = "bold",hjust=1),
        axis.text.y = element_text(size = 6,face = "bold"),
        axis.title.x = element_text(size = 6,face = "bold"),
        axis.title.y = element_text(size = 6,face = "bold"))
dev.off()

####heatmap showing suppressive score compared with TCGA dataset
####prepare user dataset
tmp <- tide.mat
rownames(tmp) <- tide.mat$sample
tmp <- subset(tmp, select = c("MDSC","CAF","TAM.M2","Exclusion","Dysfunction","TIDE"))
tmp[is.na(tmp)] <- 0
# colnames(tmp)[which(colnames(tmp) == "TAM.M2")] <- "M2"
user.mat <- rbind(t(tmp), CTL = expr.norm.ctl[rownames(tmp)])
user.mat.norm <- as.matrix(t(apply(user.mat, 1, function(x) (2*(x-min(x))/(max(x)-min(x)))-1)))
user.mat.norm <- user.mat.norm[c("MDSC","CAF", "TAM.M2","Dysfunction","Exclusion","CTL","TIDE"),]

####prepare TCGA dataset
# read in CTL score 
CTL <- get(load("static/tide/Results/TCGA_CTL_score.Rdata"))
# read in TIDE results
TIDE.res <- get(load("static/tide/Results/TCGA_TIDE_res.Rdata"))
# prepare dysfunction and exclusion score
# tcga.path <- dir(path = "~/Documents/pipeline_test/RIMA_V1/rnaseq_pipeline/static/tide/Results/Tumor_Dysf_Excl_scores", pattern = "OS_base$", full.names = TRUE)
# compare.path <- grep(paste("TCGA.",cc,".*OS_base", sep = ""),tcga.path,value = TRUE)
if (cc %in% names(TIDE.res)){
  pats <- intersect(rownames(TIDE.res[[cc]]), names(CTL[[cc]]))
  tcga.tide <- cbind(TIDE.res[[cc]][pats, ], CTL = CTL[[cc]][pats])
  # tcga.tide <- read.table(compare.path, sep = "\t", row.names = 1, header = TRUE)
  tcga.mat <- t(tcga.tide)
  tcga.mat.norm <- as.matrix(t(apply(tcga.mat, 1, function(x) (2*(x-min(x))/(max(x)-min(x)))-1)))
  tcga.mat.norm <- tcga.mat.norm[c("MDSC","CAF", "TAM.M2","Dysfunction","Exclusion","CTL","TIDE"),]
  ###top tissue annotation
  ha1 = HeatmapAnnotation(User_Data = c(rep("User_Data",ncol(user.mat))), 
                          col = list(User_Data = c("User_Data" = "pink")),
			  show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(TCGA_Data = c(rep("TCGA_Data",ncol(tcga.mat))), 
                          col = list(TCGA_Data = c("TCGA_Data" = "royalblue")),
			  show_annotation_name = FALSE)
  
  ##heatmap
  mx <- max(user.mat, tcga.mat)
  mn <- min(user.mat, tcga.mat)
  col_fun <- colorRamp2(c(-1,0,1), c("#377EB8", "white", "#E41A1C"))
  if (ncol(tcga.mat.norm)/ncol(user.mat.norm)>5){
    Pwidth <- 2
  } else {
    Pwidth <- 4*ncol(user.mat.norm)/ncol(tcga.mat.norm)
  } 
  ht_list = Heatmap(user.mat.norm,  name = "Score",col = col_fun,
                    show_row_names = TRUE,
                    show_column_names = FALSE,
                    show_column_dend = FALSE,
                    show_row_dend = TRUE,
                    cluster_rows = FALSE,
                    cluster_columns = TRUE,
                    top_annotation = ha1,
                    width = unit(Pwidth,"cm"),
                    heatmap_legend_param = list(direction = "horizontal" ), #title_position = "lefttop-rot"
                    row_names_gp = gpar(fontsize = 6)) +
    Heatmap(tcga.mat.norm,  name = "Score",col = col_fun,
            show_row_names = TRUE,
            show_column_names = FALSE,
            show_column_dend = FALSE,
            show_row_dend = TRUE,
            cluster_rows = FALSE,
            cluster_columns = TRUE,
            top_annotation = ha2,
            width = unit(8,"cm"),
            #heatmap_legend_param = list(direction = "horizontal" ), #title_position = "lefttop-rot"
	    heatmap_legend_param = list( legend_direction="horizontal", legend_width=unit(3,"mm")),
            row_names_gp = gpar(fontsize = 6))
  png(file=paste(outdir,"TIDE-TCGA_mqc.png",sep = ""), res = 300, width = 2000, height = 700)
  draw(ht_list,annotation_legend_side = "bottom", heatmap_legend_side = "bottom",merge_legends = TRUE)
  dev.off()
} else{
  print("No comparable predicted TCGA tide result")
  png(file=paste(outdir,"TIDE-TCGA_mqc.png",sep = ""), res = 300, width = 2000, height = 700)
  print(ggplot() + theme_void())
  dev.off()
}
