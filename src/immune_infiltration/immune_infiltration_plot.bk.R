###################Script to process immune infiltration modules###########
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(corrplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(optparse))
suppressMessages(library(reshape2))
suppressMessages(library(dendextend))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))

###make option list and parse command line
option_list <- list(
  make_option(c("-m","--meta"), type = "character",
              help = "path of meta file"),
  make_option(c("-i","--input_dir"), type = "character",
              help = "directory of immunedeconv output"),
  make_option(c("-o","--output_dir"), type = "character",
              help = "output path"),
  make_option(c("-c","--clinic_col"), type = "character",
              help = "clinic phenotype you want to compare")

)
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

###set parameters
inputdir<-opt$input_dir
outdir<-opt$output_dir
col <- opt$clinic_col


###read in meta file
meta <- read.csv(opt$meta, header = TRUE, row.names = 1)
###set methods in immunedeconv
methods <- c("cibersort_abs","timer","mcp_counter","quantiseq","xcell","epic")
###set cell types needing to be compared
comps <- c("B_cell", "CD4_T_cell", "CD8_T_cell", "DC","Macrophage", "Treg", "NK", "Neutrophil")

####format specific cell types from all methods
deconv.list <- list()
deconv.raw <- list()
for (m in methods){
  tmp <- read.table(paste(inputdir,m,".txt", sep = ""),header = TRUE, sep = "\t", row.names = 1, check.names = TRUE)
  if (m == "cibersort_abs"){
    tmp.ciber <- tmp %>%
      mutate(B_cell = B.cells.naive + B.cells.memory + Plasma.cells) %>%
      mutate(CD4_T_cell = T.cells.CD4.naive + T.cells.CD4.memory.resting + T.cells.CD4.memory.activated + T.cells.follicular.helper) %>%
      mutate(CD8_T_cell = T.cells.CD8) %>%
      mutate(DC = Dendritic.cells.resting + Dendritic.cells.activated) %>%
      mutate(Macrophage = Macrophages.M0 + Macrophages.M1 + Macrophages.M2) %>%
      mutate(Treg = T.cells.regulatory..Tregs.) %>%
      mutate(NK = NK.cells.resting + NK.cells.activated) %>%
      mutate(Neutrophil = Neutrophils)
    tmp.ciber <- subset(tmp.ciber, select = intersect(comps, colnames(tmp.ciber)))
    rownames(tmp.ciber) <- rownames(tmp)
    deconv.list[[m]] <- tmp.ciber
    deconv.raw[[m]] <- tmp[,1:22]
  }
  else if(m == "timer"){
    rownames(tmp) <- gsub("\\+","",gsub(" ","_", rownames(tmp)))
    tmp.timer <- data.frame(t(tmp)) %>%
      mutate(B_cell = B_cell) %>%
      mutate(CD4_T_cell = T_cell_CD4) %>%
      mutate(CD8_T_cell = T_cell_CD8) %>%
      mutate(DC = Myeloid_dendritic_cell) %>%
      mutate(Macrophage = Macrophage) %>%
      mutate(Neutrophil = Neutrophil)
    tmp.timer <- subset(tmp.timer, select = intersect(comps, colnames(tmp.timer)))
    rownames(tmp.timer) <- colnames(tmp)
    deconv.list[[m]] <- tmp.timer
    deconv.raw[[m]] <- data.frame(t(tmp))
  }
  else if(m == "mcp_counter"){
    rownames(tmp) <- gsub("\\+|\\/","",gsub(" ","_", rownames(tmp)))
    tmp.mcp <- data.frame(t(tmp)) %>%
      mutate(B_cell = B_cell) %>%
      mutate(CD8_T_cell = T_cell_CD8) %>%
      mutate(DC = Myeloid_dendritic_cell) %>%
      mutate(Macrophage = MacrophageMonocyte) %>%
      mutate(NK = NK_cell) %>%
      mutate(Neutrophil = Neutrophil)
    tmp.mcp <- subset(tmp.mcp, select = intersect(comps, colnames(tmp.mcp)))
    rownames(tmp.mcp) <- colnames(tmp)
    deconv.list[[m]] <- tmp.mcp
    deconv.raw[[m]] <- data.frame(t(tmp))
  }
  else if(m == "quantiseq"){
    rownames(tmp) <- gsub("\\(|\\)","",gsub("\\+|\\-","",gsub(" ","_", rownames(tmp))))
    tmp.quantiseq <- data.frame(t(tmp)) %>%
      mutate(B_cell = B_cell) %>%
      mutate(CD4_T_cell = T_cell_CD4_nonregulatory) %>%
      mutate(CD8_T_cell = T_cell_CD8) %>%
      mutate(DC = Myeloid_dendritic_cell) %>%
      mutate(Macrophage = Macrophage_M1 + Macrophage_M2) %>%
      mutate(Treg = T_cell_regulatory_Tregs) %>%
      mutate(Neutrophil = Neutrophil) %>%
      mutate(NK = NK_cell)
    tmp.quantiseq <- subset(tmp.quantiseq, select = intersect(comps, colnames(tmp.quantiseq)))
    rownames(tmp.quantiseq) <- colnames(tmp)
    deconv.list[[m]] <- tmp.quantiseq
    deconv.raw[[m]] <- data.frame(t(tmp))

  }
  else if(m == "epic"){
    rownames(tmp) <- gsub("\\+","",gsub(" ","_", rownames(tmp)))
    tmp.epic <- data.frame(t(tmp)) %>%
      mutate(B_cell = B_cell) %>%
      mutate(CD4_T_cell = T_cell_CD4) %>%
      mutate(CD8_T_cell = T_cell_CD8) %>%
      mutate(Macrophage = Macrophage) %>%
      mutate(NK = NK_cell)
    tmp.epic <- subset(tmp.epic, select = intersect(comps, colnames(tmp.epic)))
    rownames(tmp.epic) <- colnames(tmp)
    deconv.list[[m]] <- tmp.epic
    deconv.raw[[m]] <- data.frame(t(tmp))
  }
  else {
    rownames(tmp) <- gsub("\\(|\\)","",gsub("\\+|\\-","",gsub(" ","_", rownames(tmp))))
    tmp.xcell <- data.frame(t(tmp)) %>%
      mutate(B_cell = B_cell + Classswitched_memory_B_cell + B_cell_memory + B_cell_naive + B_cell_plasma) %>%
      mutate(CD4_T_cell = T_cell_CD4_memory + T_cell_CD4_naive + T_cell_CD4_nonregulatory + T_cell_CD4_central_memory + T_cell_CD4_effector_memory + T_cell_CD4_Th1 + T_cell_CD4_Th2) %>%
      mutate(CD8_T_cell = T_cell_CD8_naive + T_cell_CD8 + T_cell_CD8_central_memory + T_cell_CD8_effector_memory) %>%
      mutate(DC = Myeloid_dendritic_cell + Plasmacytoid_dendritic_cell + Myeloid_dendritic_cell_activated) %>%
      mutate(Macrophage = Macrophage_M1 + Macrophage_M2) %>%
      mutate(Treg = T_cell_regulatory_Tregs) %>%
      mutate(Neutrophil = Neutrophil) %>%
      mutate(NK = NK_cell)
    tmp.xcell <- subset(tmp.xcell, select = intersect(comps, colnames(tmp.xcell)))
    rownames(tmp.xcell) <- colnames(tmp)
    deconv.list[[m]] <- tmp.xcell
    deconv.raw[[m]] <- data.frame(t(tmp)[,1:36])
  }
}

###########------------plot1: compare method performance among 6 methods-----------------###########
###correlation between any two methods for each cell type
for (cell in comps){
  print (cell)
  col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0","#FFFFFF","#FDDBC7", "#F4A582", "#D6604D",  "#B2182B", "#67001F"))
  cell.deconv <- plyr::compact(lapply(deconv.list, function(x) if(cell %in% colnames(x)) subset(x, select = cell)))
  cell.deconv.new <- do.call("cbind", cell.deconv)
  colnames(cell.deconv.new) <- names(cell.deconv)
  methods.corr <- cor(cell.deconv.new, method = "spearman")
  png(paste(outdir,cell,"_corr.png",sep = ""), res = 350, width = 2000, height = 2200)
  # corrplot(methods.corr, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, col = col2(200),mar=c(0,0,2,0), title = cell, tl.cex = 1.2)
  print(corrplot.mixed(methods.corr, tl.pos = "lt", upper = "circle", lower = "number", diag = "u",tl.col = "black",lower.col = col2(200), upper.col = col2(200),tl.srt = 45,mar=c(0,3,1,1), title = cell))
  dev.off()
}


###########------------plot2: show immune cell abundance among 6 methods-----------------###########
###merge all data and max-min normalization
merge.df <- t(do.call("cbind", deconv.raw))
print(colnames(merge.df))
merge.df <- na.omit(merge.df)
merge.df.norm <- as.matrix(t(apply(merge.df, 1, function(x) (x-mean(x))/sd(x))))
print(merge.df.norm)
phenos <-sort(unique(as.character(meta[,col])),decreasing = TRUE)
print(phenos)
###clustering for each group
gr <- lapply(phenos, function(p){
  tmp.df <- merge.df.norm[,rownames(meta)[which(meta[,col] == p)]]
  tmp.cl <- hclust(dist(t(tmp.df),method = "euclidean"))
  hcl <- order.hclust(tmp.cl)
  ss <- rownames(meta)[which(meta[,col] == p)][hcl]
  return (ss)
})
###re-order samples based on clustering result
merge.df <- t(do.call("cbind", deconv.raw))[,unlist(gr)]
merge.df.norm <- as.matrix(t(apply(merge.df, 1, function(x) (x-mean(x))/sd(x))))# #as.matrix(t(apply(merge.df, 1, function(x) ((x-min(x))/(max(x)-min(x)))))) #

###common heatmap params for showing results from different methods
method_col <- brewer.pal(n = 6, name = "Set1")
names(method_col) <- names(deconv.raw)
method_split <- unlist(lapply(rownames(merge.df), function(x) strsplit(x, "\\.")[[1]][1]))
###annotation legend parameters
ht_opt(
  legend_title_gp = gpar(fontsize = 14, fontface = "bold"),
  legend_labels_gp = gpar(fontsize = 14),
  heatmap_column_names_gp = gpar(fontsize = 14),
  heatmap_column_title_gp = gpar(fontsize = 14),
  heatmap_row_title_gp = gpar(fontsize = 20),
  heatmap_row_names_gp = gpar(fontsize = 14)
)
Pheight <- length(method_split)*110
png(paste(outdir,"ImmuneDeconv_heatmap.png",sep = ""), res = 450, width = 7000, height =Pheight)
###compute significance level only when phenotypes >= 2
if(length(phenos) > 1){
  if(length(phenos) == 2){
    merge.df.melt <- melt(merge.df) %>% mutate(Group = as.character(meta[Var2,col]))
    df.tmp <- merge.df.melt %>% group_by(Var1,Group) %>% dplyr::summarise(Mean = mean(value))
    df.tmp.comp <- unlist(lapply(as.character(unique(df.tmp$Var1)), function(c){
      gap <- subset(df.tmp, Var1 == c & Group == phenos[1], select = "Mean") - subset(df.tmp, Var1 == c & Group == phenos[2], select = "Mean")
      gap <- gap[1,1]
      if(gap >0) {return(paste("Higher_in_",phenos[1], sep = ""))} else if(gap < 0){return (paste("Higher_in_",phenos[2], sep = ""))} else{return ("equal")}
    }))
    names(df.tmp.comp) <- as.character(unique(df.tmp$Var1))
    comp.sig <- data.frame(compare_means(value~Group, data = merge.df.melt, group.by = "Var1", paired = FALSE)) %>%
      mutate(Means = df.tmp.comp[as.character(Var1)]) %>%
      mutate(Means_sig = ifelse(p.signif != "ns", paste(Means, p.signif, sep = ""), p.signif))
    rownames(comp.sig) <- as.character(comp.sig$Var1)
    ###add annotation for significance level
    sig_col <- c("#FED976","#FD8D3C","#E31A1C","#800026","lightgray","#C7E9B4","#41B6C4","#225EA8","#081D58")
    names(sig_col) <- c(paste(paste("Higher_in_",phenos[1], sep = ""), c("*","**","***","****"), sep = ""), "ns",
                        paste(paste("Higher_in_",phenos[2], sep = ""), c("*","**","***","****"), sep = ""))
    cell_sig <- as.character(comp.sig[rownames(merge.df),"Means_sig"])
  } else {
    merge.df.melt <- melt(merge.df) %>% mutate(Group = as.character(meta[Var2,col]))
    df.tmp <- merge.df.melt %>% group_by(Var1,Group) %>% dplyr::summarise(Mean = mean(value))
    comp.sig <- data.frame(compare_means(value~Group, data = merge.df.melt, group.by = "Var1", paired = FALSE, method = "kruskal.test"))
    rownames(comp.sig) <- as.character(comp.sig$Var1)
    sig_col <- c("#FED976","#FD8D3C","#E31A1C","#800026","lightgray")
    names(sig_col) <- c("*","**","***","****","ns")
    cell_sig <- as.character(comp.sig[rownames(merge.df),"p.signif"])
  }
  ###top anntation for samples grouping
  Tissue_col <- structure(names = phenos, brewer.pal(n = 6, name = "Set1")[1:length(phenos)])
  ha.g3 <- columnAnnotation(
    Tissue = as.character(meta[colnames(merge.df),col]),
    col = list(Tissue = Tissue_col),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 14)
  )
  ###heatmap when phenos >= 2
  ht_list = Heatmap(merge.df.norm, name = "Estimation",
                    # col = col_fun,
                    # clustering_distance_columns = "euclidean",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    # show_column_dend = TRUE,
                    # show_row_dend = FALSE,# rect_gp = gpar(col= "white"),
                    show_column_names = FALSE,#row_split = g.group,
                    row_split = method_split,
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 12),
                    row_names_max_width = unit(10,"cm"),
                    # left_annotation = ha.g4,
                    # right_annotation = ha.g2,
                    # bottom_annotation = ha.g0,
                    top_annotation = ha.g3)+
    Heatmap(method_split, name = "Method", col = method_col, column_names_side = "bottom")

  ht_list <- ht_list + Heatmap(cell_sig, name = "Signif", col = sig_col, column_names_side = "bottom")

} else{
  ht_list = Heatmap(merge.df.norm, name = "Estimation",
                    # col = col_fun,
                    # clustering_distance_columns = "euclidean",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    # show_column_dend = TRUE,
                    # show_row_dend = FALSE,# rect_gp = gpar(col= "white"),
                    show_column_names = FALSE,#row_split = g.group,
                    row_split = method_split,
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 12),
                    row_names_max_width = unit(10,"cm")
                    )+
    Heatmap(method_split, name = "Method", col = method_col, column_names_side = "bottom")
}
####output figrue
#draw(ht_list, padding = unit(c(20, 2, 2, 2), "mm"),
     merge_legends = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right")
#ht_opt(RESET = TRUE)
#dev.off()


#preprocess the data for multiqc
#
#order_ta <- function(x, col) {
  #col_order
#  meta <- meta[order(meta[[col]], decreasing = TRUE),]
#
#  orde_col <- data.frame("samples" = rownames(meta), num = 1:nrow(meta), pheo = meta[[col]])
#  re_col <- data.frame("samples" = colnames(x), re_num = 1:nrow(meta))
#
#  orde_col <- merge(orde_col, re_col, by = 1)
#  orde_col <- orde_col[order(orde_col[["num"]], decreasing = FALSE),]
#
  #cluster by the column
#  t <- rownames(x)
#  for (sa in unique(orde_col$pheo)) {
#    tmp_sa <- subset(orde_col, pheo == sa)
#    tmp_1 <- x[tmp_sa[["re_num"]]]
#    nw_order <- data.frame("ncol" = 1:ncol(tmp_1), "colsum" = colSums(tmp_1))
#    nw_order <- nw_order[order(nw_order$colsum, decreasing = TRUE),]
#    tmp_1 <- tmp_1[nw_order[["ncol"]]]
#    t <- cbind(t, tmp_1)
#  }
#  t <- t[-1]

  #row_order based on expression
 # t$ave <- rowSums(t)
 # t <- t[order(t$ave, decreasing = TRUE),]
 # t <- t[-ncol(t)]
 # return(t)
#}


#check the significant level
#t_1 <- comp.sig[c("Var1","p.signif")]
#t_2 <- merge.df
#t_2 <- cbind(rownames(t_2), t_2)
#colnames(t_2)[1] <- "Var1"
#test <- merge(t_1, t_2, by = "Var1", all = FALSE, stringsAsFactors = FALSE)

#mutiqc <- opt$multiqc

#preprocess the table for each tools
#for (i in methods) {
#  tmp <- test[str_detect(test$Var1, i),]
#  tmp$Var1 <- gsub(paste(i,".",sep = ""),"", tmp$Var1)
#  tmp$Var1 <- paste0(tmp$p.signif, tmp$Var1, sep = "")
#  tmp$Var1 <- gsub("^ns","", tmp$Var1)
#  tmp <- tmp[-2]

#  write.table(tmp, file = paste(mutiqc, paste(i,"immune_multiqc.txt", sep = "_"), sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
#}

#files <- list.files(mutiqc)

#for (f in files) {
#  tmp <- read.delim(paste(mutiqc, f, sep = ""), row.names = 1)
#  tmp <- order_ta(tmp, col)
#  write.table(tmp, file = paste(mutiqc, f, sep = ""), quote = FALSE, col.names = NA, sep = "\t")
#}
