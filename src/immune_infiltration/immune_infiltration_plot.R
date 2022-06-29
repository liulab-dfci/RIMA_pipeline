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
  pdf(paste(outdir,cell,"_corr.pdf",sep = ""),  width = 6, height = 6)
  print(corrplot.mixed(methods.corr, tl.pos = "lt", upper = "circle", lower = "number", diag = "u",tl.col = "black",lower.col = col2(200), upper.col = col2(200),tl.srt = 45,mar=c(0,3,1,1), title = cell))
  dev.off()
}

