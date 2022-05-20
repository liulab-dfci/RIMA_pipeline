library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(optparse)

hla_oncoplot <- function(hla, exprn, meta, groups) {
  
  ### separte the classI and classII
  class1 <- c("A1", "A2", "B1", "B2", "C1", "C2")
  class2 <- c("DQA11", "DQA12", "DQB11", "DQB12", "DRB11", "DRB12")
  class <- c(class1, class2)
  #detect whether there are missing value
  hla <- hla[c("subject", class)]
  for (i in class) {
    hla[[i]] <- sub("P","", hla[[i]])
  }
  
  
  #hla_table=as.data.frame(hla)
  hla[hla==""] <- "n"
  
  #remove the duplicates
  if (nrow(meta) > length(unique(meta$PatName))){
    patein_id <- meta["PatName"]
    patein_id$subject <- rownames(patein_id)
    tmp_1 <- merge(hla, patein_id, by = "subject", all = TRUE)
    new_hla <- NULL
    for (i in unique(tmp_1$PatName)) {
      tmp <- subset(tmp_1, tmp_1$PatName == i)
      tmp$quality <- rowSums(tmp == "n")
      tmp <- subset(tmp, tmp$quality == min(tmp$quality))
      if (nrow(tmp) > 1) {
        tmp <- tmp[1,]
      }
      tmp <- tmp[c("PatName", class1, class2, "subject")]
      new_hla <- rbind(new_hla, tmp)
    }
    new_meta <- new_hla[c("subject", "PatName")]
    new_hla <- new_hla[-14]
    colnames(new_hla)[1] <- "subject"
    hla <- new_hla
  } else {
    new_meta <- data.frame(row.names = rownames(meta), subject = rownames(meta),
                           PatName = meta$PatName, stringsAsFactors = F)
    new_meta <- new_meta[as.character(hla$subject),]
    hla$subject <- plyr::mapvalues(x = hla$subject,
                                   from = rownames(new_meta), to = new_meta$PatName)
  }
  
  
  
  ###count hla allele frequency for all samples
  tem = 1
  for (i in c(1:2)) {
    if (i == 1){
      i = append("subject", class1)
    } else {
      i = append("subject", class2)
    }
    hla.tem <- hla[colnames(hla)%in%i]
    hla.melt <- melt(hla.tem, id.vars = c("subject")) %>% 
      mutate_at(.vars = c("subject","variable"), .funs = as.character) 
    hla.cast <- dcast(hla.melt, value~subject, value.var = "variable", fun.aggregate = function(x) paste(x, collapse = ";"))
    rownames(hla.cast) <- hla.cast$value
    hla.cast <- hla.cast[,-1]
    hla.cast <- hla.cast[!row.names(hla.cast)%in%"n",]
    if (tem == 1){
      hla.class1 <- hla.cast
      ##gene order
      alter.perc <- apply(hla.class1, 1, function(x)length(x[x != ""])) #+ length(grep(";", unlist(x)))
      ge.order1 <- names(sort(alter.perc[alter.perc != 0],decreasing = TRUE))
      if (length(ge.order1 > 20)) {
        ge.order1 <- head(ge.order1, 20)
      }
      hla.class1 <- hla.class1[ge.order1,]
    } else {
      hla.class2 <- hla.cast
      ##gene order
      alter.perc <- apply(hla.class2, 1, function(x)length(x[x != ""])) #+ length(grep(";", unlist(x)))
      ge.order2 <- names(sort(alter.perc[alter.perc != 0],decreasing = TRUE))
      if (length(ge.order2 > 20)) {
        ge.order2 <- head(ge.order2, 20)
      }
      hla.class2 <- hla.class2[ge.order2,]
    }
    tem = tem + 1
  }
  
  
  
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    A1 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#A6CEE3", col = NA))
    },
    A2 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#1F78B4", col = NA))
    },
    B1 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#B2DF8A", col = NA))
    },
    B2 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#33A02C", col = NA))
    },
    C1 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#FB9A99", col = NA))
    },
    C2 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#E31A1C", col = NA))
    },
    DQA11 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#FDBF6F", col = NA))
    },
    DQA12 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#FF7F00", col = NA))
    },
    DQB11 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#CAB2D6", col = NA))
    },
    DQB12 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#6A3D9A", col = NA))
    },
    DRB11 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 0, gp = gpar(fill = "#808080", col = NA))
    },
    DRB12 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.5, vjust = 1, gp = gpar(fill = "#A9A9A9", col = NA))
    }
  )
  
  ##gene order again
  hla.cast <- rbind(hla.class1, hla.class2)
  alter.perc <- apply(hla.cast, 1, function(x)length(x[x != ""])) #+ length(grep(";", unlist(x)))
  ge.order <- names(sort(alter.perc[alter.perc != 0],decreasing = TRUE))
  
  ##sample order based on averaged HLA expression with acending order
  hla.expr <- exprn[c("HLA-A","HLA-B","HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"), ]
  hla.expr.aver <- colMeans(hla.expr)
  
  #grouping the first feature
  #fe <- groups[1]
  
  #new_meta <- new_meta[order(new_meta[[fe]], decreasing = FALSE),]
  
  #sa.order <- rownames(new_meta)
  id.order <- as.character(new_meta$subject)
  ex.order <- names(sort(hla.expr.aver))
  ex.order <- intersect(ex.order, id.order)
  
  id.order <- merge(ex.order,new_meta, by = 1, sort = FALSE)
  id.order <- as.character(id.order$PatName)
  
  #by default, grouping by the response column
  if ("Responder"%in%colnames(meta)) {
    meta$Responder <- sub("^$", "unkown", meta$Responder)
    test <- meta["Responder"]
    rownames(new_meta) <- new_meta$subject
    new_meta <- merge(new_meta, test, by = 0, all = FALSE)
    new_meta <- new_meta[order(new_meta$Responder, decreasing = FALSE),]
    ex.order <- as.data.frame(ex.order)
    ex.order$order <- c(1:nrow(ex.order))
    ex.order_2 <- NULL
    for (o in unique(new_meta$Responder)) {
      tmp_order <- subset(new_meta, new_meta$Responder == o)
      tmp_order <- merge(tmp_order, ex.order, by = 1, all = FALSE)
      tmp_order <- tmp_order[order(tmp_order$order, decreasing = FALSE),]
      ex.order_1 <- as.character(tmp_order$subject)
      ex.order_2 <- append(ex.order_2, ex.order_1)
    }
    ex.order <- ex.order_2
    id.order <- merge(ex.order,new_meta, by = 1, sort = FALSE)
    id.order <- as.character(id.order$PatName)
    rownames(new_meta) <- new_meta$PatName
    anno.df <- new_meta[id.order, "Responder", drop = FALSE]
  } else {
    ##top annotations with group informations
    rownames(new_meta) <- new_meta$subject
    new_meta <- new_meta[-2]
    new_meta <- merge(new_meta, meta, by = 0, all = FALSE)
    rownames(new_meta) <- new_meta$PatName
    anno.df <- new_meta[id.order, groups]
    
  }
  
  
  DefaulfColorPalette <- c(
    "#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C",
    "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#A5AA99", "#BCBD22",
    "#B279A2", "#EECA3B", "#17BECF", "#FF9DA6", "#778AAE", "#1B9E77",
    "#A6761D", "#526A83", "#B82E2E", "#80B1D3", "#68855C", "#D95F02",
    "#BEBADA", "#AF6458", "#D9AF6B", "#9C9C5E", "#625377", "#8C785D"
  )
  col.list <- list()
  for(col in groups){
    col.len <- length(unique(meta[,col]))
    cols <- DefaulfColorPalette[1:col.len]
    names(cols) <- as.character(unique(meta[,col]))
    col.list[[col]] <- cols
    DefaulfColorPalette <- DefaulfColorPalette[-c(1:col.len)]
  }
  # col <- "Treatment"
  # col.assign<- setNames(brewer.pal(8, "Set2"), unique(as.character(meta[,col])))
  ##legends
  col_fun = structure(names = c("A1","A2","B1", "B2", "C1", "C2", "DQA11", "DQA12", "DQB11", "DQB12", "DRB11", "DRB12"),
                      brewer.pal(n=12, name = "Paired"))
  col_fun[c(11,12)] <- c("#808080", "#C0C0C0")
  
  ht_list <- oncoPrint(hla.cast[ge.order, id.order], get_type = function(x) strsplit(x, ";")[[1]],
                       alter_fun = alter_fun, col = col_fun, row_order = ge.order, column_order = id.order, width = unit(nrow(hla)*0.5, "cm"),
                       row_names_gp = gpar(fontsize = 10,fontface = "bold"), pct_gp = gpar(fontsize = 10,fontface = "bold"),
                       remove_empty_columns = FALSE,#remove_empty_rows = TRUE,
                       show_column_names = FALSE, 
                       top_annotation = HeatmapAnnotation(df = anno.df,
                                                          col = col.list,
                                                          HLA_Expression = anno_points(hla.expr.aver[ex.order])), # 
                       heatmap_legend_param = list(title = "Allele", at = names(col_fun), labels = names(col_fun)))
  return(ht_list)
}
