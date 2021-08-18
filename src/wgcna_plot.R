library(WGCNA)
library(reshape2)
library(stringr)
library(optparse)

# make option list and parse command line
option_list <- list(  
  make_option(c("-e", "--expression_dat"), type="character", 
              help="Input path of expression file. [Required]"),
  make_option(c("-t", "--traits_data"), type="character",
              help="Input path of corresponding batch file[Required]"),
  make_option(c("-c", "--clinical_phenotype"), type="character",
              help="names of phenotypes want to associate with gene modules"),
  make_option(c("-n", "--top_n"), type="numeric",  
              help="top n genes with higher mad"),
  make_option(c("-m", "--minModuleSize"), type="numeric",  
              help="covariates needs to be adjusted for"),
  make_option(c("-a", "--nThreads"), type="numeric",  
              help="number of threads to allow"),
  make_option(c("-o","--output",type="character", help="Output files [Required]"))
)
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# paramenter checking
if(is.null(opt$expression_dat)) stop('Expression file required.')  ###if not provide batch file output log expression matrix

# n <- 10000
# minModuleSize <- 30
# trait.col <- 2
# outdir <- "~/Documents/rnaseq/data/"
# 
# ##########--------------read in expression data-------------########
# exprMat <- "~/Documents/rnaseq/data/tpm_convertID.txt"
# dataExpr <- read.table(exprMat,sep=',', row.names=1, header=TRUE, quote="", comment="", check.names=F)
# dim(dataExpr)
# ##########--------------read in clinical information--------------##########
# datTraits <- read.table("~/Documents/rnaseq/data/metasheet.csv",header = TRUE, sep = ",",row.names = 1)
# rownames(datTraits) <- rownames(dataExpr)
# head(datTraits)


#########read in expression data
dataExpr <- read.table(opt$expression_dat,sep=',', row.names=1, header=TRUE, quote="", comment="", check.names=FALSE)
#########read in clinical information
datTraits <- read.table(opt$traits_data,header = TRUE, sep = ",",row.names = 1)
####overlap samples
tmp <- intersect(colnames(dataExpr),rownames(datTraits))
dataExpr <- dataExpr[,tmp]
#dataExpr <- dataExpr[,rownames(datTraits)]

#####other params
n <- opt$top_n
minModuleSize <- opt$minModuleSize
all.traits <- opt$clinical_phenotype
outdir <- opt$output
nThreads <- opt$nThreads

###for weighted gene network, point represents genes, edge represents the correlation between gene x and gene y, 
##"unsigned": undirected weighted network(don't care about positie or negative correlation); 
##"signed": directed weighted network(care about positive or negative correlation)
type = "signed"  
corType = "pearson"
# open multiple process
enableWGCNAThreads(nThreads = nThreads)

#=====================================================================================
#
#  Code chunk 1: Removing outlier genes and samples
#
#=====================================================================================
# filterign genes with more than 75% mad
m.mad <- apply(dataExpr,1,mad)
dataExpr <- dataExpr[order(m.mad, decreasing = TRUE)[1:n],]
m.mad <- m.mad[rownames(dataExpr)]
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
###detecting missing value 
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
# ###check outlier samples
# sampleTree = hclust(dist(dataExpr), method = "average")
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")




#=====================================================================================
#
#  Code chunk 2: Choose a set of soft-thresholding powers
#
#=====================================================================================
png(paste(outdir,"soft_thresholding_power_plot.png",sep = ""), res = 300, width = 3000, height = 2000)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5, networkType = type)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
####get power(empirical power value will be chosed if no proper power)
power = sft$powerEstimate
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
power
dev.off()

#=====================================================================================
#
#  Code chunk 3: One-step network construction and module detection 
#
#=====================================================================================
###set one module contains at least 30 genes
png(paste(outdir,"block_module_plot.png",sep = ""), res = 300, width = 2000, height = 1200)
minModuleSize <- minModuleSize
net = blockwiseModules(dataExpr, power = power,
                       TOMType = type, corType = corType,
                       minModuleSize = minModuleSize, maxBlockSize = nGenes,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "WGCNA_TOM",verbose=5)
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# #=====================================================================================
# #
# #  Code chunk 4: construct gene co-expression network 
# #
# #=====================================================================================
png(paste(outdir,"co-expression_network_plot.png",sep = ""), res = 350, width = 2500, height = 2500)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
## loading TOM data from above net
load(net$TOMFiles[1], verbose=TRUE)
# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
# TOM plot???clustering for rows and columns
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA  # Set diagonal to NA for a nicer plot
TOMplot(TOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

# #=====================================================================================
# #
# #  Code chunk 5: relate modules to external information
# #
# #=====================================================================================
RelateModTriat <- function(moduleTraitCor, MEs, design, textMatrix){
  # Display the correlation values within a heatmap plot
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}



# #=====================================================================================
# #
# #  Code chunk 6: figure out the relationship between interested modules, traits and genes
# #
# #=====================================================================================
RelateModTriatGene <- function(dataExpr, MEs, TOM, Mod, treatment, design, nSamples){
  #calculate correlation matrix for module and genes
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  # calculate correlation matrix for traits and genes(continous varible,or binary variable)
  treat = as.data.frame(design[,treatment])
  names(treat) = treatment
  geneTraitSignificance = as.data.frame(cor(dataExpr, treat, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(treat), sep="");
  names(GSPvalue) = paste("p.GS.", names(treat), sep="");
  
  #combine two correlation matrix and analyze for interested module
  module = Mod
  column = match(module, modNames)
  moduleGenes = moduleColors==module;
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for ",treatment,sep = ""),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,col = module)
  
  # extract genes in interested module
  probes = colnames(dataExpr) 
  inModule = (moduleColors==module);
  modProbes = probes[inModule]
  # export module to cytoscape
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  )
}



# #=====================================================================================
# #
# #  Code chunk 7: find the key drivers in interesting module
# #
# #=====================================================================================
### Intramodular connectivity, module membership, and screening for intramodular hub genes

ModuleHubGenes <- function(dataExpr, Mod, treatment, design, moduleColors){
  # Intramodular connectivity
  connet=abs(cor(dataExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  
  #  Relationship between gene significance and intramodular connectivity
  which.module = Mod
  treat = as.data.frame(design[,treatment]); # change specific 
  names(treat) = treatment
  GS1 = as.numeric(cor(treat,dataExpr, use="p"))
  GeneSignificance=abs(GS1)
  colorlevels=unique(moduleColors)
  whichmodule = Mod
  restrict1 = (moduleColors==whichmodule);
  
  #     par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  #     par(mar = c(4,5,3,1))
  #     for (i in c(1:length(colorlevels)))
  #     {
  #       whichmodule=colorlevels[[i]];
  #       restrict1 = (moduleColors==whichmodule);
  #       verboseScatterplot(Alldegrees1$kWithin[restrict1],
  #                          GeneSignificance[restrict1], col=moduleColors[restrict1],
  #                          main=whichmodule,
  #                          xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  #     }
  
  # Generalizing intramodular connectivity for all genes on the array
  datKME=signedKME(dataExpr, MEs, outputColumnName="")
  ##Finding genes with high gene significance and high intramodular connectivity in interesting modules
  # abs(GS1)> .9 ????????????????????????????????????
  # abs(datKME$MM.black)>.8 ???????????? >0.8
  FilterGenes= abs(GS1)> .9 & abs(datKME[,Mod])>.8
  table(FilterGenes)
  hubgenes = rownames(datKME)[FilterGenes]
  datKME$Genes <- rownames(datKME)
  hub <- datKME[FilterGenes,c("Genes",Mod)]
  colnames(hub) <- c("Hub_genes","Connectivity")
  write.table(as.data.frame(hub),file = paste0(outdir,"/",trait.col,"_",names(Mod),"_",Mod,"_hubgenes.txt"),quote = FALSE,col.names = TRUE,row.names =FALSE,sep = '\t')
  pr = paste("There is ",length(hubgenes)," hub genes in this module: ", sep = "")
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=paste(pr, whichmodule, sep = "\n"),
                     xlab = "Connectivity",
                     ylab = "Gene Significance", abline = TRUE)
}


###invoke function in chunk5, chunk 6 and chunk 7 for different traits
if(!is.null(datTraits) && !is.null(all.traits)){
  for(trait.col in unlist(strsplit(all.traits, ","))){
    phenos <- unique(datTraits[,trait.col])
    if(length(phenos) > 1){
      ##generate design matrix for a clinical trait
      design <- data.frame(model.matrix(~0+ datTraits[,trait.col]))
      colnames(design) <- phenos
      ####module-traits relationships
      moduleColors <- labels2colors(net$colors)
      # Recalculate MEs with color labels
      MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
      MEs = orderMEs(MEs0) ##ME matrix for different modules
      moduleTraitCor = cor(MEs, design , use = "p")
      moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
      #display correlations and their p-values
      textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                         signif(moduleTraitPvalue, 1), ")", sep = "");
      dim(textMatrix) = dim(moduleTraitCor)
      textMatrix
      # get the module positively correlated a specific trait most
      corr.max <- apply(moduleTraitCor,2,max)
      ModHighCorr <- lapply(names(corr.max), function(x){
        mod.pos <- which(moduleTraitCor[,x] == corr.max[x])
        if( corr.max[x] > 0 &&  moduleTraitPvalue[mod.pos,x] <= 0.05){
          return(gsub("ME","",rownames(moduleTraitCor)[mod.pos]))
        }
        else{
          return("")
        }
      })
      interestedMod <- do.call("c",ModHighCorr)
      names(interestedMod) <- names(corr.max)
      interestedMod ###get a module for each trait
      tryCatch({
        ###Code chunk 5: relate modules to external information
        Pwidth <- ncol(textMatrix)*1000
        Pheight <- nrow(textMatrix)*350
        png(paste(outdir,trait.col,"_module_traits_plot.png",sep = ""), res = 350, width = Pwidth, height = Pheight)
        RelateModTriat(moduleTraitCor, MEs, design = design, textMatrix)
        print("plot of associating module and traits")
        dev.off()
      }, error = function(e) {cat("ERROR:",conditionMessage(e),"\n")})
      for(trs in names(interestedMod)){
        mod <- interestedMod[trs]
        tryCatch({
          ###Code chunk 6: figure out the relationship between interested modules, traits and genes
          png(paste(outdir, trait.col, "_", trs, "_Module_Genes_plot.png", sep = ""), res = 300, width = 2000, height = 2000)
          RelateModTriatGene(dataExpr, MEs, TOM, Mod = mod, treatment = trs, design = design, nSamples)
          print("plot of associating traits and module genes")
          dev.off()
          ###Code chunk 7: find the key drivers in interesting module
          png(paste(outdir, trait.col, "_", trs, "_Connectivity_Genes_plot.png", sep = ""), res = 300, width = 2200, height = 2200)
          ModuleHubGenes(dataExpr, Mod = mod, treatment = trs, design = design, moduleColors)
          print("plot of hub genes")
          dev.off()
        }, error = function(e) {cat("ERROR:",conditionMessage(e),"\n")})
      }
    }
  }
}

