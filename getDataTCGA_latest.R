## install this package from Bioconductor
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGAToolbox")
library(RTCGAToolbox)
library(TCGAbiolinks)
library(dplyr)
library(limma)
library(MASS)
library(WGCNA)
#library(igraph)
library(ggplot2)
library(reshape2)

#look at this tutorial on tcga data http://genomicsclass.github.io/book/pages/tcga.html


setwd("~/PhD/Bladder_Network")

#getFirehoseRunningDates()
getwd()
getFirehoseDatasets() # This will print the available datasets. We want Bladder Cancer (BLCA, see https://tcga-data.nci.nih.gov/docs/publications/tcga/?)
setwd("") # go to the working directory where you want to download the data
#blcaData <- getFirehoseData(dataset="BLCA", runDate="20160128",gistic2_Date="20160128",forceDownload=F, Clinic =TRUE, RNAseq2_Gene_Norm  =TRUE)

blcaData <- getFirehoseData(dataset="BLCA", runDate="20160128",gistic2_Date="20160128", Clinic =TRUE, RNAseq2_Gene_Norm = TRUE,
                            Mutation=TRUE)

blca_rnaseq <- getData(blcaData,type = "RNASeq2GeneNorm")
blca_clinical<-getData(blcaData,  "clinical")

library(SummarizedExperiment)
getFirehoseDatasets()

listSamples<-c(colnames(blca_rnaseq))

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BLCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples,
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)
BLCARnaseqSE <- GDCprepare(query)

BLCAMatrix <- assay(BLCARnaseqSE,"raw_count")

clinical <- GDCquery_clinic(project = "TCGA-BLCA", type = "clinical")

dataNorm <- TCGAanalyze_Normalization(tabDF = BLCAMatrix, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

#Then use the filtered genes for the construction of co-expression network

# I'll remove all the data that is not tumour i.e.
blca_rnaseq.tumour <- dataFilt[, which(as.numeric(substr(colnames(dataFilt), 14,15)) < 10)]
#blca_rnaseq.tumour <- blca_rnaseq[, which(as.numeric(substr(colnames(blca_rnaseq), 14,15)) < 10)]

# Convert barcodes from the rnaseq dataset to sample barcodes, which identify the patients
colnames(blca_rnaseq.tumour) <- substr(colnames(blca_rnaseq.tumour), 1,12)

RNAseq = blca_rnaseq.tumour[apply(blca_rnaseq.tumour,1,function(x) sum(x==0))<ncol(blca_rnaseq.tumour)*0.8,]

bladder_genes <-data.frame(rownames(blca_rnaseq.tumour))

#blca_genes <-read.table("~/PhD/Bladder_Network/all_bladder_genes.tsv", header = F)
blca_genes<-read.table("~/PhD/Bladder_Network/blca7seeds_cor.tsv", header = F)

colnames(blca_genes)<-"genes"

RNAseq = t(RNAseq)

#obtain only the bladder genes
#RNAseq2<-RNAseq[, match(blca_genes$genes, colnames(RNAseq))]
#RNAseq2<-t(RNAseq2)
#RNAseq2<-na.omit(RNAseq2)

RNAseq_voom = voom(RNAseq)$E


####################################################################
#RNAseq2<-RNAseq[, match(blca_genes$genes, colnames(RNAseq))]
#RNAseq2<-t(RNAseq2)
#RNAseq2<-na.omit(RNAseq2)

#RNAseq_voom = voom(RNAseq2)$E
#RNAseq_voom =  t(RNAseq_voom)


#Take the first 5000 most variant genes were selected as an heuristic cut-off as two genes without variation across samples will be 
#will be hgihly correlated

WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
#WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:1447],])
#"EGFR" %in% colnames(WGCNA_matrix)
#construction of coexpression similarity matrix using the bicor measure
SubGeneNames = colnames(WGCNA_matrix)

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
#write.table(sft$fitIndices, "Sft-thresholdLatest")

#Preliminary check of the correlation values in the coexpressed network
sampleTree = hclust(dist(WGCNA_matrix), method = "average");
sizeGrWindow(9, 7)
par(mfrow = c(1,2))
par(cex = 0.9);
#plot the power choice and connectivity

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")

softPower = 3
#calclute the adjacency matrix
#adj= adjacency(WGCNA_matrix, corFnc = "bicor", type = "unsigned", power = softPower)

#adj = s^softPower

#w = 1-adj

#cor_blc <-melt(adj)
#names(cor_blc)<-c("from", "to", "weight")
#write.table(cor_blc, "blca_coexpNetLatest.tsv", sep = "\t")        
#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(WGCNA_matrix,networkType = "unsigned", TOMType = "unsigned", power = softPower)
colnames(TOM)=rownames(TOM) = SubGeneNames

save(TOM, file = "TomData_15Nov.RData")


dissTOM = 1-TOM
#Module detection using heirarchical clustering
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = hclust(as.dist(dissTOM),method="average")


#geneTree = hclust(as.dist(w),method="average")
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3)
#png("Gene clustering.png", height=12, width=9)

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)


library(ape)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))




#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')





# Set the minimum module size
minModuleSize = 30

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(WGCNA_matrix[,restGenes], power = softPower)

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=hclust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA

TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

#Extracting the modules in the network

module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

module.order <- unlist(tapply(1:ncol(WGCNA_matrix),as.factor(dynamicColors),I))
m<-t(t(WGCNA_matrix[,module.order])/apply(WGCNA_matrix[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])





#matching clinical significance and the gene coexpression modules
WGCNA_matrix2 = WGCNA_matrix[match(clinical$bcr_patient_barcode, rownames(WGCNA_matrix)),]
not.available = which(is.na(rownames(WGCNA_matrix2))==TRUE)
WGCNA_matrix2 = WGCNA_matrix2[-not.available,]
str(WGCNA_matrix2)


clinical2 = clinical[-not.available,]