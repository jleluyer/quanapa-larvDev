#!/usr/bin/Rscript

ls()
rm(list=ls())
ls()

# Installing packages bioconductor DEseq2
# Loading packages
library('DESeq2')
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("IHW")

#Global variables
setwd("")


# Load data ----
# Loading data
coldata<-read.table("data.trait.txt", sep="\t",header=T,na.string="NA",check.names = FALSE)
rownames(coldata) <-coldata$Ind
coldata <- coldata[order(rownames(coldata)),]
coldata$Ind <- NULL
rownames(coldata)


cts<-read.table("expression.matrix.txt", header=T,check.names=FALSE)
rownames(cts)<-cts$genes
cts$genes <- NULL
colnames(cts)

#cts <- cts[,order(vectorname)] 
cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))]
rownames(coldata)
colnames(cts)
coldata$Condition


#check 
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))


coldata$Condition<-as.factor(coldata$Condition)


# Prepare file ----
str(cts)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)
#name_filtered <-rownames(dds)

colSums(cts)

# filtering
dds <- estimateSizeFactors(dds)
dds$sizeFactor


#For test  Normalisation des donn?es
idx <- rowSums(counts(dds,normalized=TRUE) >= 10 ) >= 3  # chec use normalize
dds2<-dds
dds <- dds[idx,]
dim(dds)
dim(dds2)
#matrix log
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)

#matrix log
norm.counts.com <- counts(dds2, normalized=TRUE)
log.norm.counts.com <- log2(norm.counts.com + 1)
#write.table(log.norm.counts,file="log2CPM_larvaire.txt",quote=F)
log.norm.counts<-as.data.frame(log.norm.counts.com)


# Data transformation pour que l'ACP se fasse sur un gros jeu de donn?es

vsd.fast <- vst(dds, fitType='local',blind=FALSE)
vsd.assay<-assay(vst(dds, blind=TRUE))
write.table(vsd.assay, file ="vsd_assay_devlarve.txt", sep= "\t")

vsd.fast.comp <- vst(dds2, fitType='local',blind=FALSE)
vsd.assay.com<-assay(vst(dds2, blind=TRUE))
write.table(vsd.assay.com, file ="vsd_assay_devlarve_comp.txt", sep= "\t")

# plot PCA ----
plotData<-plotPCA(vsd.fast,intgroup=c("Condition"), returnData=TRUE)
plotData$Condition<-factor(plotData$Condition,levels=c("D","velyger","umbo","oeillee"))
percentVar <- round(100 * attr(plotData, "percentVar"))
PCA<-ggplot(plotData, aes(PC1, PC2, fill=Condition)) + geom_point(size=6,shape=21,alpha=0.8)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed()
PCA

theme_glob<-theme(axis.text.x=element_text(colour="black",size=14),
                  axis.text.y=element_text(colour="black",size=14),
                  axis.title.x=element_text(colour="black",size=14),
                  axis.title.y=element_text(colour="black",size=14),
                  panel.background = element_blank(),
                  panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"),
                  legend.title=element_blank(),
                  aspect.ratio=1,
                  legend.background = element_rect(size=0.3, linetype="solid", colour ="black"),
                  legend.position = c(0.85,0.2),
                  legend.key = element_rect(fill=NA),
                  legend.text = element_text(colour="black",size=12)) 
couleurs=c("turquoise","dodgerblue2","orange2","red") 
pdf(file = "PCA_larv.pdf", width = 6, height = 6)
set.seed(1)
graph.pca<-PCA +theme_glob + scale_fill_manual(values=couleurs,labels=c("D-shaped","Velyger","Umbo","Eye-spot"))+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")
graph.pca
dev.off()



# Differential expression contrasts----
dds <- DESeq(dds, test="LRT", reduced=~1) # Likelihood ratio test
res <- results(dds,alpha=0.01)
res <- res[order(res$padj),]
resultsNames(dds)

# D vs velyger
resDvel<- results(dds, contrast=c("Condition","D","velyger")) # FDR of 0.01
res2andVel<-resDvel
summary(resDvel)
resDvel <- resDvel[order(resDvel$pvalue),]
resDvel<-as.data.frame(resDvel)
resDvel$genes<-rownames(resDvel)
rownames(resDvel) <-NULL
resDvel <- resDvel[,c(7,2)]
D2Vel<-resDvel
res2andVel <- res2andVel[order(res2andVel$pvalue),]
res2andVel<-subset(res2andVel,abs(log2FoldChange)>2)
res2andVel<-subset(res2andVel,padj<0.01)
res2andVel<-as.data.frame(res2andVel)
write.table(rownames(res2andVel),file='DEG_D2Vel.txt',quote=F,sep=",",row.names=F,col.names = F)
write.table(res2andVel,file='D2Vel.txt',quote=F,sep=",")

# velyger vs umbp
resVelUmb<- results(dds, contrast=c("Condition","velyger","umbo")) # FDR of 0.01
summary(resVelUmb)
resVelUmb <- resVelUmb[order(resVelUmb$pvalue),]
resVelUmb<-as.data.frame(resVelUmb)
resVelUmb$genes<-rownames(resVelUmb)
rownames(resVelUmb) <-NULL
resVelU <- resVelUmb[,c(7,2)]
Vel2Umb<-resVelU

V2Umb<-resVelUmb
V2Umb <- V2Umb[order(V2Umb$pvalue),]
V2Umb<-subset(V2Umb,abs(log2FoldChange)>2)
V2Umb<-subset(V2Umb,padj<0.01)
V2Umb<-as.data.frame(V2Umb)
write.table(rownames(V2Umb),file='DEG_V2Umb.txt',quote=F,sep=",",row.names=F,col.names = F)
#write.table(resultsvelUmb,file='resultVelUmbl.csv',quote=F,sep=",")

# unmbo vs oeille
resUmbOeil<- results(dds, contrast=c("Condition","umbo","oeillee")) # FDR of 0.01
summary(resUmbOeil)
resUmbOeil <- resUmbOeil[order(resUmbOeil$pvalue),]
resUmbOeil<-as.data.frame(resUmbOeil)
resUmbOeil$genes<-rownames(resUmbOeil)
rownames(resUmbOeil) <-NULL
resUmbOeil <- resUmbOeil[,c(7,2)]
Umb2Oeil<-resUmbOeil

UvsO<-resUmbOeil
UvsO <- UvsO[order(UvsO$pvalue),]
UvsO<-subset(UvsO,abs(log2FoldChange)>2)
UvsO<-subset(UvsO,padj<0.01)
UvsO<-as.data.frame(UvsO)
write.table(rownames(UvsO),file='DEG_U2O.txt',quote=F,sep=",",row.names=F,col.names = F)
write.table(resultsUmbOeil,file='resultUmbOeil.csv',quote=F,sep=",")


# plot single genes for comparison

#gene test: evm.TU.scaffold118size383775.19 LogFC2 -8 D.2.Velyger
#gene test: evm.TU.scaffold495size234749.12 LogFC2 6 D.2.Velyger
# Si > 0, Velyger > D
#Si < 0, D > Velyger

library(ggplot2)
merge.norm<-merge(t(vsd.assay),coldata,by=0)

ggplot(merge.norm,aes(x=Condition,y=evm.TU.scaffold495size234749.1))+geom_boxplot()
       
# KOGmwu ----

library(KOGMWU)
library(pheatmap)
library(RColorBrewer)


gene2kog<-read.table("../genome_annotation/keggs/pmarg_gene2kogClass.TU.tab",header=F,sep='\t') # kog class annotations for a.millepora, included with KOGMWU package
head(gene2kog)

# launch:
d2vel.lth=kog.mwu(D2Vel,gene2kog)
vel2umb.lth=kog.mwu(Vel2Umb,gene2kog)
umb2oeil.lth=kog.mwu(Umb2Oeil,gene2kog)


#Make table
kogtable=makeDeltaRanksTable(list("D.2.Velyger"=d2vel.lth,"Velyger.2.Umbo"=vel2umb.lth,
                                      "Umbo.2.Eyed"=umb2oeil.lth))



# pheatmap
library(pheatmap)
mat<-as.matrix(kogtable)
mat <- mat[complete.cases(mat), ]
summary(mat)

mat<-mat[rownames(mat) != "",]

pdf(file="pheatmap_kog.pdf")
seed(10)
pheatmap(mat,clustering_distance_cols="correlation",
         treeheight_row=15,treeheight_col=15,
         border_color="white",
         cellheight=14, cellwidth = 14)
dev.off()
 
# exploring correlations between datasets
pairs(kogtable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(kogtable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

save(kogtable,file="kogtable.Rda")

library(ggplot2)
library(ggh4x)
load("kogtable.Rda")
# individual KOG class correlations

theme_cor<-theme(axis.text.y   = element_text(size=20),
                 axis.text.x   = element_text(size=20),
                 axis.title.y  = element_text(size=20),
                 axis.title.x  = element_text(size=20),
                 plot.title = element_text(hjust = 0.5,size=20),
                 panel.background = element_blank(),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 panel.border = element_blank(),
                 aspect.ratio=1
)

cor.test(kogtable$D.2.Velyger,kogtable$Velyger.2.Umbo)
plot.DVelVelUmb <-ggplot(kogtable,aes(x=D.2.Velyger,y=Velyger.2.Umbo)) +geom_point(size=3,shape=21) +theme_cor +
  guides(x = "axis_truncated", y = "axis_truncated") +
  geom_smooth(method = lm, se = FALSE,color="red")+
  ggtitle("r = 0.34, p = 0.10")


cor.test(kogtable$Umbo.2.Eyed,kogtable$Velyger.2.Umbo)
plot.VelUmb_UmbEyed <-ggplot(kogtable,aes(y=Umbo.2.Eyed,x=Velyger.2.Umbo)) +geom_point(size=3,shape=21) +theme_cor +
  guides(x = "axis_truncated", y = "axis_truncated") +
  geom_smooth(method = lm, se = FALSE,color="red")+
  ggtitle("r = -0.38, p = 0.07")


# plot mutliple graphs
library(ggpubr)
pdf(file="kog.correlations.pdf",width=12,height=6)
ggarrange(plot.DVelVelUmb,plot.VelUmb_UmbEyed,ncol=2,nrow=1)
dev.off()

corrPlot(x="D.2.Velyger",y="Velyger.2.Umbo",kogtable)
corrPlot(x="Velyger.2.Umbo",y="Umbo.2.Eyed",kogtable)
corrPlot(x="D.2.Velyger",y="Umbo.2.Eyed",kogtable)


# Plot interest heatmap related genes ----
## import names
library(pheatmap)
library(timeOmics)
library(PCAtools)

#data.names<-read.table("vsd_assay_devlarve_comp.txt",header=T,sep="\t",check.names = F)
data.names_t<-as.data.frame(t(vsd.assay.com))
data.names_t[1:5,1:5]
coldata[1:5,1:2]

data.combine <-merge(data.names_t,coldata,by=0)
rownames(data.combine)<-data.combine$Row.names
data.combine$Row.names <- NULL
data.combine$Timing <- NULL

## select only biomineralization genes
biomin<-read.table("gomwu/biomineralization_genes.txt",check.names=F)
data.biomin<-data.combine[,colnames(data.combine) %in% biomin$V1 ]
data.biomin <- data.biomin[,order(colnames(data.biomin))]

library(dplyr)
data.biomin %>%
  summarize_if(is.numeric, sum, na.rm=TRUE)

##PCA
## total
mat.tot<-(as.matrix(t(data.biomin)))
p <- pca(mat.tot, metadata = coldata)
biplot(p,colby = "Condition")
plotloadings(p, labSize = 2)
loading.biomin<-as.data.frame(p$loadings)
loading.biomin[rownames(loading.biomin) == "evm.TU.scaffold2500size167958.1",]

#gChange gene names for protein
gene.names<-read.table("biomin_annot.txt",header=F,check.names=F)
rownames(gene.names)<-gene.names$V1
gene.names <- gene.names[order(rownames(gene.names)),]
colnames(data.biomin) = make.names(gene.names$V2,unique=T)

# Remove low variance
#biomin.sub<-remove.low.cv(data.biomin, cutoff = 0.05)
#colnames(biomin.sub)
#write.table(colnames(biomin.sub),file="list_dopa_devLarv.txt",quote=F,col.names = F,row.names = F)
#change name
coldata2<-coldata
coldata2$Group <- as.character(coldata2$Condition)
coldata2$Group
coldata2$Group[coldata2$Group == "umbo"] <- "Umbo"
coldata2$Group[coldata2$Group == "oeillee"] <- "Eye-spot"
coldata2$Group[coldata2$Group == "D"] <- "D-shape"
coldata2$Group[coldata2$Group == "velyger"] <- "Velyger"
coldata2$Group

aka2 = data.frame(Stages = factor(rep(c("Umbo","Eye-spot","D-shape","Velyger"), each=3)))
rownames(aka2)<-rownames(data.biomin)
aka2
aka3 = list(Stages = c(Umbo = "orange2", 'Eye-spot'="red", 'D-shape'="turquoise",Velyger="dodgerblue2"))
aka3[1]

library(pheatmap)

set.seed(49)
#pdf(file="immune_test.pdf",height=20,width=5)

pheatmap(as.matrix(t(data.biomin)),
          annotation_col = aka2,
          annotation_colors = aka3[1],
          annotation_legend = T,
          show_colnames = T, show_rownames = T, cluster_rows = T,
          cluster_cols = T, legend = T,
         clustering_distance_rows = "euclidean",
         cellwidth = 10,border_color = FALSE)



library("assertthat")
library("scales")
library("WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

logcpm_test<-read.table("vsd_assay_devlarve.txt", header=T,check.names = FALSE)
head(logcpm_test)

temp<-read.table("data.trait.txt",sep="\t", header=T,na.strings="NA")
temp$Ind<-NULL
vector<-temp$Ind
vector$Ind<-NULL
logcpm_test<-logcpm_test[, names(logcpm_test) %in% vector]

# Take a quick look at what is in the data set:
dim(logcpm_test)


datExpr0 = as.data.frame(t(logcpm_test))
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 220, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 220, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nSamples
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)


# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0.05 # TEST for 0.05
table(KeepGenes)
datExpr=datExpr[, KeepGenes]

name_datExpr <-colnames(datExpr)

allTraits = read.table("data.trait.txt",sep="\t", header=T,na.strings="NA");
names(allTraits)

summary(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Ind);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
str(datTraits)

collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits,signed= FALSE);
# Plot the sample dendrogram and the colors underneath.
#pdf("dendo_heatmap_subset.pdf",width=12,height=9)
par(mar=c(1, 10, 1, 1))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

dev.off()

save(datExpr, datTraits, file = "dataInput_wgcna_global.Rda")

#setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "dataInput_wgcna_global.Rda");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 15, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")
View(sft$fitIndices)
# Plot the results:
#load("sft_signed_cell_culture_clam.Rda")
#tiff(file = "wgcna_mean_connectivity_cell.tiff", width = 20, height = 20, units = "cm", res = 300)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


save(sft,file="sft_signed_global.Rda")


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#setting is important, do not omit.
# Load the data saved in the first part
lnames = load(file = "dataInput_wgcna_global.Rda");
#The variable lnames contains the names of loaded variables.
lnames


load("dataInput_wgcna_global.Rda")
load("sft_signed_global.Rda")
message("data loaded")

softPower = 14; #low% but recommanded
adjacency = adjacency(datExpr, power = softPower,type="signed");
save(adjacency,file="adjacency_global.Rda")
#load("adjacency_global.Rda")
TOM = TOMsimilarity(adjacency,TOMType = "signed");
save(TOM,file="TOM_global.Rda")
#load("TOM_global.Rda")
dissTOM = 1-TOM

#save(dissTOM,file="disTOM.Rda")
message("dissimilarity done")
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf("genetree_subset.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf("genetree_dyn_subset.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf("cluster_modules_subset.pdf",width=8,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
dev.off()

pdf(file = "geneDendro-3_subset.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels,dynamicColors, dynamicMods, moduleColors, geneTree, file = "network_global.Rda")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "dataInput_wgcna_global.Rda");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "network_global.Rda");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# load SFT  ### not made in the proper package
lnames=load(file="sft_signed_global.Rda")
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(moduleTraitCor)
head(moduleTraitPvalue)

sizeGrWindow(25,15)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");



dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textAdj = c(0.5, 0.5),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# reformat matrix
par(mar = c(1,1, 1 ,1));
library(corrplot)
moduleTraitPvalue[moduleTraitPvalue > 0.01]<-0
moduleTraitCor[moduleTraitCor > -0.32 & moduleTraitCor < 0.32]<-0
df_module<-as.data.frame(moduleTraitPvalue)
df_module<-as.data.frame(moduleTraitCor)
new_mat <- as.data.frame(moduleTraitCor[ rownames(moduleTraitCor) %in% c("MEcyan","MEdarkgrey","MElightgreen","MEmidnightblue","MEturquoise",
                                                                         "MEred","MEblue","MEbrown","MEtan","MEdarkgreen",
                                                                         "MEblack","MEdarkred","MEpink","MEgrey60","MEyellow"),])
row.names(new_mat)<-gsub("ME", "", row.names(new_mat))
row.names(new_mat)<-c("cyan (341)","darkgrey (85)","lightgreen (159)","midnightblue (291)","turquoise (9,863)",
                      "red (1,239)","blue (6,183)","brown (3,659)","tan (493)","darkgreen (105)",
                      "black (1,083)","darkred (113)","pink (770)","grey60 (181)","yellow (2,199)")
mat<-as.matrix(new_mat)
mat<-mat[,c(1,2,3)]

#pdf(file = "correlation_matrix.global.pdf")
par(mar=c(2,5,6,5))

corrplothost<-corrplot(mat, method = "number",tl.col = "black",cl.pos = "n", col= colorRampPalette(c("darkblue","white", "darkred"))(100))
corrplothost
dev.off()


## Stage timing  ----
#### Get info
Timing= as.data.frame(datTraits$Timing);
names(Timing) = "Timing"
str(Timing)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
str(geneModuleMembership)


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
str(MMPvalue)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Timing, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Timing), sep="");
names(GSPvalue) = paste("p.GS.", names(Timing), sep="");

## plot correlation
### significance by module

module = "darkmagenta"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(file = "module_darkmagenta.pdf", width = 6, height = 6)
plot<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership ", module, "module"),
                         ylab = "Gene significance for Timing",
                         main = NULL, cex.lab = 1.2, cex.axis = 1.2, col = module)
plot + text(0.7, 0.78, "abs(cor) = 0.68\np-value < 0.001",cex = 1.5,font=2)
dev.off()



# export geneinfo
annot = read.csv(file = "annotation_genome_pmarg.110621.csv",sep=";",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$Sp[probes2annot],
                       LocusLinkID = annot$ID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for temperature
modOrder = order(-abs(cor(MEs, Timing, use = "p")));

# get relevant modules
# Get the corresponding Locuis Link IDs
allLLIDs = annot$name[probes2annot];
# $ Choose interesting modules
intModules = c("cyan","darkgrey","lightgreen","midnightblue","turquoise",
               "red","blue","brown","tan","darkgreen","black","darkred",
               "pink","grey60","yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("module_global_", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote=FALSE)
}

# A
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Timing));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneinfo_Timing.csv")

## Plot MEs
dataME<-read.table("module.eigenvalue.txt",header=T)
themeL<-theme(panel.border=element_rect(fill=NA),
             panel.background=element_rect(fill="grey94"),
             axis.text.x=element_text(colour="black",size=12),
             axis.text.y=element_text(colour="black",size=12),
             axis.title.x=element_blank(),
             axis.title.y=element_text(colour="black",size=12),
             axis.ticks=element_line(colour="black"),
             plot.title = element_text(hjust=0.1,vjust=-15,face="bold",size=14),
             aspect.ratio=1)
themeR<-theme(panel.border=element_rect(fill=NA),
              panel.background=element_rect(fill="grey94"),
              axis.text.x=element_text(colour="black",size=12),
              axis.text.y=element_text(colour="black",size=12),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour="black",size=12),
              axis.ticks=element_line(colour="black"),
             plot.title = element_text(hjust=0.9,vjust=-15,face="bold",size=14),
              aspect.ratio=1)
themeM<-theme(panel.border=element_rect(fill=NA),
              panel.background=element_rect(fill="grey94"),
              axis.text.x=element_text(colour="black",size=12),
              axis.text.y=element_text(colour="black",size=12),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour="black",size=12),
              axis.ticks=element_line(colour="black"),
              plot.title = element_text(hjust=0.1,vjust=-15,face="bold",size=14),
              aspect.ratio=1)

themeD<-theme(panel.border=element_rect(fill=NA),
              panel.background=element_rect(fill="grey94"),
              axis.text.x=element_text(colour="black",size=12),
              axis.text.y=element_text(colour="black",size=12),
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour="black",size=12),
              axis.ticks=element_line(colour="black"),
              plot.title = element_text(hjust=0.1,vjust=-50,face="bold",size=14),
              aspect.ratio=1)
dataME$Stage<-factor(dataME$Stage,levels=c("D","Velyger","Umbo","Eyed"))

sien<-ggplot(dataME,aes(x=Stage,y=MEsienna3))+geom_boxplot(fill="Sienna3") + themeL + ylab("Module eigenvalue")+ggtitle("Sienna3")
blu<-ggplot(dataME,aes(x=Stage,y=MEblue))+geom_boxplot(fill="Blue") + themeL + ggtitle("Blue")+ ylab("Module eigenvalue")
brow<-ggplot(dataME,aes(x=Stage,y=MEbrown))+geom_boxplot(fill="Brown") + themeR + ggtitle("Brown")+ ylab("Module eigenvalue")
dtur<-ggplot(dataME,aes(x=Stage,y=MEdarkturquoise))+geom_boxplot(fill="Darkturquoise") + themeL + ggtitle("Darkturquoise")+ ylab("Module eigenvalue")
mag<-ggplot(dataME,aes(x=Stage,y=MEmagenta))+geom_boxplot(fill="Magenta") + themeL + ggtitle("Magenta")+ ylab("Module eigenvalue")
rblu<-ggplot(dataME,aes(x=Stage,y=MEroyalblue))+geom_boxplot(fill="Royalblue") + themeL + ggtitle("Royalblue")+ ylab("Module eigenvalue")
ygree<-ggplot(dataME,aes(x=Stage,y=MEyellowgreen))+geom_boxplot(fill="Yellowgreen") + themeD + ggtitle("Yellowgreen")+ ylab("Module eigenvalue")
dmag<-ggplot(dataME,aes(x=Stage,y=MEdarkmagenta))+geom_boxplot(fill="Darkmagenta") + themeR + ggtitle("Darkmagenta")+ ylab("Module eigenvalue")
lcy<-ggplot(dataME,aes(x=Stage,y=MElightcyan))+geom_boxplot(fill="Lightcyan") + themeM + ggtitle("Lightcyan")+ ylab("Module eigenvalue")
lgre<-ggplot(dataME,aes(x=Stage,y=MElightgreen))+geom_boxplot(fill="Lightgreen") + themeR + ggtitle("Lightgreen")+ ylab("Module eigenvalue")
drgy<-ggplot(dataME,aes(x=Stage,y=MEdarkgrey))+geom_boxplot(fill="Darkgrey") + themeR + ggtitle("Darkgrey")+ ylab("Module eigenvalue")


library(ggpubr)
#cluster1
dmag
sien
lcy
ygree

pdf(file="cluster1.wgcna.pdf)")
ggarrange(dmag,sien,lcy,ygree,ncol=2,nrow=2)
dev.off()

#cluster2
blu
mag
rblu

pdf(file="cluster2.wgcna.pdf)")
ggarrange(blu,mag,rblu,ncol=2,nrow=2)
dev.off()

#cluster3
dtur
lgre
brow
drgy

pdf(file="cluster3.wgcna.pdf)")
ggarrange(dtur,lgre,brow,drgy,ncol=2,nrow=2)
dev.off()
