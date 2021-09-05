
# This scripts performs  hierarchical clustering
# on the homologous genes in E. coli strain K12 that 
# correspond to the genes in M. tuberculosis

# Package requirements & installation
########################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")


library("BiocManager")
BiocManager::install("org.EcK12.eg.db")
BiocManager::install("GOSemSim")
BiocManager::install("WGCNA")
BiocManager::install("cutreeStaticColor")
BiocManager::install()
BiocManager::available()

install.packages("fastcluster")
install.packages("moduleColor")
install.packages("dendextend")
install.packages("plotDendroAndColors")
installed.packages()
########################################

rm(list = ls())
getwd()
dev.off(dev.list()["RStudioGD"])
options(stringsAsFactors = FALSE)

library(GOSemSim)
library(WGCNA)
library(dynamicTreeCut)
library(cluster)


Ecoli_homolog_evalue001 <- c("aceE","acnA","frlD","alaS","argC","argG","aroC","asnB","alaC","aspS",
                             "can","cstA","cysS","dapB","dapE","dapF","def","folA","sucB","dnaE",
                             "dnaK","dnaN","dxr","emrY","eno","fadD","fadD","fabB","fbaA","fmt",
                             "folB","ftsK","ftsW","fumC","gmd","gcvP","gcvH","glcB","glgB","glmU",
                             "gltB","gpmA","pnp","ispB","guaB","gyrA","gyrB","hemA","hemB","hemC",
                             "hemL","hemH","hisI","hsdM","ilvC","infB","iscS","fabF","fabF","leuA",
                             "dmlA","leuS","lipA","aes","lysA","manA","rfbA","mdh","menD","menE",
                             "ubiE","metE","metG","yehL","baeS","cysS","murA","murC","murD","murG",
                             "ndh","nrdF","nrdH","nusA","otsA","ycjT","minD","cca","lnt","prfA",
                             "proB","proC","proS","purD","purL","purM","pyrB","pyrH","putA","rpoB",
                             "rpsM","secY","serC","sucC","thyA","topA","tpiA","trmD","pabA","trpS",
                             "tyrS","torR","rarA","ftsK","glpC","ytjC","ccmG","rbn","ftsK","mepH",
                             "argA","dnaX")

ecoli <- godata('org.EcK12.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)
sim <- mgeneSim(Ecoli_homolog_evalue001, semData=ecoli, measure="Wang",drop = NULL, combine="BMA", verbose=TRUE)
#write.csv(sim, file = "confusion_homolog.csv")
dissim <- 1-sim

########################################
require(fastcluster)
require(graphics)
require(dendextend)

#Hierarchical clustering
hc <- hclust(dist(dissim), method = "average")
hc.dend <- as.dendrogram(hc) # create dendrogram object
nleaves(hc.dend) # number of leaves in tree
nnodes(hc.dend)  # number of nodes (=leaves + joins) in tree


##Cut tree algorithm
cut_tree <- cutree(hc,h=max(hc$height)/2)
subcls <- sapply(unique(cut_tree),function(x)names(cut_tree[cut_tree==x]))
similarities <- mclusterSim(subcls,semData=ecoli, measure="Wang",drop = NULL, combine="BMA")

plot(hc, main = "genes clustering", sub="", xlab="")
plotDendroAndColors(hc, labels2colors(cut_tree),"Tree cut", dendroLabels = FALSE, hang = 0.05,
                    abHeight = max(hc$height)/2,
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )

##Dynamic cut tree
gene.names=hc$labels
#1 (Minimum cluster size = 5)
DynamicHybridADJ_minsize5 = cutreeDynamic(hc,distM= dissim,
                               method="hybrid", minClusterSize = 5,
                               pamRespectsDendro = FALSE)
table(DynamicHybridADJ_minsize5)
names(DynamicHybridADJ_minsize5) <- gene.names
subcls_minsize5 <- sapply(unique(DynamicHybridADJ_minsize5),function(x)names(DynamicHybridADJ_minsize5[DynamicHybridADJ_minsize5==x]))
similarities_minsize5 <- mclusterSim(subcls_minsize5,semData=ecoli, measure="Wang",drop = NULL, combine="BMA")


#2 (Minimum cluster size = 10)
DynamicHybridADJ_minsize10 = cutreeDynamic(hc,distM= dissim,
                                          method="hybrid", minClusterSize = 10,
                                          pamRespectsDendro = FALSE)
table(DynamicHybridADJ_minsize10)
names(DynamicHybridADJ_minsize10) <- gene.names
subcls_minsize10 <- sapply(unique(DynamicHybridADJ_minsize10),function(x)names(DynamicHybridADJ_minsize10[DynamicHybridADJ_minsize10==x]))
similarities_minsize10 <- mclusterSim(subcls_minsize10,semData=ecoli, measure="Wang",drop = NULL, combine="BMA")


# Plot results of all module detection methods together:
sizeGrWindow(40,8)
plotDendroAndColors(dendro = hc,
                    cbind(labels2colors(DynamicHybridADJ_minsize5),labels2colors(DynamicHybridADJ_minsize10)), 
                    c("min Cluster Size: 5","min Cluster Size: 10"),
                    dendroLabels = FALSE, 
                    main = "Gene dendrogram (complete) and module colors (Dynamic TreeCut)")



## Dendrogram and colors
library(gplots)     # heatmap.2
library(dendextend) # make and color dendrogram
library(colorspace) # diverge_hcl / rainbow_hcl / heat_hcl color palettes
library(RColorBrewer)

color.scheme <- rev(brewer.pal(10,"RdBu"))


par(cex.main=1)                   # adjust font size of titles
heatmap.2(sim, main = 'Gene similarity matrix',
          dendrogram = "none",        # no dendrogram for columns
          col = color.scheme,         # color pattern of the heatmap
          trace="none",              # hide trace
          cexRow=0.8, cexCol = 0.8,      # size of row / column labels
          xlab = "genes", ylab = "genes",
          Rowv = NULL, 
          Colv = NULL
)
dev.off(dev.list()["RStudioGD"])


## plot the heatmap with the dendrogram above ##
png("BPclusters.png",width = 20, height = 18, res = 300)
sizeGrWindow(40,8)
par(cex.main=1)                   # adjust font size of titles
heatmap.2(sim, main = 'Hierarchical clustering on Gene similarity matrix',
          dendrogram = "both",        # no dendrogram for columns
          col = color.scheme,         # color pattern of the heatmap
          trace="none",              # hide trace
          cexRow=0.8, cexCol = 0.8,      # size of row / column labels
          xlab = "genes", ylab = "genes",
          Rowv = ladderize(hc.dend), 
          Colv = ladderize(hc.dend), 
          ColSideColors=labels2colors(DynamicHybridADJ_minsize5),
          RowSideColors=labels2colors(DynamicHybridADJ_minsize5)
          #cbind(labels2colors(DynamicHybridADJ_minsize5))
)
dev.off(dev.list()["RStudioGD"])

