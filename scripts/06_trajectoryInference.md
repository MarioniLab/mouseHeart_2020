---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '27 November, 2019'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: true
    toc_float: 
      collapsed: false
---



### Trajectory inference

We begin with the normalised and batch-corrected data from `02_batchCorrections.Rmd`.


```r
## normalised, batch corrected counts
sce.corr <- readRDS(paste0(dir, "data/sce_goodQual.NORM.batchCorrected.Rds"))

## HVGs
hvgs <- read.table(paste0(dir, "results/HVGs_minMean1_FDR0.05.tsv"), stringsAsFactors = FALSE)
hvgs <- hvgs$V1
```

We have clustered the cells into different cell types, from all three germ layers, and used the expression of marker genes to annotate them.


```r
## clustering results
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))

plot(reducedDim(sce)$x, reducedDim(sce)$y, pch=16, col=sce$clusterCol, xlab="", ylab="", axes=FALSE)
box(bty="l")
legend("topright", legend = names(cols), col=cols, cex=0.75, pch=16)
text(7, 14, labels = "ectoderm", col="firebrick", cex=0.75, font=2)
text(9, -3, labels = "endoderm", col="steelblue4", cex=0.75, font=2)
text(-3.5, 5, labels = "cardiac mesoderm", col="darkolivegreen", cex=0.75, font=2)
text(-9, 8, labels = "endothelium", col="wheat3", cex=0.75, font=2)
text(6, 1.5, labels = "blood", col="burlywood3", cex=0.75, font=2)
```

![](06_trajectoryInference_files/figure-html/clusters-1.png)<!-- -->

Both the turquoise (Me7-8) and green (Me5) clusters from the cardiac mesoderm subpopulation are *progenitor-like*, while the blue (Me3) subcluster contains the most *mature cardiomyocytes*. The orange (Me4) and yellow (Me6) clusters are intermediates between the cardiomyocytes and progenitor cells.

Interestingly, expression of markers of the first heart field (*Hand1*) are predominant in the green and orange clusters, while second heart field genes (*Isl1*) are biased towards the turquoise and yellow clusters. 


```r
plotGeneOnUMAP <- function(umap=umap, data=data, gene=gene){
  df <- data.frame(x=umap$x, y=umap$y, expr=logcounts(data)[which(rowData(data)$gene == gene),])
  df <- df[order(df$expr),]
  p <- ggplot(df, aes(x,y)) + geom_point(aes(colour=expr), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(gene) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5) + guides(colour = guide_colorbar(title.position = "bottom"))
  return(p)
}

plots <- list()
plots[[1]] <- plotGeneOnUMAP(umap = reducedDim(sce), data = sce, gene = "Hand1")
plots[[2]] <- plotGeneOnUMAP(umap = reducedDim(sce), data = sce, gene = "Isl1")

ggarrange(plotlist = plots, ncol = 2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/plot_markers-1.png)<!-- -->

We will use a diffusion map approach to study the differentiation trajectories from the progenitor populations towards cardiomyocytes.


```r
## select only cells in cardiac clusters
sce.cardiac <- sce[,sce$clusterAnn %in% paste0("Me", c(3:8))]
sce.corr.cardiac <- sce.corr[,colnames(sce.cardiac)]
# plot(reducedDim(sce.cardiac)$x, reducedDim(sce.cardiac)$y, pch=16, col=sce.cardiac$clusterCol)
```

Using all cells in the cardiac mesoderm clusters, the diffusion map embedding arranges cells in a triangular structure, with the two progenitor populations and the cardiomyocytes at the tips; the yellow and orange clusters are intermediates, with the orange stretching out to both progenitor types.


```r
## batch corrected data for HVGs
data <- assay(sce.corr.cardiac)[hvgs,]

## diffusion maps
pca <- prcomp(t(data))
set.seed(123)
diffMap <- DiffusionMap(pca$x[,1:50], sigma="local", distance="euclidean")
# qplot(y = eigenvalues(diffMap)) + theme_minimal() + labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')

dm <- as.data.frame(diffMap@eigenvectors[,1:5])
row.names(dm) <- colnames(data)

## diffusion pseudotime
dpt <- DPT(diffMap)

plots <- list()
plots[[1]] <- ggplot(dm, aes(DC1, DC2, colour=sce.cardiac$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th
plots[[2]] <- plot(dpt, col_by="branch", root=2, path=c(1,3)) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/DM-1.png)<!-- -->

And using the diffusion pseudotime we can define branches. The turquoise subclusters (Me7 and Me8) form one branch, and so does the green cluster (Me5); the yellow, orange and blue clusters are all merged into the third branch.

This suggests that the two progenitor types converge to a common differentiation trajectory to give rise to cardiomyocytes.

#### Me5 to cardiomyocyte trajectory

To investigate in more detail the differentiation dynamics of the Me5 progenitors to mature cardiomyocytes, we now compute the diffusion map excluding the turquoise cells, which comprise a distinct branch.

Once again, cells are arranged in a triangular structure, with the progenitors from Me5 and the mature cardiomyocytes at opposite ends; the yellow cluster segregates towards the other tip.


```r
## exclude Me7 and Me8
sce.Me5 <- sce[,sce$clusterAnn %in% paste0("Me", c(3:6))]
sce.corr.Me5 <- sce.corr[,colnames(sce.Me5)]

## diffusion map
data <- assay(sce.corr.Me5)[hvgs,]
pca <- prcomp(t(data))
set.seed(123)
diffMap.me5 <- DiffusionMap(pca$x[,1:50], sigma="local", distance="euclidean")
# qplot(y = eigenvalues(diffMap.me5)) + theme_minimal() + labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')

dm.me5 <- as.data.frame(diffMap.me5@eigenvectors[,1:5])
row.names(dm.me5) <- colnames(data)

## diffusion pseudotime
dpt.me5 <- DPT(diffMap.me5)

plots <- list()
plots[[1]] <- ggplot(dm.me5, aes(DC1, DC2, colour=sce.Me5$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th
plots[[2]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,3)) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/me5-1.png)<!-- -->

And defining branches based on diffusion pseudotime separates the cardiomyocytes from the yellow cluster and the orange + green progenitors. 

The cells from branch 1 (orange + green) can be further split with two small subclusters at the tip branching out. In the third diffusion component is evident that the cells in purple differentiate towards red, which in turn link to the turquoise branch of cardiomyocytes. But the cells in the blue branch deviate in a different trajectory, that doesn't progress towards cardiomyocytes.


```r
plots <- list()
plots[[1]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,1), divide=1, dcs=c(1,2)) + th
plots[[2]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,3), divide=1, dcs=c(1,3)) + th
ggarrange(plotlist = plots, ncol = 2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/branch2-1.png)<!-- -->

The cells that deviate from the trajectory towards cardiomyocytes (blue branch) express higher levels of several keratins and *Hoxb6*.


```r
## get cells in the three branches of of the orange+green clusters
branch1.a <- row.names(dpt.me5@branch)[which(dpt.me5@branch[,2] == 5)]
branch1.b <- row.names(dpt.me5@branch)[which(dpt.me5@branch[,2] == 4)]
branch1.c <- row.names(dpt.me5@branch)[which(dpt.me5@branch[,2] == 6)] ## away from Me3

## get cells in branches 4:6 and define branches as clusters
sce.branched <- sce[,c(branch1.a, branch1.b, branch1.c)]
clst <- rep(c("A","B","C"), c(length(branch1.a), length(branch1.b), length(branch1.c)))
names(clst) <- c(branch1.a, branch1.b, branch1.c)
stopifnot(identical(colnames(sce.branched), names(clst)))

## find genes expressed higher in the branch that doesn't go to Me3
keep <- rowMeans(logcounts(sce.branched)) > 0.1
de <- findMarkers(sce.branched, clst, block=sce.branched$batch, direction="up", subset.row=keep, pval.type="all")

## plot the top DE genes
tmp <- logcounts(sce.branched)[row.names(de[[3]])[1:20],]
tmp <- t(apply(tmp, 1, function(x) (x-mean(x))/sd(x) )) # z-score 
row.names(tmp) <- rowData(sce.cardiac)[row.names(tmp),]$gene

## clusters (branches) to colours
c <- c("indianred1", "plum3", "steelblue3")
names(c) <- c("A","B","C")

heatmap.2(tmp, trace="none", ColSideColors = c[clst], col=colorRampPalette(colors = c("steelblue","white","indianred")), key.xlab="z-score")
```

![](06_trajectoryInference_files/figure-html/branch_de-1.png)<!-- -->

#### Me7/8 to cardiomyocyte trajectory

Now do the same to define the alternative trajectory from the turquoise progenitors towards cardiomyocytes.

Again, cells are arranged in a triangular structure, with the Me7-8 progenitors and the cardiomyocytes on opposite ends; and in this case, the cells from the orange cluster are in the other tip. 


```r
## exclude Me5
sce.Me7 <- sce[,sce$clusterAnn %in% paste0("Me", c(3:4,6:8))]
sce.corr.Me7 <- sce.corr[,colnames(sce.Me7)]

## diffusion map
data <- assay(sce.corr.Me7)[hvgs,]
pca <- prcomp(t(data))

set.seed(123)
diffMap.me7 <- DiffusionMap(pca$x[,1:50], sigma="local", distance="euclidean")
# qplot(y = eigenvalues(diffMap.me7)) + theme_minimal() + labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')

dm.me7 <- as.data.frame(diffMap.me7@eigenvectors[,1:5])
row.names(dm.me7) <- colnames(data)

## diffusion pseudotime
dpt.me7 <- DPT(diffMap.me7)

plots <- list()
plots[[1]] <- ggplot(dm.me7, aes(DC1, DC2, colour=sce.Me7$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th
plots[[2]] <- plot(dpt.me7, col_by="branch", root=1, path=c(2,3)) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/me7-1.png)<!-- -->

Consistently, there is a branching point that separates the most distinct orange cells; but some orange cells are part of the branch containing the turquoise and yellow progenitors.

Further splitting branch 2, which contains the turquoise progenitors, also separates a small number of cells at the tip into two branches, suggesting further heterogeneity.


```r
plots <- list()
plots[[1]] <- plot(dpt.me7, col_by="branch", root=1, path=c(2,3), divide=2, dcs=c(1,2)) + th
plots[[2]] <- plot(dpt.me7, col_by="branch", root=1, path=c(2,3), divide=2, dcs=c(1,4)) + th
ggarrange(plotlist = plots, ncol = 2, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/branch2_me7-1.png)<!-- -->

#### Recap

All together, we can identify two different trajectories leading to mature cardiomyocytes:

- A path starting in the green cells (Me5), that goes through the orange cluster (Me4).
    - The cells from the yellow cluster (Me6) form a different branch from this path.
    - A small number of cells at the tip of Me5 branch out towards a different trajectory, away from cardiomyocytes.
- Another path starting in the turquoise cells from Me7 and Me8, that goes through the yellow (Me6) and some cells from the orange (Me4) clusters.

#### Dynamic gene expression

Based on the two trajectories computed above, we can define the pseudotime progression of cells in each, considering only the branch going to mature cardiomyocytes, and ignoring the cells branching off from the yellow and orange clusters in the Me5 and Me7 trajectories respectively.

Many of the mature cardiomyocytes are shared between the two trajectories, but the intermediate cells linking them to the progenitors are different.


```r
## define pseudotime for the Me5->Me3 trajectory
## consider only cells in branches 1 and 2
cells.me5 <- row.names(dpt.me5@branch)[which(dpt.me5@branch[,1] %in% 1:2)]

# find which tip corresponds to Me5
# which(dpt.me5@tips[,1])
# plot(dpt.me5, col_by="DPT1135")
diffPseudotime_me5 <- dpt.me5$DPT1135
names(diffPseudotime_me5) <- colnames(sce.corr.Me5)
diffPseudotime_me5 <- diffPseudotime_me5[cells.me5]


## define pseudotime for the Me7->Me3 trajectory
## consider only cells in branches 1 and 2
cells.me7 <- row.names(dpt.me7@branch)[which(dpt.me7@branch[,1] %in% 1:2)]

# find which tip corresponds to Me7
# which(dpt.me7@tips[,1])
# plot(dpt.me7, col_by="DPT1161")
diffPseudotime_me7 <- dpt.me7$DPT1161
names(diffPseudotime_me7) <- colnames(sce.corr.Me7)
diffPseudotime_me7 <- diffPseudotime_me7[cells.me7]

plots <- list()
plots[[1]] <- ggplot(dm, aes(DC1, DC2, colour=sce.cardiac$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + ggtitle("Clusters") + theme(plot.title = element_text(face="bold", hjust = 0.5)) + th
plots[[2]] <- ggplot(dm, aes(DC1, DC2)) + geom_point(data=dm[setdiff(row.names(dm), names(diffPseudotime_me5)),], col="grey", size=1) + geom_point(data=dm[names(diffPseudotime_me5),], aes(col = diffPseudotime_me5)) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Greens")[-c(1,9)]) + labs(colour="DPT") + ggtitle("Me5 -> Me3 trajectory") + theme(plot.title = element_text(face="bold", hjust = 0.5)) + th
plots[[3]] <- ggplot(dm, aes(DC1, DC2)) + geom_point(data=dm[setdiff(row.names(dm), names(diffPseudotime_me7)),], col="grey", size=1) + geom_point(data=dm[names(diffPseudotime_me7),], aes(col = diffPseudotime_me7)) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Blues")[-c(1,9)]) + labs(colour="DPT") + ggtitle("Me7 -> Me3 trajectory") + theme(plot.title = element_text(face="bold", hjust = 0.5)) + th

ggarrange(plotlist = plots, ncol = 3, nrow = 1)
```

![](06_trajectoryInference_files/figure-html/pseudotime-1.png)<!-- -->

Based on these cell orderings, we can identify genes that are expressed dynamically along the trajectories. For this, we follow the method used in *Ibarra-Soria, Wajid et al., Nat Cell Biol, 2018*.


```r
## Me5->Me3 trajectory
diffPseudotime_me5 <- diffPseudotime_me5[order(diffPseudotime_me5)]
data <- logcounts(sce)[,names(diffPseudotime_me5)]
data <- data[rowSums(data)>0,] ## remove non-expressed genes

## only test genes expressed at decent levels
keep <- rowMeans(data) > 0.1
data <- data[keep,]

## fit polynomials of degree 0 and 2 with local regression, and compute the difference in their AIC
## we use linear insted of logistic regression (as in the Nat Cell Bio paper) because it is SmartSeq2 data instead of 10X
smooth = 1
deltaAIC_me5 <- sapply(row.names(data), function(x){
  df <- data.frame(diffPseudotime_me5, data[x,])
  delta <- aic(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)['aic'] - aic(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)['aic']
  return(delta)
})
names(deltaAIC_me5) <- row.names(data)
deltaAIC_me5 <- deltaAIC_me5[order(deltaAIC_me5)]

# plot(density(deltaAIC_me5))
# length(deltaAIC_me5[deltaAIC_me5 < -100])  # 2316
# length(deltaAIC_me5[deltaAIC_me5 < -200])  # 1143
# length(deltaAIC_me5[deltaAIC_me5 < -300])  # 726

dynGenes_me5 <- deltaAIC_me5[deltaAIC_me5 < -300]

# for(gene in names(deltaAIC_me5[701:710])){
#   df <- data.frame(diffPseudotime_me5, data[gene,])
#   fit0 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
# 
#   plot(df[,1], df[,2], pch=16, col=rgb(0,0,0,0.3), bty="l", xlab="Me5 ---> Me3", ylab=expression('log'[2]*' expression'), main=rowData(sce)[gene,1])
#   lines(diffPseudotime_me5, predict(fit0, diffPseudotime_me5), lwd=3, col="steelblue2")
#   lines(diffPseudotime_me5, predict(fit2, diffPseudotime_me5), lwd=3, col="indianred2")
# }
```

We identify 726 genes dynamically expressed from the Me5 green progenitors down to cardiomyocytes.


```r
## Me7->Me3 trajectory
diffPseudotime_me7 <- diffPseudotime_me7[order(diffPseudotime_me7)]
data <- logcounts(sce)[,names(diffPseudotime_me7)]
data <- data[rowSums(data)>0,] ## remove non-expressed genes

## only test genes expressed at decent levels
keep <- rowMeans(data) > 0.1
data <- data[keep,]

## fit polynomials of degree 0 and 2 with local regression, and compute the difference in their AIC
## we use linear insted of logistic regression (as in the Nat Cell Bio paper) because it is SmartSeq2 data instead of 10X
smooth = 1
deltaAIC_me7 <- sapply(row.names(data), function(x){
  df <- data.frame(diffPseudotime_me7, data[x,])
  delta <- aic(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)['aic'] - aic(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)['aic']
  return(delta)
})
names(deltaAIC_me7) <- row.names(data)
deltaAIC_me7 <- deltaAIC_me7[order(deltaAIC_me7)]

# plot(density(deltaAIC_me7))
# length(deltaAIC_me7[deltaAIC_me7 < -100])  # 3281
# length(deltaAIC_me7[deltaAIC_me7 < -200])  # 1747
# length(deltaAIC_me7[deltaAIC_me7 < -300])  # 1110

dynGenes_me7 <- deltaAIC_me7[deltaAIC_me7 < -300]

# gene  <- names(deltaAIC_me7[472])
# df <- data.frame(diffPseudotime_me7, data[gene,]); colnames(df) <- c("DPT", "expr")
# fit0 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)
# fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
# ggplot(df, aes(DPT, expr)) + geom_point(colour=sce[,names(diffPseudotime_me7)]$clusterCol, alpha=0.5) + geom_line(aes(diffPseudotime_me7, predict(fit0, diffPseudotime_me7)), lwd=1, col="steelblue2") + geom_line(aes(diffPseudotime_me7, predict(fit2, diffPseudotime_me7)), lwd=1, col="indianred2") + th
```

And 1110 genes dynamic genes from the Me7 turquoise progenitors.


```r
## diffusion map
saveRDS(diffMap, paste0(dir, "results/diffusionMap_cardiacMesoderm.Rds"))
saveRDS(diffMap.me5, paste0(dir, "results/diffusionMap_Me5.Rds"))
saveRDS(diffMap.me7, paste0(dir, "results/diffusionMap_Me7.Rds"))

## DPT
saveRDS(dpt, paste0(dir, "results/diffusionPseudotime_cardiacMesoderm.Rds"))
saveRDS(dpt.me5, paste0(dir, "results/diffusionPseudotime_Me5.Rds"))
saveRDS(dpt.me7, paste0(dir, "results/diffusionPseudotime_Me7.Rds"))

## diffusion pseudotime for each trajectory
write.table(diffPseudotime_me5, file=paste0(dir, "results/diffusionPseudotime_Me5.tsv"), quote = FALSE, sep="\t", col.names = FALSE)
write.table(diffPseudotime_me7, file=paste0(dir, "results/diffusionPseudotime_Me7.tsv"), quote = FALSE, sep="\t", col.names = FALSE)

## dynamic genes.
write.table(deltaAIC_me5, file = paste0(dir, "results/deltaAIC_Me5.tsv"), quote = FALSE, sep="\t", col.names = FALSE)
write.table(deltaAIC_me7, file = paste0(dir, "results/deltaAIC_Me7.tsv"), quote = FALSE, sep="\t", col.names = FALSE)
```




```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] locfit_1.5-9.1              gplots_3.0.1.1             
##  [3] RColorBrewer_1.1-2          ggpubr_0.2.3               
##  [5] magrittr_1.5                ggplot2_3.2.1              
##  [7] destiny_3.0.0               scran_1.14.1               
##  [9] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
## [11] DelayedArray_0.12.0         BiocParallel_1.20.0        
## [13] matrixStats_0.55.0          Biobase_2.46.0             
## [15] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [17] IRanges_2.20.0              S4Vectors_0.24.0           
## [19] BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0         colorspace_1.4-1        
##   [3] ggsignif_0.6.0           RcppEigen_0.3.3.5.0     
##   [5] class_7.3-15             rio_0.5.16              
##   [7] XVector_0.26.0           RcppHNSW_0.2.0          
##   [9] BiocNeighbors_1.4.1      rstudioapi_0.10         
##  [11] proxy_0.4-23             hexbin_1.28.0           
##  [13] RSpectra_0.15-0          ranger_0.11.2           
##  [15] codetools_0.2-16         robustbase_0.93-5       
##  [17] knitr_1.25               scater_1.14.1           
##  [19] zeallot_0.1.0            compiler_3.6.1          
##  [21] ggplot.multistats_1.0.0  dqrng_0.2.1             
##  [23] backports_1.1.5          assertthat_0.2.1        
##  [25] Matrix_1.2-17            lazyeval_0.2.2          
##  [27] limma_3.42.0             BiocSingular_1.2.0      
##  [29] htmltools_0.4.0          tools_3.6.1             
##  [31] rsvd_1.0.2               igraph_1.2.4.1          
##  [33] gtable_0.3.0             glue_1.3.1              
##  [35] GenomeInfoDbData_1.2.2   dplyr_0.8.3             
##  [37] ggthemes_4.2.0           Rcpp_1.0.2              
##  [39] carData_3.0-2            cellranger_1.1.0        
##  [41] vctrs_0.2.0              gdata_2.18.0            
##  [43] DelayedMatrixStats_1.8.0 lmtest_0.9-37           
##  [45] xfun_0.10                laeken_0.5.0            
##  [47] stringr_1.4.0            openxlsx_4.1.3          
##  [49] lifecycle_0.1.0          irlba_2.3.3             
##  [51] gtools_3.8.1             statmod_1.4.32          
##  [53] edgeR_3.28.0             DEoptimR_1.0-8          
##  [55] zlibbioc_1.32.0          MASS_7.3-51.4           
##  [57] zoo_1.8-6                scales_1.0.0            
##  [59] VIM_4.8.0                pcaMethods_1.78.0       
##  [61] hms_0.5.2                yaml_2.2.0              
##  [63] curl_4.2                 gridExtra_2.3           
##  [65] stringi_1.4.3            knn.covertree_1.0       
##  [67] e1071_1.7-2              TTR_0.23-5              
##  [69] caTools_1.17.1.2         boot_1.3-23             
##  [71] zip_2.0.4                rlang_0.4.1             
##  [73] pkgconfig_2.0.3          bitops_1.0-6            
##  [75] evaluate_0.14            lattice_0.20-38         
##  [77] purrr_0.3.3              labeling_0.3            
##  [79] cowplot_1.0.0            tidyselect_0.2.5        
##  [81] R6_2.4.0                 withr_2.1.2             
##  [83] pillar_1.4.2             haven_2.2.0             
##  [85] foreign_0.8-72           xts_0.11-2              
##  [87] scatterplot3d_0.3-41     abind_1.4-5             
##  [89] RCurl_1.95-4.12          sp_1.3-2                
##  [91] nnet_7.3-12              tibble_2.1.3            
##  [93] crayon_1.3.4             car_3.0-4               
##  [95] KernSmooth_2.23-16       rmarkdown_1.16          
##  [97] viridis_0.5.1            grid_3.6.1              
##  [99] readxl_1.3.1             data.table_1.12.6       
## [101] forcats_0.4.0            vcd_1.4-4               
## [103] digest_0.6.22            tidyr_1.0.0             
## [105] munsell_0.5.0            beeswarm_0.2.3          
## [107] viridisLite_0.3.0        smoother_1.1            
## [109] vipor_0.4.5
```

