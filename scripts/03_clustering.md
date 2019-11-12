---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '12 November, 2019'
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



### Clustering

We start from the batch corrected, normalised counts produced in the `02_batchCorrection.Rmd` script.


```r
## normalised, batch corrected counts
sce.corr <- readRDS(paste0(dir, "data/sce_goodQual.NORM.batchCorrected.Rds"))

## HVGs
hvgs <- read.table(paste0(dir, "results/HVGs_minMean1_FDR0.05.tsv"), stringsAsFactors = FALSE)
hvgs <- hvgs$V1

## UMAP
umap <- read.table(paste0(dir, "results/umapCoords_corrected.tab"))
```

The next step is to cluster the cells to define the cellular diversity of the dataset. We use the distances between cells based on expression of HVGs on the batch corrected data. Clusters are defined through hierarchical clustering and a dynamic tree cut. 


```r
## use the distance between cells to identify clusters
dat <- assay(sce.corr)[hvgs,]
test.dist <- dist(t(dat))

## define clusters by hierarchical clustering and dynamic tree cut
test.clust <- hclust(test.dist, method="average")
cut <- cutreeDynamic(test.clust, distM=as.matrix(test.dist), minClusterSize=40, method="hybrid", deepSplit = 1, verbose = 0)
sce.corr$cluster <- cut

names(cut) <- colnames(dat)
write.table(cut, paste0(dir, "results/clusters_average_min40.tsv"), quote = FALSE, sep="\t", col.names = FALSE)

stopifnot(identical(names(cut), row.names(umap)))
umap$cluster <- cut
o <- order(umap$cluster)
plot(umap$x[o], umap$y[o], pch=16, cex=0.75, col=umap$cluster[o], xlab="UMAP - dim1", ylab="UMAP - dim2", bty="l")
```

![](03_clustering_files/figure-html/cluster-1.png)<!-- -->

This procedure results in 12 clusters.


```r
table(cut)[-1]
```

```
## cut
##   1   2   3   4   5   6   7   8   9  10  11  12 
## 713 514 405 355 260 221 194 184  89  65  59  45
```

Importantly, clusters are comprised of cells from all different batches. 


```r
tmp <- sce.corr[,-which(sce.corr$cluster==0)]
table(batch=tmp$batch, cluster=tmp$cluster)
```

```
##          cluster
## batch       1   2   3   4   5   6   7   8   9  10  11  12
##   batch_1  21   1   6   0   4   1   0   1   0   4   1   0
##   batch_2  83  57  34  23  16  10  26  16  17   9   0   1
##   batch_3  63  57  16  45  11  26  42  10  13   1   4   2
##   batch_4 149 133 110  71  43  33  41  19  12  18  20   4
##   batch_5 210  71 104  69  82  61  25   4   9  21   5  18
##   batch_6 130 101  62 106  59  56  23  34  12  10   6  16
##   batch_7  57  94  73  41  45  34  37 100  26   2  23   4
```

```r
# barplot(t(t(table(batch=tmp$batch, cluster=tmp$cluster))/colSums(table(batch=tmp$batch, cluster=tmp$cluster))*100), col=1:12)
```

Finally, we check the QC statistics to make sure that there isn't a cluster that behaves abnormally.


```r
qc <- read.table(paste0(dir, "data/QCstats_allCells.tsv"))
qc <- qc[names(cut[cut!=0]),]

par(mfrow=c(2,2))
boxplot(log10(qc$libSize)~cut[cut!=0], las=2, xlab="cluster", ylab=expression('log'[10]*' library size'), col=1:12)
boxplot(qc$nGenes/1e3~cut[cut!=0], las=2, xlab="cluster", ylab="number of genes expressed x 1000", col=1:12)
boxplot(qc$mit/qc$libSize*100~cut[cut!=0], las=2, xlab="cluster", ylab="% in MT", col=1:12)
boxplot(qc$ercc/qc$libSize*100~cut[cut!=0], las=2, xlab="cluster", ylab="% in ERCC", col=1:12)
```

![](03_clustering_files/figure-html/qc-1.png)<!-- -->

### Cluster markers {.tabset}

To get an initial idea of the identity of each cluster, we use the `findMarkers` function to identify genes with large fold-changes in a particular cluster compared to the rest. We require genes to be significant in all comparisons against all other clusters; this returns genes that are most specific to each cluster. But can be problematic if there are some closely related clusters that share markers.

This approach returns significant genes for all but one cluster.


```r
## use normalised counts and block by batch instead of using batch corrected counts
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.Rds"))

## add cluster info
stopifnot(identical(colnames(logcounts(sce)), names(cut)))
sce$clusters <- cut

## remove outlier cell
sce <- sce[,sce$clusters>0]

## find markers
keep <- rowMeans(logcounts(sce)) > 0.1
markersDE <- findMarkers(sce, groups = sce$clusters, block=sce$batch, direction="up", subset.row=keep, pval.type="all")
unlist(lapply(markersDE, function(x) sum(x$FDR<0.05)))
```

```
##   1   2   3   4   5   6   7   8   9  10  11  12 
## 450   0  82 107 324  11  29  48  62  24 296 803
```

```r
saveRDS(markersDE, file=paste0(dir, "results/markerGenes_pval_all.Rds"))
```

To recover markers for cluster 2, we instead require the gene to be significantly different against at least 8 of the 12 clusters. This of course returns higher numbers of significant genes.


```r
markersDE.some <- findMarkers(sce, groups = sce$clusters, block=sce$batch, direction="up", subset.row=keep, pval.type="some", min.prop=0.75)
unlist(lapply(markersDE.some, function(x) sum(x$FDR<0.05)))
```

```
##    1    2    3    4    5    6    7    8    9   10   11   12 
## 1348  116  722  427 1013  370  244  287  191  282  715 1322
```

```r
saveRDS(markersDE, file=paste0(dir, "results/markerGenes_pval_some0.75.Rds"))
```

Below is the expression of the top 10 genes found for each cluster:


```r
th <- theme_bw() + theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.y = element_text(size=12), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5))

plotGeneOnUMAP <- function(umap=umap, data=sce, clusters=clusters, gene=gene){
  df <- data.frame(x=umap$x, y=umap$y, expr=as.numeric(logcounts(sce[gene,])), cluster=clusters)
  df <- df[df$cluster>0,]
  df <- df[order(df$expr),]
  
  p <- list()
  p[[1]] <- ggplot(df, aes(x,y)) + geom_point(aes(colour=expr), alpha = 0.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(rowData(sce)[id,1]) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(legend.position = "none") + guides(colour = guide_colorbar(title.position = "bottom"))
  p[[2]] <- ggplot(df, aes(as.factor(cluster), expr)) + geom_boxplot(fill=1:12) + ggtitle(rowData(sce)[id,1]) + xlab("") + ylab("") + th
  return(p)
}
```

#### Cluster 1


```r
cluster <- 1
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster1-1.png)<!-- -->![](03_clustering_files/figure-html/cluster1-2.png)<!-- -->![](03_clustering_files/figure-html/cluster1-3.png)<!-- -->![](03_clustering_files/figure-html/cluster1-4.png)<!-- -->![](03_clustering_files/figure-html/cluster1-5.png)<!-- -->![](03_clustering_files/figure-html/cluster1-6.png)<!-- -->![](03_clustering_files/figure-html/cluster1-7.png)<!-- -->![](03_clustering_files/figure-html/cluster1-8.png)<!-- -->![](03_clustering_files/figure-html/cluster1-9.png)<!-- -->![](03_clustering_files/figure-html/cluster1-10.png)<!-- -->

#### Cluster 2


```r
cluster <- 2
for(i in 1:10){
  id <- row.names(markersDE.some[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster2-1.png)<!-- -->![](03_clustering_files/figure-html/cluster2-2.png)<!-- -->![](03_clustering_files/figure-html/cluster2-3.png)<!-- -->![](03_clustering_files/figure-html/cluster2-4.png)<!-- -->![](03_clustering_files/figure-html/cluster2-5.png)<!-- -->![](03_clustering_files/figure-html/cluster2-6.png)<!-- -->![](03_clustering_files/figure-html/cluster2-7.png)<!-- -->![](03_clustering_files/figure-html/cluster2-8.png)<!-- -->![](03_clustering_files/figure-html/cluster2-9.png)<!-- -->![](03_clustering_files/figure-html/cluster2-10.png)<!-- -->

#### Cluster 3


```r
cluster <- 3
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster3-1.png)<!-- -->![](03_clustering_files/figure-html/cluster3-2.png)<!-- -->![](03_clustering_files/figure-html/cluster3-3.png)<!-- -->![](03_clustering_files/figure-html/cluster3-4.png)<!-- -->![](03_clustering_files/figure-html/cluster3-5.png)<!-- -->![](03_clustering_files/figure-html/cluster3-6.png)<!-- -->![](03_clustering_files/figure-html/cluster3-7.png)<!-- -->![](03_clustering_files/figure-html/cluster3-8.png)<!-- -->![](03_clustering_files/figure-html/cluster3-9.png)<!-- -->![](03_clustering_files/figure-html/cluster3-10.png)<!-- -->

#### Cluster 4


```r
cluster <- 4
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster4-1.png)<!-- -->![](03_clustering_files/figure-html/cluster4-2.png)<!-- -->![](03_clustering_files/figure-html/cluster4-3.png)<!-- -->![](03_clustering_files/figure-html/cluster4-4.png)<!-- -->![](03_clustering_files/figure-html/cluster4-5.png)<!-- -->![](03_clustering_files/figure-html/cluster4-6.png)<!-- -->![](03_clustering_files/figure-html/cluster4-7.png)<!-- -->![](03_clustering_files/figure-html/cluster4-8.png)<!-- -->![](03_clustering_files/figure-html/cluster4-9.png)<!-- -->![](03_clustering_files/figure-html/cluster4-10.png)<!-- -->

#### Cluster 5


```r
cluster <- 5
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster5-1.png)<!-- -->![](03_clustering_files/figure-html/cluster5-2.png)<!-- -->![](03_clustering_files/figure-html/cluster5-3.png)<!-- -->![](03_clustering_files/figure-html/cluster5-4.png)<!-- -->![](03_clustering_files/figure-html/cluster5-5.png)<!-- -->![](03_clustering_files/figure-html/cluster5-6.png)<!-- -->![](03_clustering_files/figure-html/cluster5-7.png)<!-- -->![](03_clustering_files/figure-html/cluster5-8.png)<!-- -->![](03_clustering_files/figure-html/cluster5-9.png)<!-- -->![](03_clustering_files/figure-html/cluster5-10.png)<!-- -->

#### Cluster 6


```r
cluster <- 6
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster6-1.png)<!-- -->![](03_clustering_files/figure-html/cluster6-2.png)<!-- -->![](03_clustering_files/figure-html/cluster6-3.png)<!-- -->![](03_clustering_files/figure-html/cluster6-4.png)<!-- -->![](03_clustering_files/figure-html/cluster6-5.png)<!-- -->![](03_clustering_files/figure-html/cluster6-6.png)<!-- -->![](03_clustering_files/figure-html/cluster6-7.png)<!-- -->![](03_clustering_files/figure-html/cluster6-8.png)<!-- -->![](03_clustering_files/figure-html/cluster6-9.png)<!-- -->![](03_clustering_files/figure-html/cluster6-10.png)<!-- -->

#### Cluster 7


```r
cluster <- 7
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster7-1.png)<!-- -->![](03_clustering_files/figure-html/cluster7-2.png)<!-- -->![](03_clustering_files/figure-html/cluster7-3.png)<!-- -->![](03_clustering_files/figure-html/cluster7-4.png)<!-- -->![](03_clustering_files/figure-html/cluster7-5.png)<!-- -->![](03_clustering_files/figure-html/cluster7-6.png)<!-- -->![](03_clustering_files/figure-html/cluster7-7.png)<!-- -->![](03_clustering_files/figure-html/cluster7-8.png)<!-- -->![](03_clustering_files/figure-html/cluster7-9.png)<!-- -->![](03_clustering_files/figure-html/cluster7-10.png)<!-- -->

#### Cluster 8


```r
cluster <- 8
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster8-1.png)<!-- -->![](03_clustering_files/figure-html/cluster8-2.png)<!-- -->![](03_clustering_files/figure-html/cluster8-3.png)<!-- -->![](03_clustering_files/figure-html/cluster8-4.png)<!-- -->![](03_clustering_files/figure-html/cluster8-5.png)<!-- -->![](03_clustering_files/figure-html/cluster8-6.png)<!-- -->![](03_clustering_files/figure-html/cluster8-7.png)<!-- -->![](03_clustering_files/figure-html/cluster8-8.png)<!-- -->![](03_clustering_files/figure-html/cluster8-9.png)<!-- -->![](03_clustering_files/figure-html/cluster8-10.png)<!-- -->

#### Cluster 9


```r
cluster <- 9
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster9-1.png)<!-- -->![](03_clustering_files/figure-html/cluster9-2.png)<!-- -->![](03_clustering_files/figure-html/cluster9-3.png)<!-- -->![](03_clustering_files/figure-html/cluster9-4.png)<!-- -->![](03_clustering_files/figure-html/cluster9-5.png)<!-- -->![](03_clustering_files/figure-html/cluster9-6.png)<!-- -->![](03_clustering_files/figure-html/cluster9-7.png)<!-- -->![](03_clustering_files/figure-html/cluster9-8.png)<!-- -->![](03_clustering_files/figure-html/cluster9-9.png)<!-- -->![](03_clustering_files/figure-html/cluster9-10.png)<!-- -->

#### Cluster 10


```r
cluster <- 10
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster10-1.png)<!-- -->![](03_clustering_files/figure-html/cluster10-2.png)<!-- -->![](03_clustering_files/figure-html/cluster10-3.png)<!-- -->![](03_clustering_files/figure-html/cluster10-4.png)<!-- -->![](03_clustering_files/figure-html/cluster10-5.png)<!-- -->![](03_clustering_files/figure-html/cluster10-6.png)<!-- -->![](03_clustering_files/figure-html/cluster10-7.png)<!-- -->![](03_clustering_files/figure-html/cluster10-8.png)<!-- -->![](03_clustering_files/figure-html/cluster10-9.png)<!-- -->![](03_clustering_files/figure-html/cluster10-10.png)<!-- -->

#### Cluster 11


```r
cluster <- 11
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster11-1.png)<!-- -->![](03_clustering_files/figure-html/cluster11-2.png)<!-- -->![](03_clustering_files/figure-html/cluster11-3.png)<!-- -->![](03_clustering_files/figure-html/cluster11-4.png)<!-- -->![](03_clustering_files/figure-html/cluster11-5.png)<!-- -->![](03_clustering_files/figure-html/cluster11-6.png)<!-- -->![](03_clustering_files/figure-html/cluster11-7.png)<!-- -->![](03_clustering_files/figure-html/cluster11-8.png)<!-- -->![](03_clustering_files/figure-html/cluster11-9.png)<!-- -->![](03_clustering_files/figure-html/cluster11-10.png)<!-- -->

#### Cluster 12


```r
cluster <- 12
for(i in 1:10){
  id <- row.names(markersDE[[cluster]])[i]
  plots <- plotGeneOnUMAP(umap = umap[-which(cut==0),], data = sce, clusters = cut[cut>0], gene = id)
  print(ggarrange(plotlist = plots, ncol=2, nrow=1))
}
```

![](03_clustering_files/figure-html/cluster12-1.png)<!-- -->![](03_clustering_files/figure-html/cluster12-2.png)<!-- -->![](03_clustering_files/figure-html/cluster12-3.png)<!-- -->![](03_clustering_files/figure-html/cluster12-4.png)<!-- -->![](03_clustering_files/figure-html/cluster12-5.png)<!-- -->![](03_clustering_files/figure-html/cluster12-6.png)<!-- -->![](03_clustering_files/figure-html/cluster12-7.png)<!-- -->![](03_clustering_files/figure-html/cluster12-8.png)<!-- -->![](03_clustering_files/figure-html/cluster12-9.png)<!-- -->![](03_clustering_files/figure-html/cluster12-10.png)<!-- -->

###


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
##  [1] ggpubr_0.2.3                magrittr_1.5               
##  [3] ggplot2_3.2.1               RColorBrewer_1.1-2         
##  [5] dynamicTreeCut_1.63-1       scran_1.14.1               
##  [7] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
##  [9] DelayedArray_0.12.0         BiocParallel_1.20.0        
## [11] matrixStats_0.55.0          Biobase_2.46.0             
## [13] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [15] IRanges_2.20.0              S4Vectors_0.24.0           
## [17] BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2               rsvd_1.0.2              
##  [3] locfit_1.5-9.1           lattice_0.20-38         
##  [5] assertthat_0.2.1         digest_0.6.22           
##  [7] R6_2.4.0                 evaluate_0.14           
##  [9] pillar_1.4.2             zlibbioc_1.32.0         
## [11] rlang_0.4.1              lazyeval_0.2.2          
## [13] rstudioapi_0.10          irlba_2.3.3             
## [15] Matrix_1.2-17            rmarkdown_1.16          
## [17] labeling_0.3             BiocNeighbors_1.4.1     
## [19] statmod_1.4.32           stringr_1.4.0           
## [21] igraph_1.2.4.1           RCurl_1.95-4.12         
## [23] munsell_0.5.0            compiler_3.6.1          
## [25] vipor_0.4.5              BiocSingular_1.2.0      
## [27] xfun_0.10                pkgconfig_2.0.3         
## [29] ggbeeswarm_0.6.0         htmltools_0.4.0         
## [31] tidyselect_0.2.5         gridExtra_2.3           
## [33] tibble_2.1.3             GenomeInfoDbData_1.2.2  
## [35] edgeR_3.28.0             viridisLite_0.3.0       
## [37] withr_2.1.2              crayon_1.3.4            
## [39] dplyr_0.8.3              bitops_1.0-6            
## [41] grid_3.6.1               gtable_0.3.0            
## [43] scales_1.0.0             dqrng_0.2.1             
## [45] stringi_1.4.3            ggsignif_0.6.0          
## [47] XVector_0.26.0           viridis_0.5.1           
## [49] limma_3.42.0             scater_1.14.1           
## [51] DelayedMatrixStats_1.8.0 cowplot_1.0.0           
## [53] tools_3.6.1              glue_1.3.1              
## [55] beeswarm_0.2.3           purrr_0.3.3             
## [57] yaml_2.2.0               colorspace_1.4-1        
## [59] knitr_1.25
```
