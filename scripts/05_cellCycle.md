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



### Cell cycle assignment

We have normalised and batch-corrected the data, and clustered the cell into 12 different populations, that have been annotated based on their expression of marker genes.


```r
## normalised data with cluster annotation
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))
```

We can also use the transcriptomic profiles to infer the phase of the cell cycle each cell is in. For this, we use the pairs classifier developed by Antonio Scialdone (Scialdone et al., *Methods*, 2015).

Cells are identified in all three phases of the cell cycle, with nearly half in G2-M phase.


```r
## read the trained data, provided with scran
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

## classify cells
set.seed(100)
cellCyclePhase <- cyclone(sce, mm.pairs)
saveRDS(cellCyclePhase, file=paste0(dir, "results/cellCyclePredictions.Rds"))

sce$cellCycle <- cellCyclePhase$phases
table(sce$cellCycle)
```

```
## 
##   G1  G2M    S 
##  634 1532  938
```

```r
#  G1  G2M    S 
# 633 1533  939 
tmp <- sce$cellCycle
names(tmp) <- colnames(sce)
write.table(tmp, file=paste0(dir, "results/cellCyclePhase.tsv"), quote = FALSE, col.names = FALSE)

plot(cellCyclePhase$score$G1, cellCyclePhase$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16, col=as.factor(sce$cellCycle), bty="l")
legend("topright", legend = levels(as.factor(sce$cellCycle)), col = 1:3, pch = 16)
```

![](05_cellCycle_files/figure-html/cyclone-1.png)<!-- -->

And when looking in each cluster, most are proliferative, with only a couple of clusters having a large proportion of cells in G1 phase.


```r
## proportion per cluster
props <- t(table(sce$clusterAnn, sce$cellCycle)[-1,]/rowSums(table(sce$clusterAnn, sce$cellCycle)[-1,]))
props <- props[,order(props[1,])]

barplot(props, las=2, col=1:3, xlim=c(0,15))
legend("topright", rev(row.names(props)), pch=15, col=3:1, cex=0.75)
```

![](05_cellCycle_files/figure-html/cellCycle-1.png)<!-- -->

Me3 corresponds to the most mature cardiomyocyte subpopulation, so it makes sense for it to contain a higher proportion of cells in G1 phase. 

#### Stage-dependent changes

Something else that is interesting to investigate, is whether the proportions of cells in each cell cycle phase change across development, for each cluster. For this, we use the sample's stage.



```r
perStage <- list()
for(cluster in unique(sce$clusterAnn)){
  tmp <- sce[,sce$clusterAnn==cluster]
  perStage[[cluster]] <- as.matrix(table(tmp$stage, tmp$cellCycle))
  perStage[[cluster]] <- perStage[[cluster]]/rowSums(perStage[[cluster]])
}

cols <- unique(sce$clusterCol)
names(cols) <- unique(sce$clusterAnn)

barplot(c(perStage[["Me5"]][,1]*100, perStage[["Me7"]][,1]*100, perStage[["Me4"]][,1]*100, c(0,perStage[["Me6"]][,1]*100), perStage[["Me3"]][,1]*100, perStage[["Me8"]][,1]*100), las=2, width = 1, space = 0, col=rep(cols[paste0("Me",c(5,7,4,6,3,8))], each=6))
mtext(side=1, line=3, text = paste0("Me",c(5,7,4,6,3,8)), at=c(3.5,9.5,15.5,21.5,27.5,33.5), col=cols[paste0("Me",c(5,7,4,6,3,8))], font=2)
```

![](05_cellCycle_files/figure-html/perStage-1.png)<!-- -->

While Me5 and Me7 -which are the most undifferentiated progenitor cells- have constantly low proportions of cells in G1, for the differentiating clusters Me4, Me6 and Me3 the fraction of cells in G1 increases at later stages.




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
##  [1] RColorBrewer_1.1-2          scran_1.14.1               
##  [3] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
##  [5] DelayedArray_0.12.0         BiocParallel_1.20.0        
##  [7] matrixStats_0.55.0          Biobase_2.46.0             
##  [9] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [11] IRanges_2.20.0              S4Vectors_0.24.0           
## [13] BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2               rsvd_1.0.2              
##  [3] locfit_1.5-9.1           lattice_0.20-38         
##  [5] assertthat_0.2.1         digest_0.6.22           
##  [7] R6_2.4.0                 evaluate_0.14           
##  [9] ggplot2_3.2.1            pillar_1.4.2            
## [11] zlibbioc_1.32.0          rlang_0.4.1             
## [13] lazyeval_0.2.2           rstudioapi_0.10         
## [15] irlba_2.3.3              Matrix_1.2-17           
## [17] rmarkdown_1.16           BiocNeighbors_1.4.1     
## [19] statmod_1.4.32           stringr_1.4.0           
## [21] igraph_1.2.4.1           RCurl_1.95-4.12         
## [23] munsell_0.5.0            compiler_3.6.1          
## [25] vipor_0.4.5              BiocSingular_1.2.0      
## [27] xfun_0.10                pkgconfig_2.0.3         
## [29] ggbeeswarm_0.6.0         htmltools_0.4.0         
## [31] tidyselect_0.2.5         gridExtra_2.3           
## [33] tibble_2.1.3             GenomeInfoDbData_1.2.2  
## [35] edgeR_3.28.0             viridisLite_0.3.0       
## [37] crayon_1.3.4             dplyr_0.8.3             
## [39] bitops_1.0-6             grid_3.6.1              
## [41] gtable_0.3.0             magrittr_1.5            
## [43] scales_1.0.0             dqrng_0.2.1             
## [45] stringi_1.4.3            XVector_0.26.0          
## [47] viridis_0.5.1            limma_3.42.0            
## [49] scater_1.14.1            DelayedMatrixStats_1.8.0
## [51] tools_3.6.1              glue_1.3.1              
## [53] beeswarm_0.2.3           purrr_0.3.3             
## [55] yaml_2.2.0               colorspace_1.4-1        
## [57] knitr_1.25
```

