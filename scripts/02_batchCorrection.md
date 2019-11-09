---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '09 November, 2019'
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




### Batch correction

We begin with the good-quality, normalised data produced in the `01_QC_and_normalisation.Rmd` script.


```r
## sce object
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.Rds"))
```

Now we need to check for batch effects. For this, we begin by visualising the data in low dimensional space, with UMAP.


```r
set.seed(497)
sce <- runUMAP(sce)

umap <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap) <- c("x", "y")
write.table(umap, paste0(dir,"results/umapCoords.tab"), quote=FALSE, sep="\t")

## colour by batch
umap$batch <- sce$batch
## randomise order per batch
tmp <- umap[sample(row.names(umap), nrow(umap), replace = FALSE),]
plot(tmp$x, tmp$y, pch=16, cex=0.75, col=tmp$batch, xlab="tSNE-dim1", ylab="tSNE-dim2", bty="l")
legend("topright", legend = levels(tmp$batch), col=1:7, pch=16, cex=0.85)
```

![](02_batchCorrection_files/figure-html/umap-1.png)<!-- -->

There is strong segregation of the cells based on their batch of origin, especially for batch_7. Thus, we need to correct for this effect before proceeding.

----

To correct the batch effects we use the method proposed by Haghverdi et al. (*Nature Biotechnology*, 2018), based on mutual nearest neighbours (MNNs), but use the quicker implementation in `fastMNN`. 

We restrict the calculations to highly varaible genes, defined using the method proposed by Brennecke et al. (*Nature Meethods*, 2015). We retain the 2000 genes with largest variation, but remove any mitochondrial and sexually dimorphic genes.


```r
hvgs <- modelGeneCV2(sce, block=sce$batch)
# plot(hvgs$mean, hvgs$total, log="xy", col=ifelse(hvgs$FDR<0.05, 'red', 'black'))
# points(hvgs$mean, hvgs$trend, col="blue", pch=16, cex=0.5)

hvgs <- getTopHVGs(hvgs, var.field="ratio", n=2000)

## we remove mitochondrial and sexually dimorphic genes from the list of HVGs
remove <- rowData(sce)[rowData(sce)$chr %in% c("Y", "MT"),]
remove <- rbind(remove, rowData(sce)[rowData(sce)$gene=='Xist',])

hvgs <- setdiff(hvgs, row.names(remove))
write.table(hvgs, file=paste0(dir, "results/HVGs_minMean1_FDR0.05.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
```

For batch correction, the order of the batches is important; we should use as *reference* the most diverse batch, to ensure cells in other batches are able to find their MNNs. We use the number of different stages per batch as a proxy for diversity.


```r
## The order of the batches alters the final result. The batch with greatest diversity
## should be used as reference (first) to maximise the chances of finding MNNs. 
## Assuming that batches with more stages have more diversity...
# table(sce$batch, sce$stage)
#          -1   0   1   2   3 LHT
# batch_1   0   0   0   0   0  40
# batch_2   0   0   0 159   0 133
# batch_3   0 173 117   0   0   0
# batch_4   0   0  71 262 320   0
# batch_5   0   0 328 227 124   0
# batch_6   0   0 313 302   0   0
# batch_7 227 309   0   0   0   0

## specifiy which cells corrspond to which batch; levels indicate the order for merging
batch <- factor(sce$batch, levels=paste0("batch_", c(5,4,6,3,2,7,1)))
sce.corr <- fastMNN(sce, batch=batch, subset.row = hvgs, correct.all = TRUE)
saveRDS(sce.corr, paste0(dir, "data/sce.goodQual.NORM.batchCorrected.Rds"))

dataCorr <- assay(sce.corr, "reconstructed")
saveRDS(dataCorr, file = paste0(dir, "data/heartData_unbiased.goodQual.NORM.batchCorrected.Rds"))
```

Finally, we visualise the corrected data.


```r
## visualise corrected data, which is already cosine normalised
set.seed(837)
sce.corr <- runUMAP(sce.corr, subset_row = hvgs, exprs_values = "reconstructed")

umap <- as.data.frame(reducedDim(sce.corr, "UMAP"))
colnames(umap) <- c("x", "y")
write.table(umap, paste0(dir,"results/umapCoords_corrected.tab"), quote=FALSE, sep="\t")

## colour by batch
umap$batch <- sce$batch
## randomise order per batch
tmp <- umap[sample(row.names(umap), nrow(umap), replace = FALSE),]
plot(tmp$x, tmp$y, pch=16, cex=0.75, col=tmp$batch, xlab="UMAP - dim1", ylab="UMAP - dim2", bty="l")
legend("topright", legend = levels(tmp$batch), col=1:7, pch=16, cex=0.85)
```

![](02_batchCorrection_files/figure-html/umap_corrected-1.png)<!-- -->

Cells from all batches are now well mixed and the clusters should reflect different biological populations.


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
##  [1] RColorBrewer_1.1-2          batchelor_1.2.1            
##  [3] scater_1.14.1               ggplot2_3.2.1              
##  [5] scran_1.14.1                SingleCellExperiment_1.8.0 
##  [7] SummarizedExperiment_1.16.0 DelayedArray_0.12.0        
##  [9] BiocParallel_1.20.0         matrixStats_0.55.0         
## [11] Biobase_2.46.0              GenomicRanges_1.38.0       
## [13] GenomeInfoDb_1.22.0         IRanges_2.20.0             
## [15] S4Vectors_0.24.0            BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2               rsvd_1.0.2              
##  [3] locfit_1.5-9.1           lattice_0.20-38         
##  [5] FNN_1.1.3                assertthat_0.2.1        
##  [7] digest_0.6.22            RSpectra_0.15-0         
##  [9] R6_2.4.0                 evaluate_0.14           
## [11] pillar_1.4.2             zlibbioc_1.32.0         
## [13] rlang_0.4.1              lazyeval_0.2.2          
## [15] rstudioapi_0.10          irlba_2.3.3             
## [17] Matrix_1.2-17            rmarkdown_1.16          
## [19] BiocNeighbors_1.4.1      statmod_1.4.32          
## [21] stringr_1.4.0            uwot_0.1.4              
## [23] igraph_1.2.4.1           RCurl_1.95-4.12         
## [25] munsell_0.5.0            compiler_3.6.1          
## [27] vipor_0.4.5              BiocSingular_1.2.0      
## [29] xfun_0.10                pkgconfig_2.0.3         
## [31] ggbeeswarm_0.6.0         htmltools_0.4.0         
## [33] tidyselect_0.2.5         gridExtra_2.3           
## [35] tibble_2.1.3             GenomeInfoDbData_1.2.2  
## [37] edgeR_3.28.0             viridisLite_0.3.0       
## [39] withr_2.1.2              crayon_1.3.4            
## [41] dplyr_0.8.3              bitops_1.0-6            
## [43] grid_3.6.1               gtable_0.3.0            
## [45] magrittr_1.5             scales_1.0.0            
## [47] dqrng_0.2.1              RcppParallel_4.4.4      
## [49] stringi_1.4.3            XVector_0.26.0          
## [51] viridis_0.5.1            limma_3.42.0            
## [53] DelayedMatrixStats_1.8.0 tools_3.6.1             
## [55] glue_1.3.1               beeswarm_0.2.3          
## [57] purrr_0.3.3              yaml_2.2.0              
## [59] colorspace_1.4-1         knitr_1.25
```
