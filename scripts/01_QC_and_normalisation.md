---
title: "Single-cell RNA-seq of mouse early heart development"
date: '05 November, 2019'
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



## QC and normalisation

We have generated single-cell RNA-seq data from *unbiased* sampling of the developing heart in mouse embryos. the sampling covers the earliest stages of heart development, when the cardiac crescent becomes evident, up to the linear heart tube (LHT) stage. Embryos were staged depending on their cardiac crescent morphology and classified as stages 0 to 3, or LHT, as defined in Tyser et al., eLife, 2016. An additional sample was collected just before the left and right portions of the prospective cardiac crescent fuse (stage -1). All data were collected across seven different batches.

Data are from SMART-seq2 protocol, sequenced in an Illumina HiSeq 2500, to generate 125bp paired-end fragments. The first batch is 40 cells, used as a pilot to ensure the protocols for library prep were working properly. Raw data from the pilot can be accessed at ArrayExpress [E-MTAB-7403](https://www.ebi.ac.uk/arrayexpress/...) and from batches 2 to 7 at the European Nucleotide Archive under project [PRJEB14363](https://www.ebi.ac.uk/ena/data/view/PRJEB14363).

Sequencing data were aligned to the mouse reference genome (`mm10` supplemented with ERCC spike-in sequences) with `GSNAP` and the fragments aligned to each gene annotated in `Ensembl v87` were quantified with `HTSeq`. The counts for each cell were compiled into a single count matrix that can be downloaded from: `add gitlab link`. 

Sample metadata is provided as `Supplementary Table 1` with the paper. Download from [here](add_link).

Gene information is available in the `data` folder.

----

First, we load the count data, sample metadata and gene information.


```r
## count matrix
data <- readRDS(paste0(dir, "data/heartData_unbiased.RAW.Rds"))

## sample metadata
meta <- read.table(paste0(dir, "data/SupplementaryTable1.tab"), header = TRUE, sep="\t", stringsAsFactors = FALSE)
stopifnot(identical(colnames(data), meta$cell))

## gene information (from Ensembl version 87)
ann <- read.table(paste0(dir, "data/Mus_musculus.GRCm38.87.tsv"), sep="\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(ann) <- c("gene", "chr", "start", "end", "strand")
```

Cells from different stages were collected across seven different batches; all batches but batch 1 (pilot) contain data from at least two different stages, and all stages but stage -1 are represented in at least two different batches.


```r
table(batch=meta$batch, stage=meta$stage)
```

```
##      stage
## batch  -1   0   1   2   3 LHT
##     1   0   0   0   0   0  40
##     2   0   0   0 379   0 379
##     3   0 247 132   0   0   0
##     4   0   0  87 303 350   0
##     5   0   0 379 241 136   0
##     6   0   0 357 398   0   0
##     7 374 378   0   0   0   0
```

### Quality control

To assess the quality of the data, we use library size, general mapping statistics such as proportion of reads in mitochondrial genes and spike-ins, and the total number of genes detected.


```r
## remove the mapping stats from count data
mapping.stats <- data[49764:49768,]
data <- data[1:49763,]

## compute general QC metrics
qc <- data.frame(libSize=colSums(data), 
                 nGenes=apply(data, 2, function(x) length(which(x>0))),
                 mit=colSums(data[row.names(ann[ann$chr=="MT",]),]),
                 ercc=colSums(data[grep("ERCC",row.names(data)),]))

## plot
plots <- list()
## total counts in genes + spike-ins
plots[[1]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=log10(libSize+1))) + geom_violin() + geom_boxplot(width=0.05) + theme_classic() + theme(legend.position="none") + ylab(expression('log'[10]*' library size')) + xlab("batch") + geom_hline(yintercept = log10(50000), lwd=0.25, lty=2, col="grey") + ggtitle("total reads in genes") + theme(plot.title = element_text(face="bold", hjust=0.5))

## genes detected
plots[[2]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=nGenes)) + geom_violin() + geom_boxplot(width=0.05) + theme_classic() + theme(legend.position="none") + ylab("total genes") + xlab("batch") + geom_hline(yintercept = 6000, lwd=0.25, lty=2, col="grey") + ggtitle("number of genes detected") + theme(plot.title = element_text(face="bold", hjust=0.5))

## mitochondrial %
plots[[3]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=mit/libSize*100)) + geom_violin() + geom_boxplot(width=0.05) + theme_classic() + theme(legend.position="none") + ylab("% reads in MT genes") + xlab("batch") + geom_hline(yintercept = 15, lwd=0.25, lty=2, col="grey") + ggtitle("% reads in mitochondrial genes") + theme(plot.title = element_text(face="bold", hjust=0.5))

## spike-ins %
plots[[4]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=ercc/libSize*100)) + geom_violin(scale="width") + geom_boxplot(width=0.05) + theme_classic() + theme(legend.position="none") + ylab("% reads in spike-ins") + xlab("batch") + geom_hline(yintercept = 30, lwd=0.25, lty=2, col="grey") + ggtitle("% reads in ERCC spike-ins") + theme(plot.title = element_text(face="bold", hjust=0.5))

ggarrange(plotlist = plots, ncol = 2, nrow = 2)
```

![](01_QC_and_normalisation_files/figure-html/qc-1.png)<!-- -->

The pilot data was sequenced at higher depth, given the small number of samples. Batch 2 -the first of the large-scale collections- is of lower quality than the rest, with smaller libraries and fewer genes detected. Also, the dilution of the ERCC spike-ins was incorrect and they were not detected in this batch.

To pass quality-control, cells need to:

- Have more than 50,000 reads mapped to annotated genes.
- Have more than 6,000 genes detected.
- Have less than 15% of their reads assigned to mitochondrial genes.
- Have less than 30% of reads mapped to ERCC spike-ins.


```r
badQual <- which(qc$libSize < 50e3 | qc$mit/(qc$libSize+1)*100 >= 15 | qc$ercc/(qc$libSize+1)*100 >= 30 | qc$nGenes <= 6000)
# length(badQual) # 1075 (25.72%)
```

With this criteria, 1075 (25.72%) cells fail and are removed from downstream analyses.


```r
## remove bad quality samples from count matris
stopifnot(identical(row.names(qc), colnames(data)))
data <- data[,-badQual]
data <- data[rowSums(data)>0,] ## remove non-expressed genes

## and from metadata
stopifnot(identical(row.names(qc), meta$cellID))
meta <- meta[-badQual,]

stopifnot(identical(colnames(data), meta$cellID))

## save QC stats for future reference
qc$pass <- ifelse(row.names(qc) %in% colnames(data), 1, 0)
write.table(qc, file=paste0(dir, "data/QCstats_allCells.tsv"), quote = FALSE, sep="\t")
```

The clean dataset now consist of 3105 cells that collectively express 37621 genes.

Batch 2 is the one affected the most, loosing more than half of its cells; but the other batches retain around 85% of the data.


```r
table(batch=meta$batch, stage=meta$stage)
```

```
##      stage
## batch  -1   0   1   2   3 LHT
##     1   0   0   0   0   0  40
##     2   0   0   0 159   0 133
##     3   0 173 117   0   0   0
##     4   0   0  71 262 320   0
##     5   0   0 328 227 124   0
##     6   0   0 313 302   0   0
##     7 227 309   0   0   0   0
```

### Normalisation

The next thing we need before analysing this clean, high-quality dataset, is to normalise the counts to account for differences in sequencing depth and other composition biases.

We use the method implemented in `scran` to estimate size factors for each cell. Spike-ins are treated separately and normalised by total counts in spike-ins.


```r
## set up SingleCellExperiment object
genes <- grep("ENSMUSG", row.names(data), value = TRUE) # 37529
spikes <- grep("ERCC", row.names(data), value = TRUE) # 92

# metadata row.names should be counts col.names
m <- meta
row.names(m) <- m$cellID
m$cellID <- NULL
stopifnot(identical(row.names(m), colnames(data)))

# add also gene info
ann <- ann[genes,]
stopifnot(identical(row.names(ann), row.names(data[genes,])))
# need info for spikes
# spikes <- data.frame(gene=spikes, chr="ERCC", start=1, end=1, strand=1, row.names = spikes)
# ann <- rbind(ann, spikes)
# stopifnot(identical(row.names(ann), row.names(data)))

## sce object
sce <- SingleCellExperiment(assays = list(counts=as.matrix(data[genes,])), colData = m, rowData = ann)

## specify spike ins
# is.spike <- grepl("^ERCC", rownames(sce))
# isSpike(sce, "ERCC") <- is.spike

## normalisation
# plot(density(log10(rowMeans(data[genes,]))), bty="l", main="")
# abline(v=log10(1), lty=2)
## the default filter of mean>=1 is appropriate

## pre-cluster the data to protect the size factor estimation from too many DEGs
set.seed(0)
clusters  <- quickCluster(sce, min.size = 100, method = "igraph")
saveRDS(clusters, file=paste0(dir, "data/quickCluster.Rds"))

## estimate size factors
sce  <- computeSumFactors(sce, cluster = clusters, min.mean = 1)
sf <- sizeFactors(sce)
names(sf) <- colnames(counts(sce))
write.table(sf, file=paste0(dir, "data/sizeFactors_unbiasedDataset_minMean1.tsv"), quote = FALSE, sep = "\t")

plot(sf, colSums(counts(sce))/1e6, pch=16, xlab="size factors", ylab="library size (millions)", bty="l")
abline(lm((colSums(counts(sce))/1e6)~sf))
```

![](01_QC_and_normalisation_files/figure-html/sizeFactors-1.png)<!-- -->

```r
## spike-ins have their own normalisation factors, following scran's method.
# sce <- computeSpikeFactors(sce, type="ERCC", general.use = FALSE)
```

Finally, we use the size factors to normalise endogenous gene counts.


```r
sce <- normalize(sce)
saveRDS(sce, file=paste0(dir, "data/sce_goodQual.NORM.Rds"))

## endogenouse genes
dataNorm <- logcounts(sce)
saveRDS(dataNorm, file = paste0(dir, "data/heartData_unbiased.goodQual.NORM.Rds"))
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
##  [1] RColorBrewer_1.1-2          ggpubr_0.2.3               
##  [3] magrittr_1.5                scater_1.12.2              
##  [5] ggplot2_3.2.1               scran_1.12.1               
##  [7] SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.1
##  [9] DelayedArray_0.10.0         BiocParallel_1.18.1        
## [11] matrixStats_0.55.0          Biobase_2.44.0             
## [13] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
## [15] IRanges_2.18.3              S4Vectors_0.22.1           
## [17] BiocGenerics_0.30.0        
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.2               rsvd_1.0.2              
##  [3] locfit_1.5-9.1           lattice_0.20-38         
##  [5] assertthat_0.2.1         digest_0.6.22           
##  [7] R6_2.4.0                 dynamicTreeCut_1.63-1   
##  [9] evaluate_0.14            pillar_1.4.2            
## [11] zlibbioc_1.30.0          rlang_0.4.1             
## [13] lazyeval_0.2.2           rstudioapi_0.10         
## [15] irlba_2.3.3              Matrix_1.2-17           
## [17] rmarkdown_1.16           labeling_0.3            
## [19] BiocNeighbors_1.2.0      statmod_1.4.32          
## [21] stringr_1.4.0            igraph_1.2.4.1          
## [23] RCurl_1.95-4.12          munsell_0.5.0           
## [25] vipor_0.4.5              compiler_3.6.1          
## [27] BiocSingular_1.0.0       xfun_0.10               
## [29] pkgconfig_2.0.3          ggbeeswarm_0.6.0        
## [31] htmltools_0.4.0          tidyselect_0.2.5        
## [33] gridExtra_2.3            tibble_2.1.3            
## [35] GenomeInfoDbData_1.2.1   edgeR_3.26.8            
## [37] viridisLite_0.3.0        withr_2.1.2             
## [39] crayon_1.3.4             dplyr_0.8.3             
## [41] bitops_1.0-6             grid_3.6.1              
## [43] gtable_0.3.0             scales_1.0.0            
## [45] dqrng_0.2.1              stringi_1.4.3           
## [47] ggsignif_0.6.0           XVector_0.24.0          
## [49] viridis_0.5.1            limma_3.40.6            
## [51] DelayedMatrixStats_1.6.1 cowplot_1.0.0           
## [53] tools_3.6.1              glue_1.3.1              
## [55] beeswarm_0.2.3           purrr_0.3.3             
## [57] yaml_2.2.0               colorspace_1.4-1        
## [59] knitr_1.25
```
