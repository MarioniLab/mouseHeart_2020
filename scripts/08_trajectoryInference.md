---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '13 June, 2020'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: true
    toc_depth: 4
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

o <- sample(1:ncol(sce), ncol(sce), replace = FALSE)
plot(reducedDim(sce)$x[o], reducedDim(sce)$y[o], pch=16, col=sce$clusterCol[o], xlab="", ylab="", axes=FALSE)
box(bty="l")
legend("bottomright", legend = names(cols), col=cols, cex=0.5, pch=16)
text(4, 0, labels = "ectoderm", col="firebrick", cex=0.75, font=2)
text(11.5, 3, labels = "endoderm", col="steelblue4", cex=0.75, font=2)
text(-7, 4, labels = "cardiac\nmesoderm", col="darkolivegreen", cex=0.75, font=2)
text(-4, 7.5, labels = "endothelium", col="wheat3", cex=0.75, font=2)
text(-5, 11, labels = "blood", col="burlywood3", cex=0.75, font=2)
```

![](08_trajectoryInference_files/figure-html/clusters-1.png)<!-- -->

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

![](08_trajectoryInference_files/figure-html/plot_markers-1.png)<!-- -->

We have also observed that the clusters of progenitor cells have transcriptional signatures consistent with different types of mesoderm (see `07.2_classificationMesodermCells.Rmd`). Particularly, cells from the `Me8` cluster mostly resemble cranial and presomitic mesoderm, instead of cardiac mesoderm.


```r
labels <- read.table(paste0(dir, "results/classesUnbiasedMe3-8.randForest.tsv"))
probs <- read.table(paste0(dir, "results/classesUnbiasedMe3-8.randForest.probs.tsv"))

calls <- data.frame(class=labels, max.prob=apply(probs, 2, max), closest=apply(probs, 2, function(x) max(x[x!=max(x)]) ))
calls$diff <- calls$max.prob-calls$closest
calls$pass <- ifelse(calls$diff > 0.15, 1, 0)

umap <- reducedDim(sce)
umap$class <- labels[match(row.names(umap), labels$V1),2]
umap$pass <- ifelse(row.names(umap) %in% row.names(calls[calls$pass==1,]), 1, 0)

order <- sample(row.names(umap[umap$pass==1,]))
plot(umap$x, umap$y, pch=16, col="lightgrey", axes=FALSE, xlab="", ylab="", xlim=c(-12,0), ylim=c(-7,4))
points(umap[umap$pass==0 & !is.na(umap$class),]$x, umap[umap$pass==0 & !is.na(umap$class),]$y, pch=16, col="grey")
points(umap[umap$pass==1,][order,]$x, umap[umap$pass==1,][order,]$y, pch=16, col=umap[umap$pass==1,][order,]$class)
box(bty="l")
legend("topleft", legend = levels(as.factor(umap$class)), pch=16, col=1:6)
```

![](08_trajectoryInference_files/figure-html/mesodermTypes-1.png)<!-- -->

We can use diffusion maps to study the differentiation trajectories from the progenitor populations towards cardiomyocytes. We exclude the `Me8` cells, given their distinct, non-cardiac, phenotype, and the cells from `Me7` that are assigned a cranial mesoderm signature.


```r
## select only cells in cardiac clusters, minus Me8
sce.cardiac <- sce[,sce$clusterAnn %in% paste0("Me", c(3:7))]
## remove Me7 cranial cells
sce.cardiac <- sce.cardiac[,-which(colnames(sce.cardiac) %in% row.names(calls[calls$class.V2 == "cranialMesoderm",]))]

sce.corr.cardiac <- sce.corr[,colnames(sce.cardiac)]
# plot(reducedDim(sce.cardiac)$x, reducedDim(sce.cardiac)$y, pch=16, col=sce.cardiac$clusterCol)
```

Using all cells in the cardiac mesoderm clusters, the diffusion map embedding arranges cells in a triangular structure, with the two progenitor populations (`Me5` and `Me7`) and the cardiomyocytes (`Me3`) at the tips; the yellow and orange clusters are intermediates, with the orange stretching out to both progenitor types.


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
plots[[2]] <- plot(dpt, col_by="branch", root=1, path=c(2,3), col_path=c("red", "blue")) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](08_trajectoryInference_files/figure-html/DM-1.png)<!-- -->

And using the diffusion pseudotime we can define branches. The turquoise Me7 cluster forms one branch, and so does the green Me5 cluster; the yellow, orange and blue clusters are all merged into the third branch.

This suggests that the two progenitor types have distinct differentiation trajectories to cardiomyocytes (instead of one progenitor type giving rise to the second).

#### Me5 to cardiomyocyte trajectory

To investigate in more detail the differentiation dynamics of the Me5 progenitors to mature cardiomyocytes, we now compute the diffusion map excluding the turquoise cells, which comprise a distinct branch.

Once again, cells are arranged in a triangular structure, with the Me5 progenitors and the mature cardiomyocytes at opposite ends; the yellow cluster segregates towards the other tip.


```r
## exclude Me7
sce.Me5 <- sce.cardiac[,-which(colnames(sce.cardiac) %in% colnames(sce.cardiac)[colData(sce.cardiac)$clusterAnn == "Me7"])]
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
plots[[2]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,3), col_path=c("red", "blue")) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](08_trajectoryInference_files/figure-html/me5-1.png)<!-- -->

And defining branches based on diffusion pseudotime separates the cardiomyocytes (Me3) from the yellow cluster (Me6) and the orange (Me4) + green (Me5) progenitors. 


```r
table(dpt.me5@branch[,1], sce.Me5[,row.names(dm.me5)]$clusterAnn)
```

```
##    
##     Me3 Me4 Me5 Me6
##   1  39 200 352   2
##   2 608   0   0   4
##   3   8   0   0  44
```


#### Me7 to cardiomyocyte trajectory

Now do the same to define the alternative trajectory from the turquoise progenitors towards cardiomyocytes.

Again, cells are arranged in a triangular structure, with the Me7 progenitors and the cardiomyocytes on opposite ends; and in this case, the cells from the orange cluster are in the other tip. 


```r
## exclude Me5
sce.Me7 <- sce.cardiac[,-which(colnames(sce.cardiac) %in% colnames(sce.cardiac)[colData(sce.cardiac)$clusterAnn == "Me5"])]
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
plots[[2]] <- plot(dpt.me7, col_by="branch", root=2, path=c(1,3), col_path=c("red", "blue")) + th
ggarrange(plotlist = plots, ncol=2, nrow = 1)
```

![](08_trajectoryInference_files/figure-html/me7-1.png)<!-- -->

Consistently, there is a branching point that separates the most distinct orange cells (Me4); but some orange cells are part of the branch containing the turquoise (Me7) and yellow (Me6) progenitors.


```r
table(dpt.me7@branch[,1], sce.Me7[,row.names(dm.me7)]$clusterAnn)
```

```
##    
##     Me3 Me4 Me6 Me7
##   1  31  46  59 421
##   2 670  49   1   0
##   3   0  97   0   0
```


#### Recap

All together, we can identify two different trajectories leading to mature cardiomyocytes:

- A path starting in the green cells (Me5), that goes through the orange cluster (Me4).
    - The cells from the yellow cluster (Me6) form a different branch from this path.
- Another path starting in the turquoise cells from Me7, that goes through the yellow (Me6) and some cells from the orange (Me4) clusters.

#### Dynamic gene expression

Based on the results above, we can define the pseudotime progression of cells in each trajectory, considering only the branch going to mature cardiomyocytes, and ignoring the cells branching off from the yellow and orange clusters in the Me5 and Me7 trajectories, respectively.

Many of the mature cardiomyocytes are shared between the two trajectories, but the intermediate cells linking them to the progenitors are different.


```r
## define pseudotime, starting in mature cardiomyocytes and progressing towards both progenitors
# which(dpt@tips[,1])
# plot(dpt, col_by="DPT684")
diffPseudotime <- dpt$DPT684
names(diffPseudotime) <- colnames(sce.cardiac)
## since it makes more sense to think of progenitors having lower pseudotime values than mature cells, we flip the values, to reflect this
diffPseudotime <- -diffPseudotime

## define pseudotime for the Me5->Me3 trajectory
## consider only cells in branches 1 and 2
cells.me5 <- row.names(dpt.me5@branch)[which(dpt.me5@branch[,1] %in% 1:2)]
diffPseudotime_me5 <- diffPseudotime[cells.me5]
names(diffPseudotime_me5) <- cells.me5

## define pseudotime for the Me7->Me3 trajectory
## consider only cells in branches 1 and 2
cells.me7 <- row.names(dpt.me7@branch)[which(dpt.me7@branch[,1] %in% 1:2)]
diffPseudotime_me7 <- diffPseudotime[cells.me7]
names(diffPseudotime_me7) <- cells.me7

plots <- list()
plots[[1]] <- ggplot(dm, aes(DC1, DC2, colour=sce.cardiac$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + ggtitle("Clusters") + th + theme(plot.title = element_text(face="bold", hjust = 0.5), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
plots[[2]] <- ggplot(dm, aes(DC1, DC2)) + geom_point(data=dm[setdiff(row.names(dm), names(diffPseudotime_me5)),], col="grey", size=1) + geom_point(data=dm[names(diffPseudotime_me5),], aes(col = diffPseudotime_me5)) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Greens")[-c(1,9)]) + labs(colour="DPT") + ggtitle("Me5 -> Me3 trajectory") + th + theme(plot.title = element_text(face="bold", hjust = 0.5), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
plots[[3]] <- ggplot(dm, aes(DC1, DC2)) + geom_point(data=dm[setdiff(row.names(dm), names(diffPseudotime_me7)),], col="grey", size=1) + geom_point(data=dm[names(diffPseudotime_me7),], aes(col = diffPseudotime_me7)) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Blues")[-c(1,9)]) + labs(colour="DPT") + ggtitle("Me7 -> Me3 trajectory") + th + theme(plot.title = element_text(face="bold", hjust = 0.5), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom") 

ggarrange(plotlist = plots, ncol = 3, nrow = 1, align="hv")
```

![](08_trajectoryInference_files/figure-html/pseudotime-1.png)<!-- -->

Based on these cell orderings, we can identify genes that are expressed dynamically along the trajectories. For this, we follow the method used in *Ibarra-Soria, Wajid et al., Nat Cell Biol, 2018*.

##### Me5 -> Me3 trajectory

First, we look at the Me5 progenitors differentiating into cardiomyocytes.


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
# length(deltaAIC_me5[deltaAIC_me5 < -100])  # 1968
# length(deltaAIC_me5[deltaAIC_me5 < -200])  # 994
# length(deltaAIC_me5[deltaAIC_me5 < -300])  # 630

dynGenes_me5 <- deltaAIC_me5[deltaAIC_me5 < -200]

# for(gene in names(deltaAIC_me5[985:995])){
#   df <- data.frame(diffPseudotime_me5, data[gene,])
#   fit0 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
# 
#   plot(df[,1], df[,2], pch=16, col=rgb(0,0,0,0.3), bty="l", xlab="Me5 ---> Me3", ylab=expression('log'[2]*' expression'), main=rowData(sce)[gene,1])
#   lines(diffPseudotime_me5, predict(fit0, diffPseudotime_me5), lwd=3, col="steelblue2")
#   lines(diffPseudotime_me5, predict(fit2, diffPseudotime_me5), lwd=3, col="indianred2")
# }
```

We identify 994 genes dynamically expressed along this trajectory. These can be clustered into different profiles.


```r
curves_me5 <- t(sapply(names(dynGenes_me5), function(x){
  df <- data.frame(diffPseudotime_me5, data[x,])
  fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
  return(predict(fit2, df[,1]))
}))
curvesStd_me5 <- t(apply(curves_me5, 1, function(x) x/max(x)))  # standarise

# cluster the different profiles
test <- cor(t(curvesStd_me5), method="spearman")
test.dissim <- sqrt(0.5*((1-test)))
test.dist <- as.dist(test.dissim)
test.clust <- hclust(test.dist, method="average")

clust.genes_me5 <- cutreeDynamic(test.clust,distM=as.matrix(test.dist), minClusterSize=30, method="hybrid", deepSplit = 1, verbose = 0)
names(clust.genes_me5) <- row.names(curvesStd_me5)
# table(clust.genes_me5)

plots <- list()
for(c in 1:max(clust.genes_me5)){
  g <- names(clust.genes_me5[clust.genes_me5==c])
  test <- curvesStd_me5[g,]
  mean <- colMeans(test)
  std <- apply(test, 2, sd)
  df.sub <- data.frame(x=diffPseudotime_me5, ymin=mean-std, ymax=mean+std)
  df.sub.avg <- data.frame(x=diffPseudotime_me5, y=mean)
  plots[[c]] <- ggplot() +
            geom_line(data=df.sub.avg, aes(x=x,y=y), colour="black" , size=1) +
            geom_ribbon(data=df.sub, aes(x=x,ymin=ymin,ymax=ymax), alpha=0.2, fill="black") +
            xlab(expression("pseudotime")) +
            ylab(expression("standaradised expr")) +
            ggtitle(paste0("Cluster ", c, " (n=", length(g), ")")) + ylim(0,1.2) + th
}
ggarrange(plotlist = plots, ncol = 4, nrow = 2)
```

![](08_trajectoryInference_files/figure-html/dynGenes_me5_clusters-1.png)<!-- -->

And these genes should be enriched for processes involved in heart development. And this is the case:


```r
### test GO enrichments
universe <- row.names(data)
all <- as.factor(as.numeric(universe %in% names(dynGenes_me5)))
names(all) <- universe

go.me5 <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test.me5 <- runTest(go.me5, algorithm = "classic", statistic = "Fisher" )
go.res.me5 <- GenTable(go.me5, Fisher.classic = go.test.me5, topNodes = length(score(go.test.me5)))
go.res.me5$Fisher.classic.adj <- p.adjust(go.res.me5$Fisher.classic, "fdr")
go.res.me5[is.na(go.res.me5$Fisher.classic.adj),]$Fisher.classic.adj <- go.res.me5[is.na(go.res.me5$Fisher.classic.adj),]$Fisher.classic
go.res.me5[c(1,4,5,13,18,20,24,25,26,36,39,40,45,55,57,64,70,81,82,86,105),-6]
```

```
##          GO.ID                              Term Annotated Significant Expected
## 1   GO:0061061      muscle structure development       529         142    41.32
## 4   GO:0072359    circulatory system development       856         170    66.86
## 5   GO:0007507                 heart development       511         124    39.91
## 13  GO:0048738 cardiac muscle tissue development       219          71    17.11
## 18  GO:0030154              cell differentiation      2637         337   205.96
## 20  GO:0030029      actin filament-based process       594         120    46.39
## 24  GO:0035051        cardiocyte differentiation       147          53    11.48
## 25  GO:0007010         cytoskeleton organization      1044         171    81.54
## 26  GO:0003007               heart morphogenesis       226          66    17.65
## 36  GO:0003206     cardiac chamber morphogenesis       125          45     9.76
## 39  GO:0030239                myofibril assembly        60          30     4.69
## 40  GO:0006936                muscle contraction       193          55    15.07
## 45  GO:0032879        regulation of localization      1883         244   147.07
## 55  GO:0007154                cell communication      3370         375   263.21
## 57  GO:0016477                    cell migration       967         148    75.53
## 64  GO:0140014          mitotic nuclear division       226          57    17.65
## 70  GO:0035295                  tube development       799         127    62.41
## 81  GO:0051301                     cell division       503          91    39.29
## 82  GO:0042127  regulation of cell proliferation      1078         153    84.20
## 86  GO:0060429            epithelium development       763         119    59.59
## 105 GO:0007049                        cell cycle      1339         173   104.58
##       Fisher.classic.adj
## 1                < 1e-30
## 4                < 1e-30
## 5                < 1e-30
## 13  1.84263333333333e-23
## 18  7.17909090909091e-22
## 20  1.74862142857143e-20
## 24  6.96794117647059e-20
## 25           1.18455e-19
## 26  2.61847894736842e-19
## 36  7.35237931034483e-17
## 39       1.159871875e-15
## 40  1.48649411764706e-15
## 45  2.49378947368421e-15
## 55           1.97425e-14
## 57           2.68498e-14
## 64  9.83661403508772e-14
## 70  2.25628571428571e-13
## 81  1.38731081081081e-12
## 82  4.26022368421053e-12
## 86         6.4163125e-12
## 105 2.73977551020408e-10
```

##### Me7 -> Me3 trajectory

Now we repeat the same for the trajectory from the Me7 progenitors.


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
# length(deltaAIC_me7[deltaAIC_me7 < -100])  # 2035
# length(deltaAIC_me7[deltaAIC_me7 < -200])  # 969
# length(deltaAIC_me7[deltaAIC_me7 < -300])  # 628

dynGenes_me7 <- deltaAIC_me7[deltaAIC_me7 < -200]

# i=860
# gene  <- names(deltaAIC_me7[i])
# df <- data.frame(diffPseudotime_me7, data[gene,]); colnames(df) <- c("DPT", "expr")
# fit0 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=0), data=df)
# fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
# ggplot(df, aes(DPT, expr)) + geom_point(colour=sce[,names(diffPseudotime_me7)]$clusterCol, alpha=0.5) + geom_line(aes(diffPseudotime_me7, predict(fit0, diffPseudotime_me7)), lwd=1, col="steelblue2") + geom_line(aes(diffPseudotime_me7, predict(fit2, diffPseudotime_me7)), lwd=1, col="indianred2") + th
```

969 genes are dynamically expressed along this pathway, and can be clustered into several profiles. 


```r
curves_me7 <- t(sapply(names(dynGenes_me7), function(x){
  df <- data.frame(diffPseudotime_me7, data[x,])
  fit2 <- locfit(df[,2]~lp(df[,1], nn=smooth, deg=2), data=df)
  return(predict(fit2, df[,1]))
}))
curvesStd_me7 <- t(apply(curves_me7, 1, function(x) x/max(x)))  # standarise

# cluster the different profiles
test <- cor(t(curvesStd_me7), method="spearman")
test.dissim <- sqrt(0.5*((1-test)))
test.dist <- as.dist(test.dissim)
test.clust <- hclust(test.dist, method="average")

clust.genes_me7 <- cutreeDynamic(test.clust,distM=as.matrix(test.dist), minClusterSize=15, method="hybrid", deepSplit = 1, verbose = 0)  
names(clust.genes_me7) <- row.names(curvesStd_me7)
# table(clust.genes_me7)

plots <- list()
for(c in 1:max(clust.genes_me7)){
  g <- names(clust.genes_me7[clust.genes_me7==c])
  test <- curvesStd_me7[g,]
  mean <- colMeans(test)
  std <- apply(test, 2, sd)
  df.sub <- data.frame(x=diffPseudotime_me7, ymin=mean-std, ymax=mean+std)
  df.sub.avg <- data.frame(x=diffPseudotime_me7, y=mean)
  plots[[c]] <- ggplot() +
            geom_line(data=df.sub.avg, aes(x=x,y=y), colour="black" , size=1) +
            geom_ribbon(data=df.sub, aes(x=x,ymin=ymin,ymax=ymax), alpha=0.2, fill="black") +
            xlab(expression("pseudotime")) +
            ylab(expression("standaradised expr")) +
            ggtitle(paste0("Cluster ", c, " (n=", length(g), ")")) + ylim(0,1.2) + th
}
ggarrange(plotlist = plots, ncol = 3, nrow = 3)
```

![](08_trajectoryInference_files/figure-html/dynGenes_me7_clusters-1.png)<!-- -->

This set is also enriched in heart development related terms.


```r
### test GO enrichments
universe <- row.names(data)
all <- as.factor(as.numeric(universe %in% names(dynGenes_me7)))
names(all) <- universe

go.me7 <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test.me7 <- runTest(go.me7, algorithm = "classic", statistic = "Fisher" )
go.res.me7 <- GenTable(go.me7, Fisher.classic = go.test.me7, topNodes = length(score(go.test.me7)))
go.res.me7$Fisher.classic.adj <- p.adjust(go.res.me7$Fisher.classic, "fdr")
go.res.me7[is.na(go.res.me7$Fisher.classic.adj),]$Fisher.classic.adj <- go.res.me7[is.na(go.res.me7$Fisher.classic.adj),]$Fisher.classic
go.res.me7[c(1,3,4,6,7,11,13,16,23,26,32,49,99),-6]
```

```
##         GO.ID                                 Term Annotated Significant
## 1  GO:0061061         muscle structure development       529         131
## 3  GO:0014706   striated muscle tissue development       346          92
## 4  GO:0007507                    heart development       517         116
## 6  GO:0030029         actin filament-based process       598         124
## 7  GO:0048738    cardiac muscle tissue development       223          69
## 11 GO:0072359       circulatory system development       867         152
## 13 GO:0006936                   muscle contraction       197          63
## 16 GO:0051146 striated muscle cell differentiation       238          69
## 23 GO:0007010            cytoskeleton organization      1046         164
## 26 GO:0035051           cardiocyte differentiation       150          49
## 32 GO:0060047                    heart contraction       159          50
## 49 GO:0051301                        cell division       505          93
## 99 GO:0016477                       cell migration       980         128
##    Expected   Fisher.classic.adj
## 1     39.53              < 1e-30
## 3     25.85          9.90875e-25
## 4     38.63          9.90875e-25
## 6     44.68        1.4070425e-23
## 7     16.66 1.71751666666667e-22
## 11    64.78            7.927e-22
## 13    14.72 1.65746363636364e-21
## 16    17.78 4.86944285714286e-21
## 23    78.15 8.68195238095238e-19
## 26    11.21           3.9635e-17
## 32    11.88           8.7197e-17
## 49    37.73 2.02391489361702e-14
## 99    73.22 8.98938144329897e-09
```



#### Relationship with developmental stages

The diffusion map orders cells based on the similarity of their transcriptomes. And given the paths that have been revealed, we infer that the algorithm is capturing the progression from the most naïve progenitors to the most mature cardiomyocytes. 

But is there any relationship between the differentiation dynamics and the developmental stage of the embryos?

The cells from embryos of different stages are arranged quite uniformly along the diffusion space, with the exception of the cells from stage -1, which preferentially occupy the earliest stages of the turquoise (Me7 -> Me6(Me4) -> Me3) trajectory and are pretty much absent from the green (Me5 -> Me4 -> Me3) path.


```r
dm$stage <- sce[,row.names(dm)]$stage
ggplot(dm, aes(DC1, DC2, colour=stage)) + geom_point(size=1, alpha=0.5) + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + facet_wrap(~stage) + th + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
```

![](08_trajectoryInference_files/figure-html/dm_perStage-1.png)<!-- -->

And by uniformly we mean consistent with the density of cells observed in the different areas of the diffusion space (more cells at the tips, especially in the turquoise and blue clusters). Thus, all stages but -1 follow similar distributions of density across the first two diffusion components.


```r
plots <- list()
plots[[1]] <- ggplot(dm, aes(DC1, colour=stage)) + geom_density() + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
plots[[2]] <- ggplot(dm, aes(DC2, colour=stage)) + geom_density() + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggarrange(plotlist = plots, ncol = 2)
```

![](08_trajectoryInference_files/figure-html/density_perStage-1.png)<!-- -->

This suggests that from the formation of the cardiac crescent at stage 0 onward, the cells do not show substantial transcriptional changes based on embryo age; instead, they are part of a continuous differentiation process that stays constant across these developmental stages.

#### Visualisation

Finally, for ease of interpretation, we plot the diffusion map so that the progenitor populations are at the top and the cardiomyocytes at the bottom (by plotting -1*DC1 on the y-axis and DC2 on the x-axis), which is a more intuitive interpretation.


```r
ggplot(dm, aes(DC2, -DC1, colour=sce.cardiac$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
```

![](08_trajectoryInference_files/figure-html/vis-1.png)<!-- -->



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
write.table(clust.genes_me5, file=paste0(dir, "results/dynGenes_Me5.tsv"), quote = FALSE, sep="\t", col.names = FALSE)
write.table(deltaAIC_me7, file = paste0(dir, "results/deltaAIC_Me7.tsv"), quote = FALSE, sep="\t", col.names = FALSE)
write.table(clust.genes_me7, file=paste0(dir, "results/dynGenes_Me7.tsv"), quote = FALSE, sep="\t", col.names = FALSE)

## GO enrichment
write.table(go.res.me5, file=paste0(dir, "results/dynamicGenes_Me5_GO.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
write.table(go.res.me7, file=paste0(dir, "results/dynamicGenes_Me7_GO.tsv"), quote = FALSE, sep="\t", row.names = FALSE)
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
##  [1] org.Mm.eg.db_3.10.0         topGO_2.37.0               
##  [3] SparseM_1.77                GO.db_3.10.0               
##  [5] AnnotationDbi_1.48.0        graph_1.64.0               
##  [7] dynamicTreeCut_1.63-1       locfit_1.5-9.1             
##  [9] gplots_3.0.1.1              RColorBrewer_1.1-2         
## [11] ggpubr_0.2.5                magrittr_1.5               
## [13] ggplot2_3.3.0               destiny_3.0.0              
## [15] scran_1.14.5                SingleCellExperiment_1.8.0 
## [17] SummarizedExperiment_1.16.0 DelayedArray_0.12.0        
## [19] BiocParallel_1.20.0         matrixStats_0.55.0         
## [21] Biobase_2.46.0              GenomicRanges_1.38.0       
## [23] GenomeInfoDb_1.22.0         IRanges_2.20.1             
## [25] S4Vectors_0.24.1            BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1             backports_1.1.5          RcppEigen_0.3.3.7.0     
##   [4] igraph_1.2.4.2           sp_1.3-2                 RcppHNSW_0.2.0          
##   [7] scater_1.14.4            digest_0.6.23            htmltools_0.4.0         
##  [10] viridis_0.5.1            gdata_2.18.0             memoise_1.1.0           
##  [13] openxlsx_4.1.3           limma_3.42.0             xts_0.11-2              
##  [16] colorspace_1.4-1         blob_1.2.0               haven_2.2.0             
##  [19] xfun_0.11                dplyr_0.8.3              crayon_1.3.4            
##  [22] RCurl_1.95-4.12          hexbin_1.28.0            zeallot_0.1.0           
##  [25] zoo_1.8-6                glue_1.3.1               gtable_0.3.0            
##  [28] zlibbioc_1.32.0          XVector_0.26.0           car_3.0-5               
##  [31] BiocSingular_1.2.0       DEoptimR_1.0-8           abind_1.4-5             
##  [34] VIM_4.8.0                scales_1.1.0             ggplot.multistats_1.0.0 
##  [37] DBI_1.0.0                edgeR_3.28.0             ggthemes_4.2.0          
##  [40] Rcpp_1.0.3               viridisLite_0.3.0        laeken_0.5.0            
##  [43] dqrng_0.2.1              foreign_0.8-72           bit_1.1-14              
##  [46] rsvd_1.0.2               proxy_0.4-23             vcd_1.4-4               
##  [49] farver_2.0.1             pkgconfig_2.0.3          nnet_7.3-12             
##  [52] labeling_0.3             tidyselect_0.2.5         rlang_0.4.2             
##  [55] munsell_0.5.0            cellranger_1.1.0         tools_3.6.1             
##  [58] RSQLite_2.1.3            ranger_0.11.2            evaluate_0.14           
##  [61] stringr_1.4.0            yaml_2.2.0               knitr_1.26              
##  [64] bit64_0.9-7              zip_2.0.4                robustbase_0.93-5       
##  [67] caTools_1.17.1.3         purrr_0.3.3              compiler_3.6.1          
##  [70] rstudioapi_0.10          beeswarm_0.2.3           curl_4.3                
##  [73] e1071_1.7-3              ggsignif_0.6.0           knn.covertree_1.0       
##  [76] smoother_1.1             tibble_2.1.3             statmod_1.4.32          
##  [79] stringi_1.4.3            RSpectra_0.16-0          forcats_0.4.0           
##  [82] lattice_0.20-38          Matrix_1.2-18            vctrs_0.2.0             
##  [85] pillar_1.4.2             lifecycle_0.1.0          lmtest_0.9-37           
##  [88] BiocNeighbors_1.4.1      cowplot_1.0.0            data.table_1.12.6       
##  [91] bitops_1.0-6             irlba_2.3.3              R6_2.4.1                
##  [94] pcaMethods_1.78.0        KernSmooth_2.23-16       gridExtra_2.3           
##  [97] rio_0.5.16               vipor_0.4.5              codetools_0.2-16        
## [100] boot_1.3-23              MASS_7.3-51.4            gtools_3.8.1            
## [103] assertthat_0.2.1         withr_2.1.2              GenomeInfoDbData_1.2.2  
## [106] hms_0.5.2                grid_3.6.1               tidyr_1.0.0             
## [109] class_7.3-15             rmarkdown_1.18           DelayedMatrixStats_1.8.0
## [112] carData_3.0-3            TTR_0.23-5               scatterplot3d_0.3-41    
## [115] ggbeeswarm_0.6.0
```

