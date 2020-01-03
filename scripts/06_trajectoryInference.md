---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '18 December, 2019'
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

And defining branches based on diffusion pseudotime separates the cardiomyocytes (Me3) from the yellow cluster (Me6) and the orange (Me4) + green (Me5) progenitors. 


```r
table(dpt.me5@branch[,1], sce.Me5[,row.names(dm.me5)]$clusterAnn)
```

```
##    
##     Me3 Me4 Me5 Me6
##   1  39 200 355   2
##   2 609   0   0   4
##   3   7   0   0  44
```

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

Consistently, there is a branching point that separates the most distinct orange cells (Me4); but some orange cells are part of the branch containing the turquoise (Me7-8) and yellow (Me6) progenitors.


```r
table(dpt.me7@branch[,1], sce.Me7[,row.names(dm.me7)]$clusterAnn)
```

```
##    
##     Me3 Me4 Me6 Me7 Me8
##   1 602  35   0   0   0
##   2  40  59  60 514 184
##   3   0  90   0   0   0
```

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

We identify 726 genes dynamically expressed along this trajectory. These can be clustered into five different profiles.


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

clust.genes_me5 <- cutreeDynamic(test.clust,distM=as.matrix(test.dist), minClusterSize=50, method="hybrid", deepSplit = 1, verbose = 0)  # reducing the min cluster size doesn't split the tree further
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
ggarrange(plotlist = plots, ncol = 3, nrow = 2)
```

![](06_trajectoryInference_files/figure-html/dynGenes_me6_clusters-1.png)<!-- -->

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
go.res.me5[c(1,2,6,7,8,20,23,25,27,41,46,58,63,64,93,95),-6]
```

```
##         GO.ID                               Term Annotated Significant Expected
## 1  GO:0007507                  heart development       512         109    28.99
## 2  GO:0061061       muscle structure development       530         110    30.01
## 6  GO:0072359     circulatory system development       857         138    48.53
## 7  GO:0014706 striated muscle tissue development       345          83    19.54
## 8  GO:0048738  cardiac muscle tissue development       220          65    12.46
## 20 GO:0030239                 myofibril assembly        60          30     3.40
## 23 GO:0035051         cardiocyte differentiation       148          45     8.38
## 25 GO:0006936                 muscle contraction       194          51    10.99
## 27 GO:0003206      cardiac chamber morphogenesis       125          40     7.08
## 41 GO:0030036    actin cytoskeleton organization       525          82    29.73
## 46 GO:0045214             sarcomere organization        43          22     2.43
## 58 GO:0008015                  blood circulation       309          55    17.50
## 63 GO:0060485             mesenchyme development       214          43    12.12
## 64 GO:0007155                      cell adhesion       765          94    43.32
## 93 GO:0001944            vasculature development       549          70    31.09
## 95 GO:0008283                 cell proliferation      1318         129    74.63
##      Fisher.classic.adj
## 1               < 1e-30
## 2               < 1e-30
## 6               < 1e-30
## 7               < 1e-30
## 8           5.84156e-26
## 20 2.97543076923077e-19
## 23            7.894e-19
## 25 1.92964444444444e-18
## 27          1.46039e-17
## 41 3.48264705882353e-15
## 46 1.82169230769231e-14
## 58 2.78611764705882e-12
## 63 2.81928571428571e-11
## 64 4.43171929824561e-11
## 93 8.16937209302326e-09
## 95 1.25586363636364e-08
```

Considering significantly enriched terms (FDR < 1%) that have at least 20 dynamic genes, we look at the proportion of genes in each of the clusters defined above. Considering the terms where the proportion for one the clusters deviates considerably from the overall behaviour, we observe that terms with a much higher proportion of genes in cluster 1, which are increasing expression along pseudotime, include the development of cardiac and striated muscle, heart and circulatory system development. This makes sense since pseudotime increases as cells mature into cardiomyocytes. In contrast, terms enriched for cluster 2 genes, which are decreasing along pseudotime, include regulation of the cell cycle and proliferation. Terms with larger proportions of genes from cluster 3, which are downregulated in the middle of pseudotime progression, include several signalling pathways. 


```r
sig <- go.res.me5[go.res.me5$Fisher.classic.adj < 0.01 & go.res.me5$Significant >= 20,]

props.me5 <- matrix(nrow = length(unique(clust.genes_me5)), ncol = nrow(sig))
colnames(props.me5) <- sig$GO.ID
row.names(props.me5) <- 1:5
for(go in sig$GO.ID){
  g <- intersect(genesInTerm(go.me5, go)[[1]], names(dynGenes_me5))
  tmp <- as.matrix(table(clust.genes_me5[g])/length(g)*100)
  props.me5[row.names(tmp),go] <- tmp
}
props.me5[is.na(props.me5)] <- 0

outliers <- c()
for(i in 1:5){
  outliers <- c(outliers, which(props.me5[i,] < median(props.me5[i,])-mad(props.me5[i,]) | props.me5[i,] > median(props.me5[i,])+mad(props.me5[i,])))
}

tmp <- props.me5[,unique(outliers)]
tmp <- tmp[,order(tmp[1,])]
colnames(tmp) <- go.res.me5[match(colnames(tmp), go.res.me5$GO.ID),]$Term
# colnames(tmp)
sel <- c(1:4,17,23,26,38,53,63,70,81,85,88,89,105,113,114,115)

par(mar=c(12,2,2,2))
barplot(tmp[,sel], col=1:5, las=2, cex.names =0.75, xlim=c(0,25))
legend("topright", legend = 1:5, col = 1:5, pch = 15, title = "cluster", cex=0.85)
```

![](06_trajectoryInference_files/figure-html/go_me5_profiles-1.png)<!-- -->

##### Me7-8 -> Me3 trajectory

Now we repeat the same for the trajectory from the Me7-8 progenitors.


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

1110 genes are dynamically expressed along this pathway, and can be clustered into six profiles. 


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

clust.genes_me7 <- cutreeDynamic(test.clust,distM=as.matrix(test.dist), minClusterSize=50, method="hybrid", deepSplit = 1, verbose = 0)  ## reducing min cluster size doesn't split the tree further
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
ggarrange(plotlist = plots, ncol = 3, nrow = 2)
```

![](06_trajectoryInference_files/figure-html/dynGenes_me7_clusters-1.png)<!-- -->

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
```

And similarly, terms related to muscle and heart development have large proportion of genes that are upregulated as pseudotime increases, whilst terms related to the cell cycle and proliferation have the opposite behaviour.


```r
sig <- go.res.me7[go.res.me7$Fisher.classic.adj < 0.01 & go.res.me7$Significant >= 20,]

props.me7 <- matrix(nrow = length(unique(clust.genes_me7)), ncol = nrow(sig))
colnames(props.me7) <- sig$GO.ID
row.names(props.me7) <- 1:6
for(go in sig$GO.ID){
  g <- intersect(genesInTerm(go.me7, go)[[1]], names(dynGenes_me7))
  tmp <- as.matrix(table(clust.genes_me7[g])/length(g)*100)
  props.me7[row.names(tmp),go] <- tmp
}
props.me7[is.na(props.me7)] <- 0

outliers <- c()
for(i in 1:6){
  outliers <- c(outliers, which(props.me7[i,] < median(props.me7[i,])-mad(props.me7[i,]) | props.me7[i,] > median(props.me7[i,])+mad(props.me7[i,])))
}

tmp <- props.me7[,unique(outliers)]
tmp <- tmp[,order(tmp[1,])]
colnames(tmp) <- go.res.me7[match(colnames(tmp), go.res.me7$GO.ID),]$Term
# colnames(tmp)
sel <- c(2,4,8,30,48,51,85,87,95,102,110,113,121)

par(mar=c(12,2,2,2))
barplot(tmp[,sel], col=1:6, las=2, cex.names =0.75, xlim=c(0,16))
legend("topright", legend = 1:6, col = 1:6, pch = 15, title = "cluster")
```

![](06_trajectoryInference_files/figure-html/go_me7_profiles-1.png)<!-- -->

#### Relationship with developmental stages

The diffusion map orders cells based on the similarity of their transcriptomes. And given the paths that have been revealed, we infer that the algorithm is capturing the progression from the most naïve progenitors to the most mature cardiomyocytes. 

But is there any relationship between the differentiation dynamics and the developmental stage of the embryos?

The cells from embryos of different stages are arranged quite uniformly along the diffusion space, with the exception of the cells fro stage -1, which preferentially occupy the earliest stages of the turquoise (Me7-8 -> Me6(Me4) -> Me3) trajectory and are pretty much absent from the green (Me5 -> Me4 -> Me3) path.


```r
dm$stage <- sce[,row.names(dm)]$stage
ggplot(dm, aes(DC1, DC2, colour=stage)) + geom_point(size=1, alpha=0.5) + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + facet_wrap(~stage) + th + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
```

![](06_trajectoryInference_files/figure-html/dm_perStage-1.png)<!-- -->

And by uniformly I mean consistent with the density of cells observed in the different areas of the diffusion space (more cells at the tips, especially in the turquoise and blue clusters). Thus, all stages but -1 follow similar distributions of density across the first two diffusion components.


```r
plots <- list()
plots[[1]] <- ggplot(dm, aes(DC1, colour=stage)) + geom_density() + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
plots[[2]] <- ggplot(dm, aes(DC2, colour=stage)) + geom_density() + scale_color_manual(values = brewer.pal(n=8, "YlOrBr")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
ggarrange(plotlist = plots, ncol = 2)
```

![](06_trajectoryInference_files/figure-html/density_perStage-1.png)<!-- -->

This suggests that from the formation of the cardiac crescent at stage 0 onward, the cells do not show substantial transcriptional changes based on embryo age; instead, they are part of a continuous differentiation process that stays constant across these developmental stages.

#### Visualisation

Finally, for ease of interpretation, we plot the diffusion map so that the progenitor populations are at the top and the cardiomyocytes at the bottom (by plotting DC1 on the y-axis and DC2 on the x-axis), which is a more intuitive interpretation.


```r
ggplot(dm, aes(DC2, DC1, colour=sce.cardiac$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
```

![](06_trajectoryInference_files/figure-html/vis-1.png)<!-- -->



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
## [11] ggpubr_0.2.4                magrittr_1.5               
## [13] ggplot2_3.2.1               destiny_3.0.0              
## [15] scran_1.14.5                SingleCellExperiment_1.8.0 
## [17] SummarizedExperiment_1.16.0 DelayedArray_0.12.0        
## [19] BiocParallel_1.20.0         matrixStats_0.55.0         
## [21] Biobase_2.46.0              GenomicRanges_1.38.0       
## [23] GenomeInfoDb_1.22.0         IRanges_2.20.1             
## [25] S4Vectors_0.24.1            BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1             backports_1.1.5          RcppEigen_0.3.3.7.0     
##   [4] igraph_1.2.4.2           lazyeval_0.2.2           sp_1.3-2                
##   [7] RcppHNSW_0.2.0           scater_1.14.4            digest_0.6.23           
##  [10] htmltools_0.4.0          viridis_0.5.1            gdata_2.18.0            
##  [13] memoise_1.1.0            openxlsx_4.1.3           limma_3.42.0            
##  [16] xts_0.11-2               colorspace_1.4-1         blob_1.2.0              
##  [19] haven_2.2.0              xfun_0.11                dplyr_0.8.3             
##  [22] crayon_1.3.4             RCurl_1.95-4.12          hexbin_1.28.0           
##  [25] zeallot_0.1.0            zoo_1.8-6                glue_1.3.1              
##  [28] gtable_0.3.0             zlibbioc_1.32.0          XVector_0.26.0          
##  [31] car_3.0-5                BiocSingular_1.2.0       DEoptimR_1.0-8          
##  [34] abind_1.4-5              VIM_4.8.0                scales_1.1.0            
##  [37] ggplot.multistats_1.0.0  DBI_1.0.0                edgeR_3.28.0            
##  [40] ggthemes_4.2.0           Rcpp_1.0.3               viridisLite_0.3.0       
##  [43] laeken_0.5.0             dqrng_0.2.1              foreign_0.8-72          
##  [46] bit_1.1-14               rsvd_1.0.2               proxy_0.4-23            
##  [49] vcd_1.4-4                farver_2.0.1             pkgconfig_2.0.3         
##  [52] nnet_7.3-12              labeling_0.3             tidyselect_0.2.5        
##  [55] rlang_0.4.2              munsell_0.5.0            cellranger_1.1.0        
##  [58] tools_3.6.1              RSQLite_2.1.3            ranger_0.11.2           
##  [61] evaluate_0.14            stringr_1.4.0            yaml_2.2.0              
##  [64] knitr_1.26               bit64_0.9-7              zip_2.0.4               
##  [67] robustbase_0.93-5        caTools_1.17.1.3         purrr_0.3.3             
##  [70] compiler_3.6.1           rstudioapi_0.10          beeswarm_0.2.3          
##  [73] curl_4.3                 e1071_1.7-3              ggsignif_0.6.0          
##  [76] knn.covertree_1.0        smoother_1.1             tibble_2.1.3            
##  [79] statmod_1.4.32           stringi_1.4.3            RSpectra_0.16-0         
##  [82] forcats_0.4.0            lattice_0.20-38          Matrix_1.2-18           
##  [85] vctrs_0.2.0              pillar_1.4.2             lifecycle_0.1.0         
##  [88] lmtest_0.9-37            BiocNeighbors_1.4.1      cowplot_1.0.0           
##  [91] data.table_1.12.6        bitops_1.0-6             irlba_2.3.3             
##  [94] R6_2.4.1                 pcaMethods_1.78.0        KernSmooth_2.23-16      
##  [97] gridExtra_2.3            rio_0.5.16               vipor_0.4.5             
## [100] codetools_0.2-16         boot_1.3-23              MASS_7.3-51.4           
## [103] gtools_3.8.1             assertthat_0.2.1         withr_2.1.2             
## [106] GenomeInfoDbData_1.2.2   hms_0.5.2                grid_3.6.1              
## [109] tidyr_1.0.0              class_7.3-15             rmarkdown_1.18          
## [112] DelayedMatrixStats_1.8.0 carData_3.0-3            TTR_0.23-5              
## [115] scatterplot3d_0.3-41     ggbeeswarm_0.6.0
```

