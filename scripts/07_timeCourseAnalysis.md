---
title: "<span style='font-size: 28px'>Single-cell RNAseq of mouse heart development</style>"
date: '10 December, 2019'
output:
  html_document:
    keep_md: true
    fig_width: 5
    fig_height: 5
    fig_caption: yes
    code_folding: hide
    toc: false
    toc_float: 
      collapsed: false
---



### Time-course analysis

We have characterised the dataset and the different cell populations present. Samples were collected across several stages of development, so let's now analyse the temporal dynamics of the different cell types.

Cells were collected just before the formation of the cardiac crescent - *stage -1* -, at four different points - *stages 0-3* - from the earliest cardiac crescent to just before the transition to a linear heart tube (LHT) structure, and at the LHT stage. The number of cells for each stage is variable.


```r
## normalised data, with cluster annotation
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))
table(colData(sce)$stage)
```

```
## 
##  -1   0   1   2   3 LHT 
## 227 482 829 950 444 172
```

#### Cell proportions over development

First, we analyse the time of emergence of each subpopulation and their progression across time.

The non-mesodermal clusters are relatively stable throughout the six stages profiled. 

- The ectoderm is almost absent from stage -1, but this is down to how the dissections were done for this stage. Both ectodermal clusters are more prevalent in stage 0 compared to later stages, perhaps due to a smaller cardiac crescent region. Cluster Ec2 only has 89 cells, and thus estimating its contribution per stage is quite unstable, more susceptible to sampling biases.
- The two endodermal clusters show pretty stable contributions across development, with coordinated behaviour from both clusters.
- The blood (Me1) and endothelial (Me2) clusters are also quite stable. As with Ec2, the very small number of cells in each cluster (45 and 59 respectively) make estimating their contributions precisely quite difficult.


```r
props <- as.matrix(table(sce$clusterAnn, sce$stage))
props <- t(t(props)/colSums(props))*100

par(mfrow=c(3:2), mar=c(2,2,2,2))
for(i in 1:6){
  barplot(props[i,], col=cols[i], main=row.names(props)[i])
}
```

![](07_timeCourseAnalysis_files/figure-html/nonCardiac-1.png)<!-- -->

For the cardiac mesoderm clusters:

- Me8, one of the most undifferentiated populations of progenitor cells, is quite specific to stage -1, and then decreases sharply.
- Me7, a progenitor subpopulation, is quite stable across time, and already present abundantly at stage -1. 
- Me6, which links Me7 to the more mature cardiomyocytes of Me3, arises later than Me7, between stages 0 and 1. After, it seems to increase gradually. However, there are only 65 cells in this cluster, so these proportions might not be very accurate.


- The other differentiation trajectory, starting with progenitor cluster Me5 emerges a little later, and seems to start decreasing after stage 2.
- Me4, which links Me5 to the cardiomyocytes of Me3 shows a similar behaviour.


- Finally, mature cardiomyocytes in Me3 appear at stage 0 and increase with time.


```r
layout(matrix(c(0,0,1,1,2,2,3,3,4,4,5,5,0,6,6,0), ncol = 4, byrow = TRUE))
par(mar=c(2,2,2,2))
for(i in c(12,9,11,8,10,7)){
  barplot(props[i,], col=cols[i], main=row.names(props)[i])
}
```

![](07_timeCourseAnalysis_files/figure-html/cardiac-1.png)<!-- -->

We can also use the different replicates (batches) to estimate the mean and standard deviation of the proportions at each stage. We exclude stage -1 from this analysis because it was captured in only one batch. All together, the observations above hold.


```r
## use the differen batches as replicates, to get estimates with mean +- se
props <- list()
for(i in unique(sce$batch)){
  props[[i]] <- as.matrix(table(sce[,sce$batch == i]$clusterAnn, sce[,sce$batch == i]$stage))
  props[[i]] <- t(t(props[[i]])/colSums(props[[i]]))*100
}

## means (exclude stage -1 which has no replicates)
props.mean <- matrix(ncol = 5, nrow = length(cols))
colnames(props.mean) <- paste0("stage_",c(0:3,"LHT"))
row.names(props.mean) <- names(cols)
props.mean[,'stage_0'] <- rowMeans(cbind(props[[3]][,1], props[[7]][,2]))
props.mean[,'stage_1'] <- rowMeans(cbind(props[[3]][,2], props[[4]][,1], props[[5]][,1], props[[6]][,1]))
props.mean[,'stage_2'] <- rowMeans(cbind(c(props[[2]][1:5,1],0,props[[2]][6:11,1]), props[[4]][,2], props[[5]][,2], props[[6]][,2]))
props.mean[,'stage_3'] <- rowMeans(cbind(props[[4]][,3], props[[5]][,3]))
props.mean[,'stage_LHT'] <- rowMeans(cbind(c(props[[1]][1:2,1],0,0,0,props[[1]][3:5,1],0,props[[1]][6:8,1]), c(props[[2]][1:5,2],0,props[[2]][6:11,2])))

## sd (exclude stage -1 which has no replicates)
props.sd <- matrix(ncol = 5, nrow = length(cols))
colnames(props.sd) <- paste0("stage_",c(0:3,"LHT"))
row.names(props.sd) <- names(cols)
props.sd[,'stage_0'] <- rowSds(cbind(props[[3]][,1], props[[7]][,2]))
props.sd[,'stage_1'] <- rowSds(cbind(props[[3]][,2], props[[4]][,1], props[[5]][,1], props[[6]][,1]))
props.sd[,'stage_2'] <- rowSds(cbind(c(props[[2]][1:5,1],0,props[[2]][6:11,1]), props[[4]][,2], props[[5]][,2], props[[6]][,2]))
props.sd[,'stage_3'] <- rowSds(cbind(props[[4]][,3], props[[5]][,3]))
props.sd[,'stage_LHT'] <- rowSds(cbind(c(props[[1]][1:2,1],0,0,0,props[[1]][3:5,1],0,props[[1]][6:8,1]), c(props[[2]][1:5,2],0,props[[2]][6:11,2])))


df <- data.frame(pop = rep(row.names(props.mean),5), stage = rep(substr(colnames(props.mean),7,9), each=nrow(props.mean)), mean = c(props.mean), sd = c(props.sd))

ggplot(df, aes(x=stage, y=mean, group=pop, color=pop)) + geom_line() + geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) + scale_color_manual(values=cols) + facet_wrap(. ~ pop, ncol=4) + ylab("% of cells in cluster") + th + theme(legend.position = "none")
```

![](07_timeCourseAnalysis_files/figure-html/per_batch-1.png)<!-- -->

These dynamics can also be observed by plotting the UMAP representation, stratified by stage.


```r
umap <- reducedDim(sce)
umap$stage <- colData(sce)$stage
umap$cluster <- sce$clusterAnn
umap$col <- sce$clusterCol

stages <- c(-1:3,"LHT")
par(mfrow=c(2,3), mar=c(2,2,2,2))
for(stage in stages){
  plot(umap[umap$stage==stage,]$x, umap[umap$stage==stage,]$y, col=umap[umap$stage==stage,]$col, pch=16, ylim=c(-6.5,11), xlim=c(-10,17), main=paste("stage",stage), xlab="", ylab="", axes=FALSE)
  box(bty="l")
}
```

![](07_timeCourseAnalysis_files/figure-html/umap-1.png)<!-- -->

#### Gene expression across development

Most cell types are present across all six stages profiled. And cells from different stages are generally well mixed in the UMAP representation.


```r
order <- sample(1:nrow(umap), nrow(umap), replace = FALSE)
plot(umap[order,]$x, umap[order,]$y, pch=16, cex=0.75, col=as.factor(umap[order,]$stage), xlab="", ylab="", axes=FALSE)
box(bty="l")
legend("bottomright", legend = levels(as.factor(umap$stage)), pch=16, col=1:6, cex=0.75)
```

![](07_timeCourseAnalysis_files/figure-html/stage-1.png)<!-- -->

This suggests that the transcriptional profiles are determined by the cell type, and their progression along the developmental program (i.e. from progenitors to cardiomyocytes), and not by the particular stage of the embryo.

Nonetheless, we can test whether there are any genes differentially expressed along development. We test all clusters except Ec2, Me1 and Me2, which have very few cells per stage, per cluster. We also exclude Me8, since the vast majority of cells are from stage -1. For all others, we test for differences in expression between stages, while controlling for batch.


```r
resAdj <- list()
clust <- "En1"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
# 47  33 126  99  77  23

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:7)] <- paste0("stage", colnames(design)[-c(1:7)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 8:12)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:7], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
universe <- row.names(y$counts)
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go.res <- list()
go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```



```r
clust <- "En2"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
# 22  28  89  68  43  10

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:7)] <- paste0("stage", colnames(design)[-c(1:7)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 8:12)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:7], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```



```r
clust <- "Ec1"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
#  2  66  31  51  31  13 

## remove stage -1
y$counts <- y$counts[,-which(y$samples$stage==-1)]
y$samples <- y$samples[-which(y$samples$stage==-1),]
stopifnot(identical(row.names(y$samples), colnames(y$counts)))

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:6)] <- paste0("stage", colnames(design)[-c(1:6)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 7:10)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:6], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```



```r
clust <- "Me3"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
#  3  83 261 201  98  67 

## remove stage -1
y$counts <- y$counts[,-which(y$samples$stage==-1)]
y$samples <- y$samples[-which(y$samples$stage==-1),]
stopifnot(identical(row.names(y$samples), colnames(y$counts)))

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:7)] <- paste0("stage", colnames(design)[-c(1:7)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 8:11)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:6], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```



```r
clust <- "Me4"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
#  8  44  77  66  20   6 

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:7)] <- paste0("stage", colnames(design)[-c(1:7)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 8:12)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:7], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```


```r
clust <- "Me5"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
#  4  64  89 150  36  12 

## remove stage -1
y$counts <- y$counts[,-which(y$samples$stage==-1)]
y$samples <- y$samples[-which(y$samples$stage==-1),]
stopifnot(identical(row.names(y$samples), colnames(y$counts)))

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:6)] <- paste0("stage", colnames(design)[-c(1:6)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 7:10)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:6], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```


```r
clust <- "Me6"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
#  0   1   2   3 LHT 
#  2  19  23  13   8

## remove stage -1,0
y$counts <- y$counts[,-which(y$samples$stage %in% c(-1,0))]
y$samples <- y$samples[-which(y$samples$stage%in% c(-1,0)),]
stopifnot(identical(row.names(y$samples), colnames(y$counts)))
y$samples$batch <- droplevels(y$samples$batch)

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:6)] <- paste0("stage", colnames(design)[-c(1:6)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 7:9)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:5], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=1
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```



```r
clust <- "Me7"
tmp <- sce[,sce$clusterAnn == clust]

## filter lowly expressed genes
keep <- rowMeans(logcounts(tmp)) > 1

## convert sce object to a DGElist to use with edgeR
y <- convertTo(tmp[keep,], type="edgeR")
y$samples$batch <- droplevels(y$samples$batch)
# table(y$samples$stage)
# -1   0   1   2   3 LHT 
# 43  85 111 149  97  29

design <- model.matrix(~0+y$samples$batch+y$samples$stage)
colnames(design) <- substr(colnames(design), 16, 30)
colnames(design)[-c(1:7)] <- paste0("stage", colnames(design)[-c(1:7)])

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
test <- glmQLFTest(fit, 8:12)
resAdj[[clust]] <- as.data.frame(topTags(test, n=nrow(y$counts)))
resAdj[[clust]]$logFC <- apply(resAdj[[clust]][,3:7], 1, function(x) x[which.max(abs(x))])
# summary(resAdj[[clust]]$FDR<0.05 & abs(resAdj[[clust]]$logFC) > 1)

# i=100
# boxplot(logcounts(tmp)[row.names(resAdj[[clust]])[i],]~tmp$stage, main=resAdj[[clust]][i,1]); i <- i+1

### test GO enrichments
all <- as.factor(as.numeric(universe %in% row.names(resAdj[[clust]][resAdj[[clust]]$FDR < 0.05 & abs(resAdj[[clust]]$logFC) >1,])))
names(all) <- universe

go <- new("topGOdata", ontology="BP", allGenes = all, nodeSize=5, annot=annFUN.org, mapping="org.Mm.eg.db", ID = "ensembl")
go.test <- runTest(go, algorithm = "classic", statistic = "Fisher" )
go.res[[clust]] <- GenTable(go, Fisher.classic = go.test, topNodes = length(score(go.test)))
go.res[[clust]]$Fisher.classic.adj <- p.adjust(go.res[[clust]]$Fisher.classic, "fdr")
```

A few hundred to a few thousand genes are significantly DE and show a fold-change of 2-fold or higher.


```r
unlist(lapply(resAdj, function(x) sum(x$FDR < 0.05 & abs(x$logFC) > 1)))
```

```
##  En1  En2  Ec1  Me3  Me4  Me5  Me6  Me7 
## 1470  616  899 1147  801  714  348 1970
```

Most of these genes are significant in only one cluster, although nearly a fifth are identified in two different comparisons.


```r
DEgenes <- data.frame(gene=do.call('c', lapply(resAdj, function(x) x[x$FDR < 0.05 & abs(x$logFC) > 1,]$gene)))
DEgenes$pop <- substr(row.names(DEgenes), 1, 3)
DEgenes.mat <- matrix(nrow = length(unique(DEgenes$gene)), ncol = length(unique(DEgenes$pop)))
row.names(DEgenes.mat) <- unique(DEgenes$gene)
colnames(DEgenes.mat) <- unique(DEgenes$pop)
for(clust in colnames(DEgenes.mat)){
  DEgenes.mat[,clust] <- ifelse(row.names(DEgenes.mat) %in% DEgenes[DEgenes$pop == clust,1], 1, 0)
}

upset(as.data.frame(DEgenes.mat), nsets = 8)
```

![](07_timeCourseAnalysis_files/figure-html/upsetDE-1.png)<!-- -->

The sets of differential genes are enriched for few GO terms, and these are more often shared between clusters.


```r
GOterms <- data.frame(GO=do.call('c', lapply(go.res, function(x) x[x$Fisher.classic.adj < 0.05,]$GO.ID)))
GOterms$pop <- substr(row.names(GOterms), 1, 3)

GOterms.mat <- matrix(nrow = length(unique(GOterms$GO)), ncol = length(unique(GOterms$pop)))
row.names(GOterms.mat) <- unique(GOterms$GO)
colnames(GOterms.mat) <- unique(GOterms$pop)
for(clust in colnames(GOterms.mat)){
  GOterms.mat[,clust] <- ifelse(row.names(GOterms.mat) %in% GOterms[GOterms$pop == clust,1], 1, 0)
}
upset(as.data.frame(GOterms.mat), nsets = 8)
```

![](07_timeCourseAnalysis_files/figure-html/upseGO-1.png)<!-- -->

```r
## large proportion of terms shared by all clusters
## these relate to translation, biosynthetic/metabolic processes, ribosomal/ribonucleprotein terms; mithochondrial respiratory chain/ATP synthesis, oxidative phosphorylation
```

This suggests that, although the sets of genes are different between clusters, many tend to be involved in related processes. And these processes are related to translation and other core bio-metabolic pathways.


```r
go.res[['En1']][go.res[['En1']]$GO.ID %in% row.names(GOterms.mat[rowSums(GOterms.mat)>4,]),2]
```

```
##  [1] "translation"                                
##  [2] "peptide biosynthetic process"               
##  [3] "amide biosynthetic process"                 
##  [4] "peptide metabolic process"                  
##  [5] "cytoplasmic translation"                    
##  [6] "cellular amide metabolic process"           
##  [7] "organonitrogen compound biosynthetic pro..."
##  [8] "ribosomal small subunit assembly"           
##  [9] "ribosome assembly"                          
## [10] "cellular protein-containing complex asse..."
## [11] "NADH dehydrogenase complex assembly"        
## [12] "mitochondrial respiratory chain complex ..."
## [13] "ATP metabolic process"                      
## [14] "mitochondrial respiratory chain complex ..."
## [15] "ribonucleoside monophosphate metabolic p..."
## [16] "ribosomal small subunit biogenesis"         
## [17] "nucleoside triphosphate metabolic proces..."
## [18] "purine nucleoside monophosphate metaboli..."
## [19] "purine ribonucleoside monophosphate meta..."
## [20] "ribonucleoside triphosphate metabolic pr..."
## [21] "purine ribonucleoside triphosphate metab..."
## [22] "antimicrobial humoral immune response me..."
## [23] "nucleoside monophosphate metabolic proce..."
## [24] "ribonucleoprotein complex assembly"         
## [25] "oxidative phosphorylation"                  
## [26] "ATP synthesis coupled electron transport"   
## [27] "purine nucleoside triphosphate metabolic..."
## [28] "drug metabolic process"                     
## [29] "mitochondrial ATP synthesis coupled elec..."
## [30] "ribonucleoprotein complex subunit organi..."
## [31] "protein-containing complex subunit organ..."
## [32] "protein-containing complex assembly"        
## [33] "respiratory electron transport chain"       
## [34] "ribose phosphate metabolic process"         
## [35] "generation of precursor metabolites and ..."
```

This suggests that, over development, the metabolic activity and protein production rate is altered, but no other processes are significantly regulated.


```r
saveRDS(resAdj, paste0(dir, "results/DE_stage_analysis.Rds"))
saveRDS(go.res, paste0(dir, "results/DE_stage_analysis_GOenrichments.Rds"))
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
##  [1] org.Mm.eg.db_3.10.0         UpSetR_1.4.0               
##  [3] topGO_2.37.0                SparseM_1.77               
##  [5] GO.db_3.10.0                AnnotationDbi_1.48.0       
##  [7] graph_1.64.0                ggpubr_0.2.4               
##  [9] magrittr_1.5                ggplot2_3.2.1              
## [11] RColorBrewer_1.1-2          edgeR_3.28.0               
## [13] limma_3.42.0                scran_1.14.5               
## [15] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
## [17] DelayedArray_0.12.0         BiocParallel_1.20.0        
## [19] matrixStats_0.55.0          Biobase_2.46.0             
## [21] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [23] IRanges_2.20.1              S4Vectors_0.24.1           
## [25] BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6             bit64_0.9-7              tools_3.6.1             
##  [4] backports_1.1.5          R6_2.4.1                 irlba_2.3.3             
##  [7] vipor_0.4.5              DBI_1.0.0                lazyeval_0.2.2          
## [10] colorspace_1.4-1         withr_2.1.2              tidyselect_0.2.5        
## [13] gridExtra_2.3            bit_1.1-14               compiler_3.6.1          
## [16] BiocNeighbors_1.4.1      labeling_0.3             scales_1.1.0            
## [19] stringr_1.4.0            digest_0.6.23            rmarkdown_1.18          
## [22] XVector_0.26.0           scater_1.14.4            pkgconfig_2.0.3         
## [25] htmltools_0.4.0          rlang_0.4.2              rstudioapi_0.10         
## [28] RSQLite_2.1.3            DelayedMatrixStats_1.8.0 farver_2.0.1            
## [31] dplyr_0.8.3              RCurl_1.95-4.12          BiocSingular_1.2.0      
## [34] GenomeInfoDbData_1.2.2   Matrix_1.2-18            Rcpp_1.0.3              
## [37] ggbeeswarm_0.6.0         munsell_0.5.0            viridis_0.5.1           
## [40] lifecycle_0.1.0          stringi_1.4.3            yaml_2.2.0              
## [43] zlibbioc_1.32.0          plyr_1.8.4               grid_3.6.1              
## [46] blob_1.2.0               dqrng_0.2.1              crayon_1.3.4            
## [49] lattice_0.20-38          splines_3.6.1            locfit_1.5-9.1          
## [52] zeallot_0.1.0            knitr_1.26               pillar_1.4.2            
## [55] igraph_1.2.4.2           ggsignif_0.6.0           glue_1.3.1              
## [58] evaluate_0.14            vctrs_0.2.0              gtable_0.3.0            
## [61] purrr_0.3.3              assertthat_0.2.1         xfun_0.11               
## [64] rsvd_1.0.2               viridisLite_0.3.0        tibble_2.1.3            
## [67] beeswarm_0.2.3           memoise_1.1.0            statmod_1.4.32
```

