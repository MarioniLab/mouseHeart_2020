###### Figures ######
library(scran)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(locfit)
library(reshape2)
library(destiny)
library(gplots)

dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/mouseHeart_2020/"
out <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/WRITING/HEART/Figures/figureElements/"

th <- theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size=10), axis.title.x = element_text(size=12), axis.text.y = element_text(size=10), axis.title.y = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5))


### Figure 1 ###############

## UMAP with clusters (1C) ======
umap <- read.table(paste0(dir, "results/umapCoords_corrected.tab"))
clusters <- read.table(paste0(dir, "results/clusters_average_min40.tsv"), row.names = 1)
stopifnot(identical(row.names(umap), row.names(clusters)))

# annotate the clusters
ann <- c("En1","En2","Ec1","Ec2",paste0("Me",1:8)) # En, Ec and Me for the three germ layers
names(ann) <- c(3, 5, 7, 9, 12, 11, 1, 6, 4, 10, 2, 8) ## map the cluster number to the annotation
ann <- ann[order(as.numeric(names(ann)))]

# colors
cols <- c(Ec1 = "#ec6646", Ec2 = "#af4424", En1 = "#3c537c", En2 = "#768ba5",
          Me1 = "#bf9a77", Me2 = "#debc95", Me3 = "#556dad", Me4 = "#f28a31", 
          Me5 = "#729f3c", Me6 = "#fbba14", Me7 = "#5fa398", Me8 = "#9FD3C5")

# remove outlier
outlier <- which(clusters$V2==0)
umap <- umap[-outlier,]
clusters <- clusters[-outlier,]

# add cluster colors
clusters <- data.frame(cluster = clusters, ann = ann[clusters], row.names = row.names(umap), stringsAsFactors = FALSE)
clusters$col <- cols[clusters$ann]

order <- sample(1:nrow(umap), nrow(umap), replace = FALSE)
pdf(paste0(out, "Fig1_UMAP.pdf"), width = 7, height = 7, useDingbats = FALSE)
plot(umap$x[order], umap$y[order], 
     pch=16, col=clusters$col[order], 
     xlab="UMAP dim1", ylab="UMAP dim2", 
     axes=FALSE)
box(bty="l")
legend("bottomright", legend = ann[order(ann)], 
       col=cols[ann[order(ann)]], 
       pch=16, cex=0.5, ncol=2)
dev.off()


## marker heatmap (1D,H)======
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))
markers <- c("Ttr","Hhex","Foxa2","Cdh1","Sox17","Pax9","Sox2","Dlx5","Wnt6","Tfap2a","Sox1","Gata1","Hba-a1","Tal1","Sox7","Cdh5","Sox17","Pdgfra","Nkx2-5","Mef2c","Gata4","Tnnt2","Ttn","Myh6","Actn2",
             "Nkx2-5","Hand1","Tbx5","Hcn4","Sfrp5",
             "Isl1","Foxc2","Tbx1","Fgf8","Hoxa1","Hoxb1")
tmp <- logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene %in% markers,]),]
tmp <- t(apply(tmp, 1, function(x) x/max(x)))
row.names(tmp) <- rowData(sce)[row.names(tmp),]$gene
tmp <- tmp[match(markers, row.names(tmp)),]

palette <- brewer.pal(n=9, name="Purples")[-9]

clust <- hclust(dist(t(tmp)))
order <- colnames(tmp)[clust$order]
order <- clusters[order,]
order <- order[order(order$ann),]
order <- rbind(order[order$ann=="En2",], order[order$ann=="En1",], 
               order[order$ann=="Ec1",], order[order$ann=="Ec2",],
               order[order$ann=="Me1",], order[order$ann=="Me2",],
               order[order$ann=="Me5",], order[order$ann=="Me4",], 
               order[order$ann=="Me8",], order[order$ann=="Me7",], order[order$ann=="Me6",], 
               order[order$ann=="Me3",])

order$split <- ifelse(grepl("En", order$ann), 1, ifelse(grepl("Ec", order$ann), 2, ifelse(grepl("Me1", order$ann), 3, ifelse(grepl("Me2", order$ann), 3, 4))))

gene.split <- c(rep(1,6), rep(2,5), rep(3,3), rep(4,3), rep(5,8), 6, rep(7,4), rep(8,6))
ha = HeatmapAnnotation(df = data.frame(cluster = order$ann), col = list(cluster = cols))

pdf(paste0(out, "Fig1_heatmap_markers.pdf"), width = 7, height = 7, useDingbats = FALSE)
Heatmap(tmp[,row.names(order)], col=palette, 
        top_annotation = ha, 
        cluster_columns = FALSE, cluster_rows = FALSE, 
        show_column_names = FALSE, 
        column_split = order$split, 
        row_split = gene.split)
dev.off()


## Cell-cycle phase assignment (1G)======
cell_cycle <- read.table(paste0(dir, "results/cellCyclePhase.tsv"), 
                         stringsAsFactors = FALSE, row.names = 1)
stopifnot(identical(row.names(umap), row.names(cell_cycle)))

palette(c("#626D71", "#B38867", "#DDBC95"))

order <- sample(1:nrow(umap), nrow(umap), replace = FALSE)
pdf(paste0(out, "Fig1_cellCycle.pdf"), width = 7, height = 6, useDingbats = FALSE)
plot(umap$x[order], umap$y[order], pch=16, cex=0.75, 
     col=as.factor(cell_cycle$V2[order]), 
     bty="l", xlab="", ylab="", axes = FALSE)
box(bty="l")
legend("bottomright", legend = levels(as.factor(cell_cycle$V2[order])), 
       col=1:3, pch=16, cex=0.75)
dev.off()


### Figure 2 ###############

## Reference data UMAP (2B)======
sce.ref <- readRDS(paste0(dir, "data/sce_referenceCells_goodQual_clean.NORM.clusters.Rds"))
umap.ref <- as.data.frame(reducedDim(sce.ref))
umap.ref$cluster <- sce.ref$cluster
umap.ref$ann <- sce.ref$regionAnn

cols.ref <- c(cranialMesoderm = "#75B3E2", cardiacMesoderm = "#ED6B58", 
              dorsalMesoderm = "#A57CB5", caudalPSM = "#589E46")
cols.ref.cluster <- c(cluster3 = "#25958C", cluster1 = "#6BB288", cluster2 = "#F4CC71", 
                      cluster6 = "#DF7A24", cluster5 = "#9C5016")

pdf(paste0(out, "Fig2_UMAP_refData.pdf"), width = 10, height = 5, useDingbats = FALSE)
par(mfrow=c(1,2))
plot(umap.ref$V2, umap.ref$V1, axes=FALSE, 
     xlab="UMAP_1", ylab="UMAP_2", 
     col="black", bg=cols.ref[umap.ref$ann], pch=21)
box(bty="l")
legend("topleft", legend = names(cols.ref), pch=16, col=cols.ref)
plot(umap.ref$V2, umap.ref$V1, axes=FALSE, 
     xlab="UMAP_1", ylab="UMAP_2", 
     col="black", bg=cols.ref.cluster[umap.ref$cluster], pch=21)
box(bty="l")
legend("topleft", legend = names(cols.ref.cluster), pch=16, col=cols.ref.cluster)
dev.off()


## Markers of reference regions (2C)======
genes <- c("Cyp26c1", "Tcf15", "Hoxa1", "Meox1", "Gata5", "Ttn")
plots <- list()
for(gene in genes){
  tmp <- data.frame(x=umap.ref$V1, y=umap.ref$V2, 
                    expr=as.numeric(logcounts(sce.ref[row.names(rowData(sce.ref)[rowData(sce.ref)$gene == gene,]),])))
  plots[[gene]] <- ggplot(tmp, aes(y, x, col=expr)) + 
    geom_point() + 
    xlab("") + ylab("") + 
    ggtitle(gene) +
    scale_color_gradientn(colors=c("grey", brewer.pal(n=9, "Blues")[-1])) +
    th + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
}
pdf(paste0(out, "Fig2_importantGenesForest.pdf"), width = 12, height = 8, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=3, nrow=2, common.legend = TRUE)
dev.off()


## Classification of atlas (2D)======
labels <- read.table(paste0(dir, "results/classesUnbiasedMe3-8.randForest.tsv"))
probs <- read.table(paste0(dir, "results/classesUnbiasedMe3-8.randForest.probs.tsv"))

calls <- data.frame(class=labels, max.prob=apply(probs, 2, max), 
                    closest=apply(probs, 2, function(x) max(x[x!=max(x)]) ))
calls$diff <- calls$max.prob-calls$closest
calls$pass <- ifelse(calls$diff > 0.15, 1, 0)

umap$class <- labels[match(row.names(umap), labels$V1),2]
umap$pass <- ifelse(row.names(umap) %in% row.names(calls[calls$pass==1,]), 1, 0)
umap$cluster <- clusters[row.names(umap),2]

palette(brewer.pal(n=8, "Set2"))

order <- sample(row.names(umap[umap$pass==1,]))
pdf(paste0(out, "Fig2_atlasClassification.pdf"), useDingbats = FALSE)
plot(umap$x, umap$y, pch=16, col="lightgrey", 
     axes=FALSE, xlab="", ylab="", 
     xlim=c(-12,0), ylim=c(-7,4))
points(umap[umap$pass==0 & !is.na(umap$class),]$x, 
       umap[umap$pass==0 & !is.na(umap$class),]$y, 
       pch=16, col="grey")
points(umap[umap$pass==1,][order,]$x, umap[umap$pass==1,][order,]$y, 
       pch=16, col=umap[umap$pass==1,][order,]$class)
box(bty="l")
legend("topleft", legend = levels(as.factor(umap$class)), pch=16, col=1:6)
dev.off()


## Proportions of classes per cluster (2E)======
tmp <- umap[umap$pass==1,] ## remove low-confidence calls
freqs <- table(tmp$class, tmp$cluster)
props <- prop.table(freqs, 2)*100

freqs <- melt(freqs)
props <- melt(props)
props$Var1 <- factor(props$Var1, levels = c("cranialMesoderm", "caudalPSM_13", 
                                            "dorsalMesoderm_1", "dorsalMesoderm_2", 
                                            "cardiacMesoderm_2", "cardiacMesoderm_56"))

pdf(paste0(out, "Fig2_atlasClassification_props.pdf"), useDingbats = FALSE)
ggplot(props, aes(Var2, value, fill = Var1, label = freqs$value)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  xlab("") + ylab("") + labs(fill="") +
  th
dev.off()


## Expression of marker genes (2F)======
genes <- c("Tbx1", "Fst", "Foxc2", "Nkx2-5")

tmp <- logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene %in% genes,]), ]
row.names(tmp) <- rowData(sce)[row.names(tmp),]$gene

plots <- list()
for(gene in row.names(tmp)){
  df <- data.frame(umap, expr=tmp[gene,])
  df <- df[order(df$expr),]
  plots[[gene]] <- ggplot(df, aes(x, y, colour=expr)) + geom_point() + 
    scale_color_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + 
    xlab("") + ylab("") + xlim(-9.5,-2.5) + ylim(-6.5,3) + 
    ggtitle(gene) + th + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                               axis.ticks.y = element_blank())
}

pdf(paste0(out, "Fig2_markers.pdf"), width = 12, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 4, common.legend = TRUE)
dev.off()

## Diffusion map (2H)======
diffMap <- readRDS(paste0(dir, "results/diffusionMap_cardiacMesoderm.Rds"))
diffMap <- as.data.frame(diffMap@eigenvectors[,1:5])
diffMap$cluster <- clusters[row.names(diffMap),]$ann

pdf(paste0(out, "Fig2_diffMap.pdf"), useDingbats = FALSE)
ggplot(diffMap, aes(DC2, -DC1, colour=cluster)) + 
  geom_point() + 
  scale_color_manual(values = cols[7:12]) + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank())
dev.off()


## Diffusion pseudotime
dpt <- readRDS(paste0(dir, "results/diffusionPseudotime_cardiacMesoderm.Rds"))
diffMap$dpt <- dpt$DPT684

pdf(paste0(out, "Fig2_diffPDT.pdf"), useDingbats = FALSE)
ggplot(diffMap, aes(DC2, -DC1, colour=dpt)) + 
  geom_point() + 
  scale_color_gradientn(colours = rev(brewer.pal(n=9, "YlOrBr")[-c(1)])) + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank())
dev.off()



### Figure 5 ###############

## Diffusion pseudotime for Me5->M3
dpt.me5 <- read.table(paste0(dir, "results/diffusionPseudotime_Me5.tsv"), row.names = 1)

## Nkx2-5 expr (5B)======
df <- data.frame(x = diffMap$DC2, y = -diffMap$DC1, 
                 log2.exp = logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene == "Nkx2-5",]),
                                           row.names(diffMap)])
df$me5 <- dpt.me5[match(row.names(df), row.names(dpt.me5)),1]
fit2 <- locfit(df[!is.na(df$me5),3]~lp(df[!is.na(df$me5),4], nn=1, deg=2), data=df)
curve <- predict(fit2, df[!is.na(df$me5),4])

pdf(paste0(out, "Fig5_Nkx2-5_Me5-3.pdf"), useDingbats = FALSE)
ggplot(df[!is.na(df$me5),], aes(me5, log2.exp)) + 
  geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me5),]),]$col)) + 
  ylim(c(0,max(df$me5))) + 
  ggtitle("Nkx2-5") + 
  xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + 
  th + 
  theme(axis.text.x = element_blank()) + 
  geom_line(aes(df[!is.na(df$me5),]$me5, curve), lwd=1)
dev.off()


## Mab21l2 expr (5E)======
df <- data.frame(x = diffMap$DC2, y = -diffMap$DC1, 
                 log2.exp = logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene == "Mab21l2",]),
                                           row.names(diffMap)])
df$me5 <- dpt.me5[match(row.names(df), row.names(dpt.me5)),1]
fit2 <- locfit(df[!is.na(df$me5),3]~lp(df[!is.na(df$me5),4], nn=1, deg=2), data=df)
curve <- predict(fit2, df[!is.na(df$me5),4])

pdf(paste0(out, "Fig5_Mab21l2_Me5-3.pdf"), useDingbats = FALSE)
ggplot(df[!is.na(df$me5),], aes(me5, log2.exp)) + 
  geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me5),]),]$col)) + 
  ylim(c(0,max(df$me5))) + 
  ggtitle("Mab21l2") + 
  xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + 
  th + 
  theme(axis.text.x = element_blank()) + 
  geom_line(aes(df[!is.na(df$me5),]$me5, curve), lwd=1)
dev.off()


## Tbx18+ cells (5F)======
df <- data.frame(cluster=sce$clusterAnn, 
                 expr=logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene == "Tbx18",]),])
df$positive <- ifelse(df$expr>0, 1, 0)
props <- as.data.frame(prop.table(table(df$cluster, df$positive)[-c(1:6),], 1)*100)
props <- props[props$Var2==1,]
props



### Figure S1 ###############

## UMAP per stage (S1A)======
umap$stage <- colData(sce)$stage
umap$cluster <- sce$clusterAnn
umap$col <- sce$clusterCol

pdf(paste0(out, "FigS1_UMAPstages.pdf"), width = 10, height = 4.5, useDingbats = FALSE)
stages <- c(-1:3,"LHT")
par(mfrow=c(2,3), mar=c(2,2,2,2))
for(stage in stages){
  plot(umap[umap$stage==stage,]$x, umap[umap$stage==stage,]$y, col=umap[umap$stage==stage,]$col, 
       pch=16, ylim=c(-6.5,11), xlim=c(-10,17), main=paste("stage",stage), 
       xlab="", ylab="", axes=FALSE)
  box(bty="l")
}
dev.off()


## Proportion of cells per cluster across stages (S1B)======
# use the different batches as replicates, to get estimates with mean +- sd
props <- list()
for(i in unique(sce$batch)){
  props[[i]] <- as.matrix(table(sce[,sce$batch == i]$clusterAnn, sce[,sce$batch == i]$stage))
  props[[i]] <- t(t(props[[i]])/colSums(props[[i]]))*100
}

# means (exclude stage -1 which has no replicates)
props.mean <- matrix(ncol = 5, nrow = length(cols))
colnames(props.mean) <- paste0("stage_",c(0:3,"LHT"))
row.names(props.mean) <- names(cols)
props.mean[,'stage_0'] <- rowMeans(cbind(props[[3]][,1], props[[7]][,2]))
props.mean[,'stage_1'] <- rowMeans(cbind(props[[3]][,2], props[[4]][,1], props[[5]][,1], props[[6]][,1]))
props.mean[,'stage_2'] <- rowMeans(cbind(c(props[[2]][1:5,1],0,props[[2]][6:11,1]), props[[4]][,2], props[[5]][,2], props[[6]][,2]))
props.mean[,'stage_3'] <- rowMeans(cbind(props[[4]][,3], props[[5]][,3]))
props.mean[,'stage_LHT'] <- rowMeans(cbind(c(props[[1]][1:2,1],0,0,0,props[[1]][3:5,1],0,props[[1]][6:8,1]), c(props[[2]][1:5,2],0,props[[2]][6:11,2])))

# sd (exclude stage -1 which has no replicates)
props.sd <- matrix(ncol = 5, nrow = length(cols))
colnames(props.sd) <- paste0("stage_",c(0:3,"LHT"))
row.names(props.sd) <- names(cols)
props.sd[,'stage_0'] <- rowSds(cbind(props[[3]][,1], props[[7]][,2]))
props.sd[,'stage_1'] <- rowSds(cbind(props[[3]][,2], props[[4]][,1], props[[5]][,1], props[[6]][,1]))
props.sd[,'stage_2'] <- rowSds(cbind(c(props[[2]][1:5,1],0,props[[2]][6:11,1]), props[[4]][,2], props[[5]][,2], props[[6]][,2]))
props.sd[,'stage_3'] <- rowSds(cbind(props[[4]][,3], props[[5]][,3]))
props.sd[,'stage_LHT'] <- rowSds(cbind(c(props[[1]][1:2,1],0,0,0,props[[1]][3:5,1],0,props[[1]][6:8,1]), c(props[[2]][1:5,2],0,props[[2]][6:11,2])))

## add stage -1 without sd
props.mean <- cbind(props[["batch_7"]][,1], props.mean)
colnames(props.mean)[1] <- "stage_-1"
props.sd <- cbind(rep(0, nrow(props.sd)), props.sd)
colnames(props.sd)[1] <- "stage_-1"

# plot
df <- data.frame(pop = rep(row.names(props.mean), 6), 
                 stage = rep(substr(colnames(props.mean),7,9), each=nrow(props.mean)), 
                 mean = c(props.mean), sd = c(props.sd))

pdf(paste0(out, "FigS1_propsStage.pdf"), width = 8, height = 6, useDingbats = FALSE)
ggplot(df, aes(x=stage, y=mean, group=pop, color=pop)) + 
  geom_line() + 
  geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) + 
  scale_color_manual(values=cols) + 
  facet_wrap(. ~pop, ncol=4) + 
  ylab("% of cells in cluster") + 
  th + theme(legend.position = "none")
dev.off()


## Proportion of cells in G1 per stage (S1C)======
cell_cycle <- cbind(cell_cycle, umap[match(row.names(cell_cycle), row.names(umap)),])

props <- as.matrix(table(cell_cycle$V2, cell_cycle$cluster))
props <- t(t(props)/colSums(props))*100
props <- props[,order(props[1,])]

palette(c("#626D71", "#B38867", "#DDBC95"))

pdf(paste0(out, "FigS1_cellCycle_perCluster.pdf"), useDingbats = FALSE)
barplot(props, col=1:3, las=2, xlim=c(0,18))
legend("topright", legend = row.names(props), col=1:3, pch=15, cex=0.85)
dev.off()



### Figure S2 ###############

## QC metrics (S2A)======
qc <- read.table(paste0(dir, "data/QCstats_allCells.tsv"))
meta <- read.table(paste0(dir, "data/SupplementaryTable1.tab"), header = TRUE, row.names = 1)
stopifnot(identical(row.names(qc), row.names(meta)))

plots <- list()
plots[[1]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=log10(libSize+1), fill=as.factor(meta$batch), alpha=0.5)) + 
  geom_violin() + 
  geom_boxplot(width=0.05) + 
  scale_fill_manual(values = brewer.pal(n=8, "Set3")[-2]) + 
  ylab(expression('log'[10]*' library size')) + xlab("batch") + 
  geom_hline(yintercept = log10(50000), lwd=0.25, lty=2, col="grey") + 
  ggtitle("total reads in genes") + 
  th + 
  theme(legend.position="none")
plots[[2]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=nGenes/1e3, fill=as.factor(meta$batch), alpha=0.5)) + 
  geom_violin() + 
  geom_boxplot(width=0.05) + 
  scale_fill_manual(values = brewer.pal(n=8, "Set3")[-2]) + 
  ylab("total genes x1000") + xlab("batch") + 
  geom_hline(yintercept = 6, lwd=0.25, lty=2, col="grey") + 
  ggtitle("number of genes detected") + 
  th + 
  theme(legend.position="none")
plots[[3]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=mit/libSize*100, fill=as.factor(meta$batch), alpha=0.5)) + 
  geom_violin() + 
  geom_boxplot(width=0.05) + 
  scale_fill_manual(values = brewer.pal(n=8, "Set3")[-2]) + 
  ylab("% reads in MT genes") + xlab("batch") + 
  geom_hline(yintercept = 15, lwd=0.25, lty=2, col="grey") + 
  ggtitle("% reads in mitochondrial genes") + 
  th + 
  theme(legend.position="none")
plots[[4]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=ercc/libSize*100, fill=as.factor(meta$batch), alpha=0.5)) + 
  geom_violin(scale="width") + 
  geom_boxplot(width=0.05) + 
  scale_fill_manual(values = brewer.pal(n=8, "Set3")[-2]) + 
  ylab("% reads in spike-ins") + xlab("batch") + 
  geom_hline(yintercept = 30, lwd=0.25, lty=2, col="grey") + 
  ggtitle("% reads in ERCC spike-ins") + 
  th + 
  theme(legend.position="none")

pdf(paste0(out, "FigS2_QCstats.pdf"), width = 12, height = 3, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 4, nrow = 1)
dev.off()


## Batch correction (S2B)======
umap_bb <- read.table(paste0(dir, "results/umapCoords.tab"))
umap_bb$batch <- paste0("batch_", meta[row.names(umap_bb),]$batch)

order <- sample(row.names(umap_bb), size = nrow(umap_bb), replace = FALSE)
pdf(paste0(out, "FigS2_UMAPbefore.pdf"), useDingbats = FALSE)
ggplot(umap_bb[order,], aes(x, y, colour=as.factor(batch))) + 
  geom_point() + 
  scale_color_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + 
  labs(colour="batch") + 
  xlab("UMAP - dim 1") + ylab("UMAP - dim 2") + 
  th + 
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())
dev.off()

order <- sample(row.names(umap), size = nrow(umap), replace = FALSE)
pdf(paste0(out, "FigS2_UMAPafter.pdf"), useDingbats = FALSE)
ggplot(umap[order,], aes(x, y, colour=as.factor(batch))) + 
  geom_point() + 
  scale_color_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + 
  labs(colour="batch") + 
  xlab("UMAP - dim 1") + ylab("UMAP - dim 2") +
  th + 
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank())
dev.off()



### Figure S3 ###############

## Marker expression (S3A-D)======
stopifnot(identical(row.names(clusters), row.names(umap)))
plotGeneOnUMAP <- function(gene="Nkx2-5") {
  id <- row.names(rowData(sce)[rowData(sce)$gene==gene,])
  df <- data.frame(x = umap$x, y = umap$y, log2.exp = logcounts(sce)[id,row.names(umap)])
  df <- df[order(df$log2.exp),]
  
  plots <- list()
  plots[[1]] <- ggplot(df, aes(x,y)) + 
    geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + 
    scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + 
    ggtitle(gene) + xlab("") + ylab("") + 
    labs(colour=expression('log'[2]*' counts')) + 
    th + 
    theme(axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          legend.position = "bottom", 
          legend.text.align = 0.5, 
          legend.title.align = 0.5, 
          legend.box.margin=margin(-20,0,0,0)) + 
    guides(colour = guide_colorbar(title.position = "bottom"))
  
  df <- data.frame(log2.exp = logcounts(sce)[id,row.names(umap)], cluster=clusters$ann)
  
  plots[[2]] <- ggplot(df, aes(cluster, log2.exp, fill=cluster)) + 
    geom_boxplot() + 
    scale_fill_manual(values = cols) + 
    xlab("") + ylab(expression('log'[2]*' counts')) + 
    ggtitle(gene) + 
    th + 
    theme(legend.position = "none", 
          axis.text.x = element_text(color = cols, face="bold", size=12, angle=45, vjust=0.5), 
          plot.margin = margin(0.2, 0.2, 1, 0.2, "cm"))
  ggarrange(plotlist = plots, ncol = 2, nrow = 1, widths = c(0.6,0.4))
}

pdf(paste0(out, "FigS3_Cdh1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Cdh1")
dev.off()
pdf(paste0(out, "FigS3_Sox2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Sox2")
dev.off()
pdf(paste0(out, "FigS3_Ttn.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Ttn")
dev.off()
pdf(paste0(out, "FigS3_Nkx2-5.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Nkx2-5")
dev.off()
pdf(paste0(out, "FigS3_Actn2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Actn2")
dev.off()
pdf(paste0(out, "FigS3_Sox17.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Sox17")
dev.off()


### Figure S4 ###############

pdf(paste0(out, "FigS4_Emcn.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Emcn")
dev.off()

pdf(paste0(out, "FigS4_Tfap2a.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tfap2a")
dev.off()

pdf(paste0(out, "FigS4_Isl1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Isl1")
dev.off()


### Figure S5 ###############

## UMAP with classes used for RF (S5A)======
umap.ref$label <- as.character(umap.ref$ann)
umap.ref[umap.ref$label == "dorsalMesoderm" & umap.ref[,2]>0,]$label <- "dorsalMesoderm_1"
umap.ref[umap.ref$label == "dorsalMesoderm" & umap.ref[,2]<0,]$label <- "dorsalMesoderm_2"
umap.ref[umap.ref$label == "caudalPSM" & umap.ref[,2]>0,]$label <- "caudalPSM_13"
umap.ref[umap.ref$label == "caudalPSM" & umap.ref[,2]<0,]$label <- "caudalPSM_2"
umap.ref[umap.ref$label == "cardiacMesoderm" & umap.ref$cluster == "cluster2",]$label <- "cardiacMesoderm_2"
umap.ref[umap.ref$label == "cardiacMesoderm" & umap.ref$cluster != "cluster2",]$label <- "cardiacMesoderm_56"
umap.ref$label <- as.factor(umap.ref$label)

pdf(paste0(out, "FigS5_UMAP_refData_labels.pdf"), width = 5, height = 5, useDingbats = FALSE)
plot(umap.ref$V2, umap.ref$V1, 
     axes=FALSE, 
     xlab="UMAP_1", ylab="UMAP_2", 
     col=umap.ref$label, pch=16)
box(bty="l")
legend("topleft", legend = levels(umap.ref$label), pch=16, col=1:7)
dev.off()

## Important genes (S5B)=====
forest <- readRDS(paste0(dir, "results/randomForest_referenceCells.Rds"))
imp.genes <- row.names(forest$importance[order(forest$importance[,7], decreasing = TRUE),][1:50,])

tmp <- sce.ref[,-which(sce.ref$label == "caudalPSM_2")]
tmp$label <- droplevels(tmp$label)

dat <- logcounts(tmp)[imp.genes,]
dat <- t(apply(dat, 1, function(x) (x-mean(x))/sd(x)))
row.names(dat) <- rowData(sce.ref)[row.names(dat),]$gene

labs <- factor(tmp$label, labels = 1:6)

pdf(paste0(out, "FigS5_impGenes_heatmap.pdf"), useDingbats = FALSE)
heatmap.2(dat, trace="none", col=rev(brewer.pal(n=10, "RdYlBu")), 
          ColSideColors = as.character(labs), 
          labCol = NA, 
          key.title = "", 
          key.xlab = "z-score", 
          density.info = "none")
dev.off()


## Prediction probs (S5C)=====
calls <- data.frame(class=labels, 
                    max.prob=apply(probs, 2, max), 
                    closest=apply(probs, 2, function(x) max(x[x!=max(x)]) ))
calls$diff <- calls$max.prob-calls$closest
calls$pass <- ifelse(calls$diff > 0.15, "yes", "no")

pdf(paste0(out, "FigS5_predictionProbs.pdf"), useDingbats = FALSE)
ggplot(calls, aes(x=class.V2, y=max.prob, colour=pass)) + 
  geom_boxplot() + 
  xlab("class") + ylab("max prob") + 
  th +
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()


### Figure S6 ###############

pdf(paste0(out, "FigS6_Tbx1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tbx1")
dev.off()

pdf(paste0(out, "FigS6_Asb2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Asb2")
dev.off()


### Figure S8 ###############

pdf(paste0(out, "FigS8_Mab21l2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Mab21l2")
dev.off()

pdf(paste0(out, "FigS8_Fst.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Fst")
dev.off()

pdf(paste0(out, "FigS8_Foxc2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Foxc2")
dev.off()


### Figure S9 ###############

## Branching points (S9A)======
pdf(paste0(out, "FigS9_diffMap_branches.pdf"), useDingbats = FALSE)
plot(dpt, root=1, paths_to=c(2,3), col_by="branch", dcs=c(2,-1)) + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank())
dev.off()

## Me5->Me3 trajectory (S9B)======
diffMap.me5 <- readRDS(paste0(dir, "results/diffusionMap_me5.Rds"))
diffMap.me5 <- as.data.frame(diffMap.me5@eigenvectors[,1:2])
diffMap.me5$ann <- sce[,row.names(diffMap.me5)]$clusterAnn
diffMap.me5$dpt <- diffMap[row.names(diffMap.me5),]$dpt

dpt.me5 <- readRDS(paste0(dir, "results/diffusionPseudotime_Me5.Rds"))

plots <- list()
plots[[1]] <- ggplot(diffMap.me5, aes(DC1, DC2, colour=ann)) + 
  geom_point() + 
  scale_color_manual(values=cols) + 
  labs(colour="pop") + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

plots[[2]] <- ggplot(diffMap.me5, aes(DC1, DC2, colour=dpt)) + 
  geom_point() + 
  scale_color_gradientn(colours = rev(brewer.pal(n=9, "YlOrBr")[-c(1)])) + 
  labs(colour="pop") + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

plots[[3]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,3)) + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

pdf(paste0(out, "FigS9_DMme5.pdf"), width = 10, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=3, nrow = 1, align = "h")
dev.off()


## Me7->Me3 trajectory (S9C)======
diffMap.me7 <- readRDS(paste0(dir, "results/diffusionMap_me7.Rds"))
diffMap.me7 <- as.data.frame(diffMap.me7@eigenvectors[,1:2])
diffMap.me7$ann <- sce[,row.names(diffMap.me7)]$clusterAnn
diffMap.me7$dpt <- diffMap[row.names(diffMap.me7),]$dpt

dpt.me7 <- readRDS(paste0(dir, "results/diffusionPseudotime_Me7.Rds"))

plots <- list()
plots[[1]] <- ggplot(diffMap.me7, aes(DC1, DC2, colour=ann)) + 
  geom_point() + 
  scale_color_manual(values=cols) + 
  labs(colour="pop") + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

plots[[2]] <- ggplot(diffMap.me7, aes(DC1, DC2, colour=dpt)) + 
  geom_point() + 
  scale_color_gradientn(colours = rev(brewer.pal(n=9, "YlOrBr")[-c(1)])) + 
  labs(colour="pop") + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

plots[[3]] <- plot(dpt.me7, col_by="branch", root=2, path=c(1,3)) + 
  th + theme(axis.ticks.y = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             legend.position = "bottom")

pdf(paste0(out, "FigS9_DMme7.pdf"), width = 10, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=3, nrow = 1, align="h")
dev.off()


## relationship with stage (S9D-E)======
diffMap$stage <- sce[,row.names(diffMap)]$stage

pdf(paste0(out, "FigS9_perStage.pdf"), width = 7, height = 4.5, useDingbats = FALSE)
ggplot(diffMap, aes(DC2, -DC1, colour=stage)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_color_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + 
  facet_wrap(~stage) + 
  th + 
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
dev.off()

plots <- list()
plots[[1]] <- ggplot(diffMap, aes(-DC1, fill=stage, alpha=0.5)) + 
  geom_density() + 
  ylim(c(0,80)) + 
  scale_fill_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + 
  th + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plots[[2]] <- ggplot(diffMap, aes(DC2, fill=stage, alpha=0.5)) + 
  geom_density() + 
  ylim(c(0,50)) + 
  scale_fill_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + 
  th + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plots[[3]] <- ggplot(diffMap, aes(-DC1, fill="black", alpha=0.5)) + 
  geom_density() + 
  ylim(c(0,80)) + 
  th + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plots[[4]] <- ggplot(diffMap, aes(DC2, fill="black", alpha=0.5)) + 
  geom_density() + 
  ylim(c(0,50)) + 
  th + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pdf(paste0(out, "FigS9_perStage_density.pdf"), width = 7, height = 7, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 2, nrow = 2)
dev.off()


### Figure S10 ###############

## Combination of marker genes used to identify clusters in staining experiments
sce.cardiac <- sce[,sce$clusterAnn %in% paste0("Me",3:8)]
sce.cardiac <- sce.cardiac[,-which(reducedDim(sce.cardiac)[,2]>5)] # remove outlier cell

## markers that differentiate between all cardiac clusters
m <- c("Fst", "Tbx1", "Vsnl1", "Fsd2", "Nkx2-5", "Mab21l2", "Asb2")
m <- row.names(rowData(sce.cardiac)[rowData(sce.cardiac)$gene %in% m, ])
df <- reducedDim(sce.cardiac)
df <- cbind(df, t(logcounts(sce.cardiac)[m,]))
colnames(df)[4:10] <- rowData(sce)[colnames(df)[4:10],]$gene

cols.shiny <- hcl.colors(n=6, palette = "Viridis")
cols.shiny <- c("#6D187C","grey",cols.shiny[-1])
plots <- list()
for(i in 4:10){
  plots[[i-3]] <- ggplot(df, aes(x, y)) + 
    # geom_point(colour="lightgrey", alpha = 0.5, size = 1) +
    stat_density_2d(data=df[df[,i]>1,], aes(x, y, alpha=..level..),
                    fill = cols.shiny[i-3], 
                    size = 2, bins = 10, geom = 'polygon') +
    xlim(c(-10, -2.5)) + ylim(c(-7, 4)) +
    ggtitle(colnames(df)[i]) + xlab("") + ylab("") +
    th + theme(axis.text.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               legend.position = "none",
               panel.background = element_rect(fill = "transparent"),
               plot.background = element_rect(fill = "transparent", color = NA))
}
plots[[8]] <- ggplot(df, aes(x, y)) + 
  geom_point(colour="lightgrey", alpha = 0.5, size = 1) +
  xlim(c(-10, -2.5)) + ylim(c(-7, 4)) +
  ggtitle("All") + xlab("") + ylab("") +
  th + theme(axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position = "none",
             panel.background = element_rect(fill = "transparent"),
             plot.background = element_rect(fill = "transparent", color = NA))

pdf(paste0(out, "FigS10_clusterStainingMarkers.pdf"), width = 12, height = 6, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=4, nrow=2)
dev.off()


### Figure S11 ###############

pdf(paste0(out, "FigS11_Smarcd3.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Smarcd3")
dev.off()

pdf(paste0(out, "FigS11_Hand1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Hand1")
dev.off()

pdf(paste0(out, "FigS11_Snai1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Snai1")
dev.off()

pdf(paste0(out, "FigS11_Tbx5.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tbx5")
dev.off()

pdf(paste0(out, "FigS11_Ttn.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Ttn")
dev.off()


### Figure S14 ###############

pdf(paste0(out, "FigS14_Fsd2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Fsd2")
dev.off()

pdf(paste0(out, "FigS14_Vsnl1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Vsnl1")
dev.off()


### Figure S19 ###############

pdf(paste0(out, "FigS19_Tbx18.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tbx18")
dev.off()
