###### Figures ######
library(scran)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(ggpubr)
library(locfit)

dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/mouseHeart_earlyDev_atlas/"
out <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/WRITING/HEART/Figures/figureElements/"

th <- theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size=10), axis.title.x = element_text(size=12), axis.text.y = element_text(size=10), axis.title.y = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5))



### Figure 1

## UMAP with clusters
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
plot(umap$x[order], umap$y[order], pch=16, col=clusters$col[order], bty="l", xlab="UMAP dim1", ylab="UMAP dim2")
legend("bottomright", legend = ann[order(ann)], col=cols[ann[order(ann)]], pch=16, cex=0.5)
dev.off()


## Expression of marker genes
sce <- readRDS(paste0(dir, "data/sce_goodQual.NORM.clusters.Rds"))

# markers <- c("Ttr","Hhex","Sox17","Foxa2","Pax9","Cdh1","Sox2","Dlx5","Tfap2a","Wnt6","Sox1","Snai1","Nkx2-5","Smarcd3","Ttn","Myh6","Actn2","Gata1","Runx1","Hba-a1","Tal1","Pecam1","Tie1","Cdh5","Emcn")
markers <- c("Ttr","Hhex","Sox17","Pax9","Cdh1","Sox2","Dlx5","Wnt6","Sox1","Snai1","Isl1","Smarcd3","Nkx2-5","Mef2c","Gata4","Ttn","Myh6","Actn2","Gata1","Hba-a1","Tal1","Cdh5","Emcn")
tmp <- logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene %in% markers,]),]
tmp <- t(apply(tmp, 1, function(x) x/max(x)))
row.names(tmp) <- rowData(sce)[row.names(tmp),]$gene
tmp <- tmp[match(markers, row.names(tmp)),]

palette <- brewer.pal(n=9, name="Purples")[-9]

clust <- hclust(dist(t(tmp)))
order <- colnames(tmp)[clust$order]
order <- clusters[order,]
order <- order[order(order$ann),]
order <- rbind(order[order$ann=="En2",], order[order$ann=="En1",], order[order$ann=="Ec1",], order[order$ann=="Ec2",],
           order[order$ann=="Me8",], order[order$ann=="Me7",], order[order$ann=="Me5",], order[order$ann=="Me6",],
           order[order$ann=="Me4",], order[order$ann=="Me3",], order[order$ann=="Me1",], order[order$ann=="Me2",])

ha = HeatmapAnnotation(df = data.frame(cluster = order$ann), col = list(cluster = cols))

pdf(paste0(out, "Fig1_heatmap_markers.pdf"), width = 7, height = 7, useDingbats = FALSE)
Heatmap(tmp[,row.names(order)], col=palette, top_annotation = ha, cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE)
dev.off()


## Cell-cycle phase assignment
cell_cycle <- read.table(paste0(dir, "results/cellCyclePhase.tsv"), stringsAsFactors = FALSE, row.names = 1)
stopifnot(identical(row.names(umap), row.names(cell_cycle)))

palette(c("#626D71", "#B38867", "#DDBC95"))

order <- sample(1:nrow(umap), nrow(umap), replace = FALSE)
pdf(paste0(out, "Fig1_cellCycle.pdf"), width = 7, height = 6, useDingbats = FALSE)
plot(umap$x[order], umap$y[order], pch=16, cex=0.75, col=as.factor(cell_cycle$V2[order]), bty="l", xlab="", ylab="", axes = FALSE)
box(bty="l")
legend("topright", legend = levels(as.factor(cell_cycle$V2[order])), col=1:3, pch=16, cex=0.75)
dev.off()


## Expression of FHF and SHF markers
heartFields <- c("Tbx5", "Hand1", "Sfrp5", "Hcn4", "Tbx1", "Foxc2", "Hoxb1", "Isl1")

tmp <- logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene %in% heartFields,]), sce$clusterAnn %in% paste0("Me", 3:8)]
row.names(tmp) <- rowData(sce)[row.names(tmp),]$gene
pops <- sce[,sce$clusterAnn %in% paste0("Me", 3:8)]$clusterAnn
pops <- factor(pops, levels = c(paste0("Me",5:3), paste0("Me",6:8)))

plots <- list()
for(gene in row.names(tmp)){
  df <- data.frame(expr=tmp[gene,], cluster=pops)
  plots[[gene]] <- ggplot(df, aes(pops, expr, col=pops, alpha=0.95)) + geom_violin(scale="width") + ylim(c(0,10)) + scale_colour_manual(values=cols) + ggtitle(gene) + th + theme(legend.position = "none")
}

pdf(paste0(out, "Fig1_FSHF_markers.pdf"), width = 12, height = 7, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 4, nrow = 2)
dev.off()




### Figure 2

## Diffusion map
diffMap <- readRDS(paste0(dir, "results/diffusionMap_cardiacMesoderm.Rds"))
diffMap <- as.data.frame(diffMap@eigenvectors[,1:5])
diffMap$cluster <- clusters[row.names(diffMap),]$ann

pdf(paste0(out, "Fig2_diffMap.pdf"), useDingbats = FALSE)
ggplot(diffMap, aes(DC2, DC1, colour=cluster)) + geom_point() + scale_color_manual(values = cols[7:12]) + th
dev.off()

## Diffusion pseudotime
dpt <- readRDS(paste0(dir, "results/diffusionPseudotime_cardiacMesoderm.Rds"))
# which(dpt@tips[,1])
diffMap$dpt <- dpt$DPT770

pdf(paste0(out, "Fig2_diffPDT.pdf"), useDingbats = FALSE)
ggplot(diffMap, aes(DC2, DC1, colour=dpt)) + geom_point() + scale_color_gradientn(colours = brewer.pal(n=9, "YlOrBr")[-c(1)]) + th
dev.off()

## Add branches
pdf(paste0(out, "Fig2_diffMap_paths.pdf"), useDingbats = FALSE)
plot(dpt, root=2, paths_to=c(1,3), col_by="branch",dcs=c(2,1)) + th
dev.off()

## Nkx2-5 expr
df <- data.frame(x = diffMap$DC2, y = diffMap$DC1, log2.exp = logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene == "Nkx2-5",]),row.names(diffMap)])
df <- df[order(df$log2.exp),]

pdf(paste0(out, "Fig2_Nkx2-5_diffMap.pdf"), useDingbats = FALSE)
ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle("Nkx2-5") + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))
dev.off()

## Diffusion pseudotime for Me5->M3
dpt <- read.table(paste0(dir, "results/diffusionPseudotime_Me5.tsv"), row.names = 1)

df$me5 <- dpt[match(row.names(df), row.names(dpt)),1]
fit2 <- locfit(df[!is.na(df$me5),3]~lp(df[!is.na(df$me5),4], nn=1, deg=2), data=df)
curve <- predict(fit2, df[!is.na(df$me5),4])

pdf(paste0(out, "Fig2_Nkx2-5_Me5-3.pdf"), useDingbats = FALSE)
ggplot(df[!is.na(df$me5),], aes(me5, log2.exp)) + geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me5),]),]$col)) + ylim(c(0,max(df$me5))) + ggtitle("Nkx2-5") + xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + th + theme(axis.text.x = element_blank()) + geom_line(aes(df[!is.na(df$me5),]$me5, curve), lwd=1)
dev.off()




### Figure 3

## Mab21l2 expression in diffusion map
df <- data.frame(x = diffMap$DC1, y = diffMap$DC2, log2.exp = logcounts(sce)[row.names(rowData(sce)[rowData(sce)$gene == "Mab21l2",]),row.names(diffMap)])
df <- df[order(df$log2.exp),]

pdf(paste0(out, "Fig3_Mab21l2_diffMap.pdf"), useDingbats = FALSE)
ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle("Mab21l2") + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))
dev.off()


## Mab21l2 expressin in Me5->Me3 trajectory
df$me5 <- dpt[match(row.names(df), row.names(dpt)),1]
fit2 <- locfit(df[!is.na(df$me5),3]~lp(df[!is.na(df$me5),4], nn=1, deg=2), data=df)
curve <- predict(fit2, df[!is.na(df$me5),4])

pdf(paste0(out, "Fig3_Mab21l2_Me5-3.pdf"), useDingbats = FALSE)
ggplot(df[!is.na(df$me5),], aes(me5, log2.exp)) + geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me5),]),]$col)) + ylim(c(0,max(df$me5))) + ggtitle("Mab21l2") + xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + th + theme(axis.text.x = element_blank()) + geom_line(aes(df[!is.na(df$me5),]$me5, curve), lwd=1)
dev.off()




##########

### Figure S1

## UMAP per stage
umap$stage <- colData(sce)$stage
umap$cluster <- sce$clusterAnn
umap$col <- sce$clusterCol

pdf(paste0(out, "FigS1_UMAPstages.pdf"), width = 10, height = 4.5, useDingbats = FALSE)
stages <- c(-1:3,"LHT")
par(mfrow=c(2,3), mar=c(2,2,2,2))
for(stage in stages){
  plot(umap[umap$stage==stage,]$x, umap[umap$stage==stage,]$y, col=umap[umap$stage==stage,]$col, pch=16, ylim=c(-6.5,11), xlim=c(-10,17), main=paste("stage",stage), xlab="", ylab="", axes=FALSE)
  box(bty="l")
}
dev.off()

## Proportion of cells per cluster across stages

# use the different batches as replicates, to get estimates with mean +- se
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

# plot
df <- data.frame(pop = rep(row.names(props.mean),5), stage = rep(substr(colnames(props.mean),7,9), each=nrow(props.mean)), mean = c(props.mean), sd = c(props.sd))
pdf(paste0(out, "FigS1_propsStage.pdf"), width = 8, height = 6, useDingbats = FALSE)
ggplot(df, aes(x=stage, y=mean, group=pop, color=pop)) + geom_line() + geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) + scale_color_manual(values=cols) + facet_wrap(. ~ pop, ncol=4) + ylab("% of cells in cluster") + th + theme(legend.position = "none")
dev.off()


## Proportion of cells in G1 per stage
cell_cycle <- cbind(cell_cycle, umap[match(row.names(cell_cycle), row.names(umap)),])

props <- as.matrix(table(cell_cycle$V2, cell_cycle$cluster))
props <- t(t(props)/colSums(props))*100
props <- props[,order(props[1,])]

palette(c("#626D71", "#B38867", "#DDBC95"))

pdf(paste0(out, "FigS1_cellCycle_perCluster.pdf"), useDingbats = FALSE)
barplot(props, col=1:3, las=2, xlim=c(0,18))
legend("topright", legend = row.names(props), col=1:3, pch=15, cex=0.85)
dev.off()




### Figure S2

## QC metrics
qc <- read.table(paste0(dir, "data/QCstats_allCells.tsv"))
meta <- read.table(paste0(dir, "data/SupplementaryTable1.tab"), header = TRUE, row.names = 1)
stopifnot(identical(row.names(qc), row.names(meta)))

plots <- list()
plots[[1]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=log10(libSize+1), fill=as.factor(meta$batch), alpha=0.5)) + geom_violin() + geom_boxplot(width=0.05) + scale_fill_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + ylab(expression('log'[10]*' library size')) + xlab("batch") + geom_hline(yintercept = log10(50000), lwd=0.25, lty=2, col="grey") + ggtitle("total reads in genes") + th + theme(legend.position="none")
plots[[2]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=nGenes/1e3, fill=as.factor(meta$batch), alpha=0.5)) + geom_violin() + geom_boxplot(width=0.05) + scale_fill_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + ylab("total genes x1000") + xlab("batch") + geom_hline(yintercept = 6, lwd=0.25, lty=2, col="grey") + ggtitle("number of genes detected") + th + theme(legend.position="none")
plots[[3]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=mit/libSize*100, fill=as.factor(meta$batch), alpha=0.5)) + geom_violin() + geom_boxplot(width=0.05) + scale_fill_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + ylab("% reads in MT genes") + xlab("batch") + geom_hline(yintercept = 15, lwd=0.25, lty=2, col="grey") + ggtitle("% reads in mitochondrial genes") + th + theme(legend.position="none")
plots[[4]] <- ggplot(qc, aes(x=as.factor(meta$batch), y=ercc/libSize*100, fill=as.factor(meta$batch), alpha=0.5)) + geom_violin(scale="width") + geom_boxplot(width=0.05) + scale_fill_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + ylab("% reads in spike-ins") + xlab("batch") + geom_hline(yintercept = 30, lwd=0.25, lty=2, col="grey") + ggtitle("% reads in ERCC spike-ins") + th + theme(legend.position="none")

pdf(paste0(out, "FigS2_QCstats.pdf"), width = 12, height = 3, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 4, nrow = 1)
dev.off()

## Batch correction
umap_bb <- read.table(paste0(dir, "results/umapCoords.tab"))
umap_bb$batch <- paste0("batch_", meta[row.names(umap_bb),]$batch)

order <- sample(row.names(umap_bb), size = nrow(umap_bb), replace = FALSE)
pdf(paste0(out, "FigS2_UMAPbefore.pdf"), useDingbats = FALSE)
ggplot(umap_bb[order,], aes(x, y, colour=as.factor(batch))) + geom_point() + scale_color_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + labs(colour="batch") + xlab("UMAP - dim 1") + ylab("UMAP - dim 2") + th + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()

order <- sample(row.names(umap), size = nrow(umap), replace = FALSE)
pdf(paste0(out, "FigS2_UMAPafter.pdf"), useDingbats = FALSE)
ggplot(umap[order,], aes(x, y, colour=as.factor(batch))) + geom_point() + scale_color_manual(values = palette(brewer.pal(n=8, "Set3")[-2])) + labs(colour="batch") + xlab("UMAP - dim 1") + ylab("UMAP - dim 2") + th + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
dev.off()




### Figure S3

## Marker expression
stopifnot(identical(row.names(clusters), row.names(umap)))
plotGeneOnUMAP <- function(gene="Nkx2-5") {
  id <- row.names(rowData(sce)[rowData(sce)$gene==gene,])
  df <- data.frame(x = umap$x, y = umap$y, log2.exp = logcounts(sce)[id,row.names(umap)])
  df <- df[order(df$log2.exp),]
  
  plots <- list()
  plots[[1]] <- ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(gene) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))
  
  df <- data.frame(log2.exp = logcounts(sce)[id,row.names(umap)], cluster=clusters$ann)
  
  plots[[2]] <- ggplot(df, aes(cluster, log2.exp, fill=cluster)) + geom_boxplot() + scale_fill_manual(values = cols) + xlab("") + ylab(expression('log'[2]*' counts')) + ggtitle(gene) + th + theme(legend.position = "none", axis.text.x = element_text(color = cols, face="bold", size=12, angle=45, vjust=0.5), plot.margin = margin(0.2, 0.2, 1, 0.2, "cm"))
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




### Figure S4

pdf(paste0(out, "FigS4_Emcn.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Emcn")
dev.off()

pdf(paste0(out, "FigS4_Tfap2a.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tfap2a")
dev.off()

pdf(paste0(out, "FigS4_Isl1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Isl1")
dev.off()




### Figure S5

pdf(paste0(out, "FigS5_Tbx1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tbx1")
dev.off()

pdf(paste0(out, "FigS5_Asb2.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Asb2")
dev.off()

pdf(paste0(out, "FigS5_Nkx2-5.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Nkx2-5")
dev.off()




### Figure S6

pdf(paste0(out, "FigS6_Smarcd3.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Smarcd3")
dev.off()

pdf(paste0(out, "FigS6_Hand1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Hand1")
dev.off()

pdf(paste0(out, "FigS6_Snai1.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Snai1")
dev.off()

pdf(paste0(out, "FigS6_Tbx5.pdf"), width = 7, height = 4, useDingbats = FALSE)
plotGeneOnUMAP(gene = "Tbx5")
dev.off()




### Figure S7

## Me5->Me3 trajectory
diffMap.me5 <- readRDS(paste0(dir, "results/diffusionMap_me5.Rds"))
dpt.me5 <- readRDS(paste0(dir, "results/diffusionPseudotime_Me5.Rds"))
sce.Me5 <- sce[,sce$clusterAnn %in% paste0("Me", c(3:6))]

plots <- list()
plots[[1]] <- ggplot(diffMap.me5, aes(DC1, DC2, colour=sce.Me5$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th
plots[[2]] <- plot(dpt.me5, col_by="branch", root=2, path=c(1,3)) + th

pdf(paste0(out, "FigSX_DMme5.pdf"), width = 10, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=2, nrow = 1)
dev.off()

pdf(paste0(out, "FigSX_DMme5_dpt.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(diffMap.me5, aes(DC1, DC2, colour=dpt.me5$DPT1135)) + geom_point() + scale_color_gradientn(colours = brewer.pal(n=9, "YlOrBr")[-c(1)]) + th
dev.off()

## Me7->Me3 trajectory
diffMap.me7 <- readRDS(paste0(dir, "results/diffusionMap_me7.Rds"))
dpt.me7 <- readRDS(paste0(dir, "results/diffusionPseudotime_Me7.Rds"))
sce.Me7 <- sce[,sce$clusterAnn %in% paste0("Me", c(3:4,6:8))]

plots <- list()
plots[[1]] <- ggplot(diffMap.me7, aes(DC1, DC2, colour=sce.Me7$clusterAnn)) + geom_point() + scale_color_manual(values=cols) + labs(colour="pop") + th
plots[[2]] <- plot(dpt.me7, col_by="branch", root=2, path=c(1,3)) + th

pdf(paste0(out, "FigSX_DMme7.pdf"), width = 10, height = 4, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol=2, nrow = 1)
dev.off()

pdf(paste0(out, "FigSX_DMme7_dpt.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(diffMap.me7, aes(DC1, DC2, colour=dpt.me7$DPT1161)) + geom_point() + scale_color_gradientn(colours = brewer.pal(n=9, "YlOrBr")[-c(1)]) + th
dev.off()

## relationship with stage
diffMap$stage <- sce[,row.names(diffMap)]$stage

pdf(paste0(out, "FigSX_perStage.pdf"), width = 7, height = 4.5, useDingbats = FALSE)
ggplot(diffMap, aes(DC2, DC1, colour=stage)) + geom_point(size=1, alpha=0.5) + scale_color_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + facet_wrap(~stage) + th + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
dev.off()

plots <- list()
plots[[1]] <- ggplot(diffMap, aes(DC1, colour=stage)) + geom_density() + ylim(c(0,70)) + scale_color_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
plots[[2]] <- ggplot(diffMap, aes(DC2, colour=stage)) + geom_density() + ylim(c(0,120)) + scale_color_manual(values = brewer.pal(n=8, "Blues")[-c(1:2)]) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
plots[[3]] <- ggplot(diffMap, aes(DC1)) + geom_density() + ylim(c(0,70)) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())
plots[[4]] <- ggplot(diffMap, aes(DC2)) + geom_density() + ylim(c(0,120)) + th + theme(axis.text = element_blank(), axis.ticks = element_blank())

pdf(paste0(out, "FigSX_perStage_density.pdf"), width = 7, height = 7, useDingbats = FALSE)
ggarrange(plotlist = plots, ncol = 2, nrow = 2)
dev.off()









