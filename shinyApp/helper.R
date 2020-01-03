# helper.R
load("data/data.RData")

selectList <- as.character(ann$gene)
selectListDM <- as.character(ann[union(names(dynGenes.me5.curves), names(dynGenes.me7.curves)),]$gene)

cols <- c(Ec1 = "#ec6646", Ec2 = "#af4424", En1 = "#3c537c", En2 = "#768ba5",
          Me1 = "#bf9a77", Me2 = "#debc95", Me3 = "#556dad", Me4 = "#f28a31", 
          Me5 = "#729f3c", Me6 = "#fbba14", Me7 = "#5fa398", Me8 = "#9FD3C5")

th <- theme_bw() + theme(axis.text.x = element_text(size=12), axis.title.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.y = element_text(size=12), axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), plot.title = element_text(face="bold", hjust = 0.5))


## geneExpression
plotGeneOnUMAP <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x = clusters$x, y = clusters$y, log2.exp = countsNorm[id,row.names(clusters)])
  df <- df[order(df$log2.exp),]
  
  plots <- list()
  plots[[1]] <- ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(gene) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))

  df <- data.frame(log2.exp = countsNorm[id,row.names(clusters)], cluster=clusters$population)
  
  plots[[2]] <- ggplot(df, aes(cluster, log2.exp, fill=cluster)) + geom_boxplot() + scale_fill_manual(values = cols) + xlab("") + ylab(expression('log'[2]*' counts')) + ggtitle(gene) + th + theme(legend.position = "none", axis.text.x = element_text(color = cols, face="bold", size=12, angle=45, vjust=0.5), plot.margin = margin(0.2, 0.2, 1, 0.2, "cm"))
  ggarrange(plotlist = plots, ncol = 2, nrow = 1, widths = c(0.6,0.4))
}


## geneExpression per-stage
plotGeneOnUMAPperStage <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x = clusters$x, y = clusters$y, log2.exp = countsNorm[id,row.names(clusters)] )
  df$stage <- meta[match(row.names(df), meta$cellID),]$stage
  df <- df[order(df$log2.exp),]

  stages <- c(-1:3,"LHT")
  plots <- list()
  for(stage in stages){
    plots[[paste0("stage_",stage)]] <- ggplot(df[df$stage == stage,], aes(x,y)) + geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ylim(c(-6.5,11)) + xlim(c(-10,17)) + ggtitle(paste("stage",stage)) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))
  }  
  ggarrange(plotlist = plots, ncol = 6, nrow = 1, common.legend = TRUE, legend = "bottom")
}



## cluster markers
printClusterMarkerTable <- function(pop=NULL, type='all tests significant'){
  chosen <- as.numeric(pop)
  if(type == "all tests significant"){
    t <- markers.all[[chosen]]
  }else{
    t <- markers.some[[chosen]]
  }
  t <- data.frame(gene=as.character(ann[row.names(t),1]), fdr = t$FDR, fold_change=round(t$logFC,2), row.names = row.names(t))
  t[,2] <- format(t[,2], digits=4, width=5)
  if(sum(duplicated(t$gene))>0){
    t$gene <- as.character(t$gene)
    idx <- which(duplicated(t$gene))
    g <- t[idx,]$gene
    idx <- which(t$gene == g)
    j <- 1
    for(i in idx){
      t[i,]$gene <- paste(g, j, sep=".")
      j <- j+1
    }
  }
  row.names(t) <- t$gene
  t$gene <- NULL
  return(t)
}

retrieveClusterMarkers <- function(pop=NULL, type='all tests significant'){
  chosen <- as.numeric(pop)
  if(type == "all tests significant"){
    t <- markers.all[[chosen]]
  }else{
    t <- markers.some[[chosen]]
  }
  t <- data.frame(gene=as.character(ann[row.names(t),1]), fdr = t$FDR, fold_change=round(t$logFC,2), row.names = row.names(t))
  t[,2] <- format(t[,2], digits=4, width=5)
  return(t)
}

plotClusterMarkerGeneOnUMAP <- function(pop=NULL, sel=NULL, type='all tests significant') {
  chosen <- as.numeric(pop)
  if(type == "all tests significant"){
    t <- markers.all[[chosen]]
  }else{
    t <- markers.some[[chosen]]
  }
  id <- row.names(t[sel,])
  gene <- as.character(ann[id,1])
  df <- data.frame(x = clusters$x, y = clusters$y, log2.exp = countsNorm[id,row.names(clusters)])
  df <- df[order(df$log2.exp),]
  
  plots <- list()
  plots[[1]] <- ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(gene) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_line(colour="white"), axis.text.x = element_text(color = "white", face="bold", size=12, angle=45, vjust=0.5), axis.text.y = element_text(colour="white"), legend.position = "none") #, plot.margin = margin(0.2, 0.2, 0, 0.2, "cm"))
  
  df <- data.frame(log2.exp = countsNorm[id,row.names(clusters)], cluster=clusters$population)
  plots[[2]] <- ggplot(df, aes(cluster, log2.exp, fill=cluster)) + geom_boxplot() + scale_fill_manual(values = cols) + xlab("") + ylab(expression('log'[2]*' counts')) + th + theme(legend.position = "none", axis.text.x = element_text(color = cols, face="bold", size=12, angle=45, vjust=0.5))
  
  ggarrange(plotlist = plots, ncol = 1, nrow = 2)
}


## Differentiation trajectories
plotGeneOnDiffMap <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x=diffMap$DC2, y=diffMap$DC1, log2.exp = countsNorm[id, row.names(diffMap)], me5 = pseudotime.me5[row.names(diffMap)], me7 = pseudotime.me7[row.names(diffMap)])
  df <- df[order(df$log2.exp),]
  
  plots <- list()
  plots[[1]] <- ggplot(df, aes(x,y)) + geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + scale_colour_gradientn(colours = colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)) + ggtitle(gene) + xlab("") + ylab("") + labs(colour=expression('log'[2]*' counts')) + th + theme(axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", legend.text.align = 0.5, legend.title.align = 0.5, legend.box.margin=margin(-20,0,0,0)) + guides(colour = guide_colorbar(title.position = "bottom"))
  plots[[2]] <- ggplot(df[!is.na(df$me5),], aes(me5, log2.exp)) + geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me5),]),]$col), alpha=0.5) + ylim(c(0,max(df$me5))) + ggtitle("Me5 --> Me3 trajectory") + xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + th + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 2, 0.2, "cm"))
  plots[[3]] <- ggplot(df[!is.na(df$me7),], aes(me7, log2.exp)) + geom_point(colour=as.character(clusters[row.names(df[!is.na(df$me7),]),]$col), alpha=0.5) + ylim(c(0,max(df$me7))) + ggtitle("Me7 --> Me3 trajectory") + xlab("pseudotime") + ylab(expression('log'[2]*' counts')) + th + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 2, 0.2, "cm"))
  
  if(id %in% names(dynGenes.me5.curves)){ plots[[2]] <- plots[[2]] + geom_line(aes(pseudotime.me5, dynGenes.me5.curves[[id]]), lwd=1) }
  if(id %in% names(dynGenes.me7.curves)){ plots[[3]] <- plots[[3]] + geom_line(aes(pseudotime.me7, dynGenes.me7.curves[[id]]), lwd=1) }
  
  ggarrange(plotlist = plots, ncol = 3, nrow = 1)
}

retrieveDynamicGenes <- function(){
  df <- rbind(ann[names(dynGenes.me5.curves),], ann[names(dynGenes.me7.curves),])
  df$trajectory <- c(rep("Me5-Me3", length(dynGenes.me5.curves)), rep("Me7-Me3",length(dynGenes.me7.curves)))
  return(df)
}


###############
#### data #####
###############
# dir <- "/Users/ibarra01/OneDrive - CRUK Cambridge Institute/github/mouseHeart_earlyDev_atlas/"

## normalised counts
# countsNorm <- readRDS(paste0(dir, "data/heartData_unbiased.goodQual.NORM.Rds"))

## gene info
# ann <- read.table(paste0(dir, "data/Mus_musculus.GRCm38.87.tsv"), stringsAsFactors = FALSE, sep="\t", header = TRUE, row.names = 1)
# colnames(ann) <- c("gene", "chr", "start", "end", "strand")
# ann <- ann[row.names(countsNorm),]

## metadata
# meta <- read.table(paste0(dir, "data/SupplementaryTable1.tab"), header = TRUE)
# meta <- meta[meta$cellID %in% colnames(countsNorm),]

## umap coords
# umap <- read.table(paste0(dir, "results/umapCoords_corrected.tab"))

## clusters
# clusters <- read.table(paste0(dir, "results/clusters_average_min40.tsv"), row.names = 1)
# remove outlier
# outlier <- which(clusters$V2==0)
# umap <- umap[-outlier,]
# clusters <- clusters[-outlier,]

## cluster annotation
# clust.ann <- c("En1","En2","Ec1","Ec2",paste0("Me",1:8)) # En, Ec and Me for the three germ layers
# names(clust.ann) <- c(3, 5, 7, 9, 12, 11, 1, 6, 4, 10, 2, 8) ## map the cluster number to the annotation
# clust.ann <- clust.ann[order(as.numeric(names(clust.ann)))]

## consolidate UMAP and clusters in one dataframe
# clusters <- data.frame(x = umap$x, y = umap$y, cluster = clusters, population = clust.ann[clusters], row.names = row.names(umap), stringsAsFactors = FALSE)
# clusters$col <- cols[clusters$population]

## markers
# markers.all <- readRDS(paste0(dir, "results/markerGenes_pval_all.Rds"))
# names(markers.all) <- clust.ann
# for(m in names(markers.all)){
#   markers.all[[m]] <- markers.all[[m]][markers.all[[m]]$FDR < 0.05, -1]
#   markers.all[[m]]$logFC <- apply(markers.all[[m]], 1, function(x) x[which.max(x)])
#   markers.all[[m]] <- as.data.frame(markers.all[[m]][,c(1,ncol(markers.all[[m]]))])
# }
# 
# markers.all <- list()
# for(m in clust.ann){
#   markers.all[[m]] <- read.table(paste0(dir, "shinyApp/data/markers_all_",m,".tsv"), stringsAsFactors = FALSE, header = TRUE)
# }

# markers.some <- readRDS(paste0(dir, "results/markerGenes_pval_some0.75.Rds"))
# names(markers.some) <- clust.ann
# for(m in names(markers.some)){
#   markers.some[[m]] <- markers.some[[m]][markers.some[[m]]$FDR < 0.05, -1]
#   markers.some[[m]]$logFC <- apply(markers.some[[m]], 1, function(x) x[which.max(x)])
#   markers.some[[m]] <- as.data.frame(markers.some[[m]][,c(1,ncol(markers.some[[m]]))])
# }
# markers.some <- list()
# for(m in clust.ann){
#   markers.some[[m]] <- read.table(paste0(dir, "shinyApp/data/markers_some_",m,".tsv"), stringsAsFactors = FALSE, header = TRUE)
# }

## trajectories
# diffMap <- readRDS(paste0(dir, "results/diffusionMap_cardiacMesoderm.Rds"))
# diffMap <- as.data.frame(diffMap@eigenvectors[,1:5])
# diffMap$col <- clusters[row.names(diffMap),]$col

# x <- read.table(paste0(dir, "results/diffusionPseudotime_Me5.tsv"), stringsAsFactors = FALSE, row.names = 1)
# pseudotime.me5 <- x$V2
# names(pseudotime.me5) <- row.names(x)
# x <- read.table(paste0(dir, "results/diffusionPseudotime_Me7.tsv"), stringsAsFactors = FALSE, row.names = 1)
# pseudotime.me7 <- x$V2
# names(pseudotime.me7) <- row.names(x)

# x <- read.table(paste0(dir, "results/deltaAIC_Me5.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me5 <- x$V2
# names(dynGenes.me5) <- row.names(x)
# dynGenes.me5 <- dynGenes.me5[dynGenes.me5 < -300]
# x <- read.table(paste0(dir, "results/deltaAIC_Me7.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me7 <- x$V2
# names(dynGenes.me7) <- row.names(x)
# dynGenes.me7 <- dynGenes.me7[dynGenes.me7 < -300]
# 
# library(locfit)
# dynGenes.me5.curves <- list()
# for(gene in names(dynGenes.me5)){
#   df <- data.frame(dpt=pseudotime.me5, log2.expr=countsNorm[gene,names(pseudotime.me5)])
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=1, deg=2), data=df)
#   dynGenes.me5.curves[[gene]] <- predict(fit2, df[,1])
# }
# dynGenes.me7.curves <- list()
# for(gene in names(dynGenes.me7)){
#   df <- data.frame(dpt=pseudotime.me7, log2.expr=countsNorm[gene,names(pseudotime.me7)])
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=1, deg=2), data=df)
#   dynGenes.me7.curves[[gene]] <- predict(fit2, df[,1])
# }

# save(countsNorm, ann, meta, clusters, markers.all, markers.some, diffMap, pseudotime.me5, pseudotime.me7, dynGenes.me5.curves, dynGenes.me7.curves, file=paste0(dir, "shinyApp/data/data.RData"))
###############
