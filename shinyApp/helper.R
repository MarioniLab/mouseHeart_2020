# helper functions

library(ggplot2)
library(ggpubr)
library(grid)

load("data/data.RData")

selectList <- as.character(ann$gene)
selectListDM <- as.character(ann$gene)
selectList_ref <- as.character(ann_ref$gene)

cols <- c(Ec1 = "#ec6646", Ec2 = "#af4424", En1 = "#3c537c", En2 = "#768ba5",
          Me1 = "#bf9a77", Me2 = "#debc95", Me3 = "#556dad", Me4 = "#f28a31", 
          Me5 = "#729f3c", Me6 = "#fbba14", Me7 = "#5fa398", Me8 = "#9FD3C5")
cols.ref <- c(CrM = "#75B3E2", VM = "#ED6B58", DM = "#A57CB5", PSM = "#589E46")
cols.ref.cluster <- c(cluster1 = "#25958C", cluster2 = "#6BB288", cluster3 = "#F4CC71", 
                      cluster4 = "#DF7A24", cluster5 = "#9C5016")

th <- theme_bw() + theme(axis.text.x = element_text(size=12), 
                         axis.title.x = element_text(size=12), 
                         axis.text.y = element_text(size=12), 
                         axis.title.y = element_text(size=12), 
                         axis.ticks.x = element_blank(), 
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), 
                         axis.line = element_line(colour = "black"), 
                         panel.border = element_blank(), 
                         plot.title = element_text(face="bold", hjust = 0.5))

axes <- list(
  zeroline = FALSE,
  showline = TRUE,
  showticklabels = FALSE,
  showgrid = FALSE
)

blue_pal <- colorRampPalette(c("grey", "lightblue", "dodgerblue4", "royalblue4"))(100)

## geneExpression
plotlyGeneOnUMAP <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x = clusters$x, y = clusters$y, pop = clusters$population, log2.exp = countsNorm[id,row.names(clusters)])
  df <- df[order(df$log2.exp),]
  
  lab <- "UMAP "
  plot_ly(data = df,
          type = "scatter", x = ~x, y = ~y,
          color = ~log2.exp, colors = blue_pal,
          mode = "markers",
          marker = list(size = 4, opacity = 0.5),
          width = 425, height = 300,
          hoverinfo = 'text',
          text = ~paste('log expr:', round(log2.exp, digits=1), '\ncluster: ', pop)
  ) %>%
    layout(
      autosize  = FALSE,
      title = gene,
      xaxis = c(list(title = paste0(lab,1)), axes),
      yaxis = c(list(title = paste0(lab,2)), axes),
      margin =  list(l=10, r=20, b=10, t=30)
    ) %>%
    colorbar(title = "log2 counts", x = 1, y = 0.75)
}

expr_boxplot <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x = clusters$x, y = clusters$y, pop = clusters$population, log2.exp = countsNorm[id,row.names(clusters)])
  df <- df[order(df$log2.exp),]
  
  plot_ly(data = df,
          type = "box",
          x = ~pop, y = ~log2.exp,
          color = ~pop, colors = cols,
          width = 350, height = 325
  )  %>%
    layout(
      autosize  = FALSE,
      title = gene,
      xaxis = c(list(title = "", tickangle = 90)),
      yaxis = c(list(title = "", showline=T, range=c(0, ceiling(max(df$log2.exp))))),
      margin = list(l=10, r=10, b=0, t=30),
      showlegend = FALSE
    )
}

numberPositiveCells <- function(gene="Nkx2-5"){
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(pop = clusters$population,
                   log2.exp = countsNorm[id,row.names(clusters)])
  df$log2.exp <- ifelse(df$log2.exp>0, 1, 0)
  
  df <- cbind(as.data.frame(table(df$pop, df$log2.exp)[,ncol(table(df$pop, df$log2.exp))]),
              as.data.frame(table(df$pop)))[,-2]
  colnames(df) <- c("number_cells", "percent_of_cluster")
  df$percent_of_cluster <- round(df$number_cells/df$percent_of_cluster*100, 1)
  return(t(df))
}


## geneExpression per-stage
plotUMAPperStage <- function() {
  df <- clusters[,c(1,2,4)]
  df$stage <- paste("stage", meta[match(row.names(df), meta$cellID),]$stage)

  ggplot(df, aes(x, y, colour=population)) + 
    geom_point(alpha = 0.5, cex=1.5) + 
    scale_color_manual(values=cols) +
    xlab("") + ylab("") + 
    facet_wrap(~stage, ncol=6) +
    labs(colour=expression('log'[2]*' counts')) + 
    th + theme(axis.ticks.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.text.y = element_blank(), 
               strip.text = element_text(size=12),
               legend.position = "bottom", 
               legend.text.align = 0.5, 
               legend.title.align = 0.5, 
               legend.box.margin=margin(-20,0,0,0)) + 
    guides(colour = guide_colorbar(title.position = "bottom"))
}

plotGeneOnUMAPperStage <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x = clusters$x, y = clusters$y, log2.exp = countsNorm[id,row.names(clusters)] )
  df$stage <- paste("stage", meta[match(row.names(df), meta$cellID),]$stage)
  df <- df[order(df$log2.exp),]

  ggplot(df, aes(x,y)) + 
    geom_point(aes(colour=log2.exp), alpha = 0.5, cex=1.5) + 
    scale_colour_gradientn(colours = blue_pal) + 
    xlab("") + ylab("") + 
    facet_wrap(~stage, ncol=6) +
    labs(colour=expression('log'[2]*' counts')) + 
    th + theme(axis.ticks.y = element_blank(), 
               axis.text.x = element_blank(), 
               axis.text.y = element_blank(), 
               strip.text = element_text(size=12),
               legend.position = "bottom", 
               legend.text.align = 0.5, 
               legend.title.align = 0.5, 
               legend.box.margin=margin(-20,0,0,0)) + 
    guides(colour = guide_colorbar(title.position = "bottom", barwidth = 8))
}


## cluster markers
printClusterMarkerTable <- function(pop=NULL, type='all tests significant'){
  chosen <- as.numeric(pop)
  if(type == "all tests significant"){
    t <- markers.all[[chosen]]
  }else{
    t <- markers.some[[chosen]]
  }
  t <- data.frame(gene=as.character(ann[row.names(t),1]), 
                  fdr = t$FDR, 
                  fold_change=round(t$logFC,2), 
                  row.names = row.names(t))
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
  t <- data.frame(gene=as.character(ann[row.names(t),1]), 
                  fdr = t$FDR, 
                  fold_change=round(t$logFC,2), 
                  row.names = row.names(t))
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
  df <- data.frame(x = clusters$x, y = clusters$y, 
                   pop = clusters$population,
                   log2.exp = countsNorm[id,row.names(clusters)])
  df <- df[order(df$log2.exp),]
  
  lab <- "UMAP "
  plot1 <- plot_ly(
    data = df,
    type = "scatter", x = ~x, y = ~y,
    color = ~log2.exp, colors = blue_pal,
    mode = "markers",
    marker = list(size = 4, opacity = 0.5),
    hoverinfo = 'text',
    text = ~paste('log expr:', round(log2.exp, digits=1), '\ncluster: ', pop)
    ) %>%
    layout(
      title = gene,
      xaxis = c(list(title = ""), axes),
      yaxis = c(list(title = ""), axes),
      margin =  list(l=10, r=20, b=10, t=30),
      annotations = list(
        text = gene,
        font = list(size = 16),
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      )
    )
  
  plot2 <- plot_ly(
    data = df,
    type = "box",
    x = ~pop, y = ~log2.exp,
    color = ~pop, colors = cols
  )  %>%
    layout(
      xaxis = c(list(title = "", tickangle = 90)),
      yaxis = c(list(title = "log2 counts", showline=T, range=c(0, ceiling(max(df$log2.exp))))),
      margin = list(l=20, r=10, b=0, t=30),
      showlegend = FALSE
    )
  
  subplot(plot1, plot2, nrows = 2, titleX = TRUE, titleY = TRUE, margin = 0.05, heights = c(0.6, 0.4)) %>% 
    layout(title = "", autosize = FALSE, width = 370, height = 500
    )  %>%
    colorbar(title = "log2 counts", x = 1, y = 1)
}


## Differentiation trajectories
plotGeneOnDiffMap <- function(gene="Nkx2-5") {
  id <- row.names(ann[ann$gene==gene,])
  df <- data.frame(x=diffMap$DC2, y=-diffMap$DC1, 
                   log2.exp = countsNorm[id, row.names(diffMap)], 
                   pop = clusters[row.names(diffMap),]$pop,
                   me5 = pseudotime.me5[row.names(diffMap)],
                   me7 = pseudotime.me7[row.names(diffMap)])
  if(gene %in% dynGenes$gene){
    if(dynGenes[dynGenes$gene == gene, ]$type %in% c("both", "Me5")){
      df$me5.curve <- curves_me5[id,match(row.names(df), colnames(curves_me5))]
    }
    if(dynGenes[dynGenes$gene == gene, ]$type %in% c("both", "Me7")){
      df$me7.curve <- curves_me7[id,match(row.names(df), colnames(curves_me7))]
    }
  }
  df <- df[order(df$log2.exp),]

  lab <- "DC "
  plot1 <- plot_ly(
    data = df,
    type = "scatter", x = ~x, y = ~y,
    color = ~log2.exp, colors = blue_pal,
    mode = "markers",
    marker = list(size = 4, opacity = 0.5),
    # width = 275, height = 275,
    hoverinfo = 'text',
    text = ~paste('log expr:', round(log2.exp, digits=1), '\ncluster: ', pop)
  ) %>%
    layout(
      # autosize  = FALSE,
      xaxis = c(list(title = paste0(lab,2)), axes),
      yaxis = c(list(title = paste0(lab,1)), axes),
      margin =  list(l=10, r=20, b=10, t=30),
      annotations = list(
        text = gene,
        font = list(size = 16),
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      )
    ) %>% 
    hide_colorbar()
  
  plot2 <- plot_ly(data = df[!is.na(df$me5),],
                  type = "scatter", x = ~me5, y = ~log2.exp,
                  color = ~pop, colors = cols,
                  mode = "markers",
                  marker = list(size = 4, opacity = 0.5),
                  # width = 250, height = 275,
                  hoverinfo = 'text',
                  text = ~paste('log expr:', round(log2.exp, digits=1), '\ncluster: ', pop)
  ) %>%
    layout(
      # autosize  = FALSE,
      title = "Me5 -> Me3",
      xaxis = c(list(title = "pseudotime"), axes),
      yaxis = c(list(title = "", showticklabels=TRUE, tickcolor = toRGB("black")), axes),
      margin =  list(l=10, r=20, b=10, t=30),
      annotations = list(
        text = "Me5 -> Me3",
        font = list(size = 16),
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      ),
      showlegend = FALSE
    )
  if(gene %in% dynGenes$gene){
    if(dynGenes[dynGenes$gene == gene, ]$type %in% c("both", "Me5")){
      plot2 <- plot2 %>% 
        add_trace(y = ~me5.curve, 
                  mode="markers", 
                  marker=list(color="black", size=3), hoverinfo='skip')
    }
  }
  
  plot3 <- plot_ly(data = df[!is.na(df$me7),],
                  type = "scatter", x = ~me7, y = ~log2.exp,
                  color = ~pop, colors = cols,
                  mode = "markers",
                  marker = list(size = 4, opacity = 0.5),
                  # width = 250, height = 275,
                  hoverinfo = 'text',
                  text = ~paste('log expr:', round(log2.exp, digits=1), '\ncluster: ', pop)
  ) %>%
    layout(
      # autosize  = FALSE,
      title = "Me7 -> Me3",
      xaxis = c(list(title = "pseudotime"), axes),
      yaxis = c(list(title = "", showticklabels=TRUE, tickcolor = toRGB("black")), axes),
      margin =  list(l=10, r=20, b=10, t=30),
      annotations = list(
        text = "Me7 -> Me3",
        font = list(size = 16),
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
      ),
      showlegend = FALSE
    )
  if(gene %in% dynGenes$gene){
    if(dynGenes[dynGenes$gene == gene, ]$type %in% c("both", "Me7")){
      plot3 <- plot3 %>% 
        add_trace(y = ~me7.curve, mode="markers", 
                  marker=list(color="black", size=3), 
                  hoverinfo='skip')
    }
  }
  
  subplot(plot1, plot2, plot3, titleX = TRUE, titleY = TRUE, margin = 0.03) %>% 
    layout(title = "", autosize = FALSE, width = 800, height = 275
  )
}

printDynamicGenesTable <- function(trajectory="both", cluster="all"){
  df <- dynGenes[,c(1,4:6)]
  row.names(df) <- df$gene
  df$gene <- NULL
  df$type <- as.factor(df$type)
  df$clusterMe5 <- as.factor(df$clusterMe5)
  df$clusterMe7 <- as.factor(df$clusterMe7)
  return(df)
}

plotGeneOnDiffMap_fromTable <- function(sel=1, trajectory="all", cluster="all") {
  df <- dynGenes[,c(1,4:6)]
  
  id <- row.names(df[sel,])
  gene <- df[sel,1]
  plotGeneOnDiffMap(gene = gene)
}

retrieveDynamicGenes <- function(){
  df <- dynGenes[,c(1,4:6)]
  return(df)
}


## geneExpression_reference
plotlyGeneOnUMAP_ref <- function(gene="Nkx2-5") {
  id <- row.names(ann_ref[ann_ref$gene==gene,])
  df <- umap_ref
  df$log2.exp <- countsNorm_ref[id, row.names(df)]
  df <- df[order(df$log2.exp),]
  
  lab <- "UMAP "
  plot_ly(data = df,
          type = "scatter", x = ~y, y = ~x,
          color = ~log2.exp, colors = blue_pal,
          mode = "markers",
          marker = list(size = 6, opacity = 0.5),
          width = 425, height = 300,
          hoverinfo = 'text',
          text = ~paste('log expr:', round(log2.exp, digits=1), '\nregion: ', dissection)
  ) %>%
    layout(
      autosize  = FALSE,
      title = gene,
      xaxis = c(list(title = paste0(lab,1)), axes),
      yaxis = c(list(title = paste0(lab,2)), axes),
      margin =  list(l=10, r=20, b=10, t=30)
    ) %>%
    colorbar(title = "log2 counts", x = 1, y = 0.75)
}

expr_boxplot_ref <- function(gene="Nkx2-5", col_by = "mesodermType") {
  id <- row.names(ann_ref[ann_ref$gene==gene,])
  df <- umap_ref
  df$log2.exp <- countsNorm_ref[id, row.names(df)]
  df <- df[order(df$log2.exp),]
  
  colors <- switch(col_by,
                   mesodermType = cols.ref,
                   cluster = cols.ref.cluster)
  column <- switch(col_by,
                   mesodermType = print(~dissection),
                   cluster = print(~cluster))
  
  plot_ly(data = df,
          type = "box",
          x = column, y = ~log2.exp,
          color = column, colors = colors,
          width = 350, height = 310
  )  %>%
    layout(
      autosize  = FALSE,
      title = gene,
      xaxis = c(list(title = "")),
      yaxis = c(list(title = "", showline=T, range=c(0, ceiling(max(df$log2.exp))))),
      margin = list(l=10, r=10, b=0, t=30),
      showlegend = FALSE
    )
}

numberPositiveCells_ref <- function(gene="Nkx2-5", col_by = "mesodermType"){
  id <- row.names(ann_ref[ann_ref$gene==gene,])
  df <- umap_ref
  df$log2.exp <- countsNorm_ref[id, row.names(df)]
  df$log2.exp <- ifelse(df$log2.exp > 0, 1, 0)
  
  if(col_by == "mesodermType"){
    df <- cbind(as.data.frame(table(df$dissection, df$log2.exp)[,ncol(table(df$dissection, df$log2.exp))]),
                as.data.frame(table(df$dissection)))[,-2]
    colnames(df) <- c("number_cells", "percent_of_cluster")
    df$percent_of_cluster <- round(df$number_cells/df$percent_of_cluster*100, 1)
  }else{
    df <- cbind(as.data.frame(table(df$cluster, df$log2.exp)[,ncol(table(df$dissection, df$log2.exp))]),
                as.data.frame(table(df$cluster)))[,-2]
    colnames(df) <- c("number_cells", "percent_of_cluster")
    df$percent_of_cluster <- round(df$number_cells/df$percent_of_cluster*100, 1)
  }
  return(t(df))
}

printImportantGenesTable <- function(){
  df <- imp.genes
  row.names(df) <- imp.genes$gene
  df$gene <- NULL
  df$MeanDecreaseAccuracy <- round(df$MeanDecreaseAccuracy, 4)
  df$MeanDecreaseGini <- round(df$MeanDecreaseGini, 4)
  return(df)
}

plotlyGeneOnUMAP_ref_fromTable <- function(sel=NULL) {
  id <- row.names(imp.genes[sel,])
  df <- umap_ref
  df$log2.exp <- countsNorm_ref[id, row.names(df)]
  df <- df[order(df$log2.exp),]
  
  lab <- "UMAP "
  plot_ly(data = df,
          type = "scatter", x = ~y, y = ~x,
          color = ~log2.exp, colors = blue_pal,
          mode = "markers",
          marker = list(size = 6, opacity = 0.5),
          width = 425, height = 300,
          hoverinfo = 'text',
          text = ~paste('log expr:', round(log2.exp, digits=1), '\nregion: ', dissection)
  ) %>%
    layout(
      autosize  = FALSE,
      title = ann_ref[id,1],
      xaxis = c(list(title = paste0(lab,1)), axes),
      yaxis = c(list(title = paste0(lab,2)), axes),
      margin =  list(l=10, r=20, b=10, t=30)
    ) %>%
    colorbar(title = "log2 counts", x = 1, y = 0.75)
}

expr_boxplot_ref_fromTable <- function(sel=1, col_by = "mesodermType") {
  id <- row.names(imp.genes[sel,])
  df <- umap_ref
  df$log2.exp <- countsNorm_ref[id, row.names(df)]
  df <- df[order(df$log2.exp),]
  
  colors <- switch(col_by,
                   mesodermType = cols.ref,
                   cluster = cols.ref.cluster)
  column <- switch(col_by,
                   mesodermType = print(~dissection),
                   cluster = print(~cluster))
  
  plot_ly(data = df,
          type = "box",
          x = column, y = ~log2.exp,
          color = column, colors = colors,
          width = 350, height = 310
  )  %>%
    layout(
      autosize  = FALSE,
      title = ann_ref[id,1],
      xaxis = c(list(title = "")),
      yaxis = c(list(title = "", showline=T, range=c(0, ceiling(max(df$log2.exp))))),
      margin = list(l=10, r=10, b=0, t=30),
      showlegend = FALSE
    )
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
# diffMap <- as.data.frame(diffMap@eigenvectors[,1:2])
# diffMap$col <- clusters[row.names(diffMap),]$col

# pseudotime
# x <- read.table(paste0(dir, "results/diffusionPseudotime_Me5.tsv"), stringsAsFactors = FALSE, row.names = 1)
# pseudotime.me5 <- x$V2
# names(pseudotime.me5) <- row.names(x)
# x <- read.table(paste0(dir, "results/diffusionPseudotime_Me7.tsv"), stringsAsFactors = FALSE, row.names = 1)
# pseudotime.me7 <- x$V2
# names(pseudotime.me7) <- row.names(x)

# dynamic genes
# x <- read.table(paste0(dir, "results/deltaAIC_Me5.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me5 <- x$V2
# names(dynGenes.me5) <- row.names(x)
# dynGenes.me5 <- dynGenes.me5[dynGenes.me5 < -200]
# x <- read.table(paste0(dir, "results/deltaAIC_Me7.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me7 <- x$V2
# names(dynGenes.me7) <- row.names(x)
# dynGenes.me7 <- dynGenes.me7[dynGenes.me7 < -200]

# clusters
# x <- read.table(paste0(dir, "results/dynGenes_Me5.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me5 <- data.frame(delta = dynGenes.me5, cluster = x[names(dynGenes.me5),1])
# profiles <- c("up", "down", "concave_down", "convex_up", "dip", "convex_down", "convex_down", "concave_up")
# names(profiles) <- 1:8
# dynGenes.me5$profile <- profiles[dynGenes.me5$cluster]
# x <- read.table(paste0(dir, "results/dynGenes_Me7.tsv"), stringsAsFactors = FALSE, row.names = 1)
# dynGenes.me7 <- data.frame(delta = dynGenes.me7, cluster = x[names(dynGenes.me7),1])
# # rename clusters by profile
# profiles <- c("up", "down", "convex_up", "convex_down", "convex_down", "concave_up", "concave_down", "concave_up", "concave_down")
# names(profiles) <- 1:9
# dynGenes.me7$profile <- profiles[dynGenes.me7$cluster]

# dynGenes <- data.frame(gene=ann[union(row.names(dynGenes.me5), row.names(dynGenes.me7)),1],
#                        row.names = union(row.names(dynGenes.me5), row.names(dynGenes.me7)))
# dynGenes$me5 <- dynGenes.me5[match(row.names(dynGenes), row.names(dynGenes.me5)),1]
# dynGenes$me7 <- dynGenes.me5[match(row.names(dynGenes), row.names(dynGenes.me7)),1]
# dynGenes$type <- ifelse(is.na(dynGenes$me5), "Me7", ifelse(is.na(dynGenes$me7), "Me5", "both"))
# dynGenes$clusterMe5 <- dynGenes.me5[row.names(dynGenes),3]
# dynGenes$clusterMe7 <- dynGenes.me7[row.names(dynGenes),3]

# curves
# library(locfit)
# curves_me5 <- t(sapply(row.names(dynGenes.me5), function(x){
#   df <- data.frame(pseudotime.me5, countsNorm[x,names(pseudotime.me5)])
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=1, deg=2), data=df)
#   return(predict(fit2, df[,1]))
# }))
# colnames(curves_me5) <- names(pseudotime.me5)
# 
# curves_me7 <- t(sapply(row.names(dynGenes.me7), function(x){
#   df <- data.frame(pseudotime.me7, countsNorm[x,names(pseudotime.me7)])
#   fit2 <- locfit(df[,2]~lp(df[,1], nn=1, deg=2), data=df)
#   return(predict(fit2, df[,1]))
# }))
# colnames(curves_me7) <- names(pseudotime.me7)


## reference dataset
## umap coords and metadata
# sce <- readRDS(paste0(dir, "data/sce_referenceCells_goodQual_clean.NORM.clusters.Rds"))
# umap_ref <- as.data.frame(reducedDim(sce))
# colnames(umap_ref) <- c("x", "y")
# umap_ref$dissection <- sce$regionAnn
# umap_ref[umap_ref$dissection=="cardiacMesoderm",]$dissection <- "ventralMesoderm"
# umap_ref$cluster <- sce$cluster

# relabel dissections
# labs <- c("DM", "PSM", "CrM", "VM")
# names(labs) <- unique(umap_ref$dissection)
# umap_ref$dissection <- labs[umap_ref$dissection]

# renumber clusters as presented in paper
# labs <- paste0("cluster", 1:5)
# names(labs) <- paste0("cluster", c(3,1,2,6,5))
# umap_ref$cluster <- labs[umap_ref$cluster]

# expression data
# countsNorm_ref <- logcounts(sce)
# countsNorm_ref <- countsNorm_ref[rowSums(countsNorm_ref)>0,]

## gene info
# ann_ref <- read.table(paste0(dir, "data/Mus_musculus.GRCm38.87.tsv"), stringsAsFactors = FALSE, sep="\t", header = TRUE, row.names = 1)
# colnames(ann_ref) <- c("gene", "chr", "start", "end", "strand")
# ann_ref <- ann_ref[row.names(countsNorm_ref),]

# randomForest important genes 
# rf <- readRDS(paste0(dir, "results/randomForest_referenceCells.Rds"))
# imp.genes <- rf$importance[order(rf$importance[,7], decreasing = TRUE),][1:500, 7:8]
# imp.genes <- cbind(gene=ann_ref[row.names(imp.genes),1], as.data.frame(imp.genes))

# save(countsNorm, ann, meta, clusters,
#      markers.all, markers.some,
#      diffMap, pseudotime.me5, pseudotime.me7,
#      dynGenes, curves_me5, curves_me7,
#      countsNorm_ref, ann_ref, umap_ref, imp.genes,
#      file=paste0(dir, "shinyApp/data/data.RData"))
# ##############
