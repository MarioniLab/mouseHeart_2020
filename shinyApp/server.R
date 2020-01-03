library(shiny)
library(DT)
library(png)
library(ggplot2)
library(ggpubr)
library(grid)
source("helper.R")

shinyServer(
  function(input, output, session) {
    updateSelectizeInput(session = session, inputId = 'gene', choices = c(Choose = '', selectList), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'gene2', choices = c(Choose = '', selectList), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'geneDM', choices = c(Choose = '', selectList), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'geneDMdiff', choices = c(Choose = '', selectListDM), server = TRUE)

    ## geneExpression
    output$myImage <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("UMAP.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 300, height = 250)
    }, deleteFile = TRUE)
    
    output$UMAPplot <- renderPlot({ 
      validate(
        need( input$gene %in% ann$gene, '')
      )
      plotGeneOnUMAP(gene=input$gene)
    }, height=350, width=640)
    
    output$downloadPlot <- downloadHandler(
      filename = function() { paste0(input$gene, '_expression_UMAP.pdf') },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plotGeneOnUMAP(gene=input$gene), device = "pdf", width = 8, height = 5, units = "in")
        dev.off()
      }
    )
    
    
    ## geneExpression per-stage
    output$myImage2 <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("UMAP.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 250, height = 200)
    }, deleteFile = TRUE)
    
    output$myImageStages <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("UMAP_stages.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 850, height = 130)
    }, deleteFile = TRUE)
    
    output$UMAPplotPerStage <- renderPlot({ 
      validate(
        need( input$gene2 %in% ann$gene, '')
      )
      plotGeneOnUMAPperStage(gene=input$gene2)
    }, height=200, width=850)
    
    output$downloadPlotStage <- downloadHandler(
      filename = function() { paste0(input$gene2, '_expression_UMAP_perStage.pdf') },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plotGeneOnUMAPperStage(gene=input$gene2), device = "pdf", width = 13, height = 3, units = "in")
        dev.off()
      }
    )

    ## clusterMarkers
    output$myImage3 <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("UMAP.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 300, height = 250)
    }, deleteFile = TRUE)


    output$clusterMarkerTable <- DT::renderDataTable(
      datatable( printClusterMarkerTable(pop=input$selectCluster, type=input$test), colnames = c("gene", "FDR", "log2 fold-change")),
      server=TRUE
    )

    datasetClusterMarkers <- reactive({
      retrieveClusterMarkers(pop=input$selectCluster, type=input$test)
    })
    output$downloadClusterMarkers <- downloadHandler(
      filename = function() { 
        if(input$test == 'all tests significant'){ 
          paste0('clusterMarkers_', names(markers.all)[as.numeric(input$selectCluster)], '_all.tsv')
        }else{
          paste0('clusterMarkers_', names(markers.some)[as.numeric(input$selectCluster)], '_over0.75.tsv')
        }
      },
      content = function(file) {
        write.table(datasetClusterMarkers(), file, quote = FALSE, sep="\t")
      }
    )

    output$umapClusterMarker <- renderPlot({
      validate(
        if(input$test == "all"){
          need(input$clusterMarkerTable_row_last_clicked <= nrow(markers.all[[as.numeric(input$selectCluster)]]), '')
        }else{
          need(input$clusterMarkerTable_row_last_clicked <= nrow(markers.some[[as.numeric(input$selectCluster)]]), '')
        }
      )
      sel = input$clusterMarkerTable_row_last_clicked
      if(length(sel)){
        plotClusterMarkerGeneOnUMAP(pop=input$selectCluster, type=input$test, sel=sel)
      }
    }, height=550, width=300)


    ## trajectories
    output$myImageDM <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("diffMap.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 285, height = 285)
    }, deleteFile = TRUE)

    output$diffMapPlot <- renderPlot({
      if(input$diffOnly){
        validate(
          need( input$geneDMdiff %in% ann$gene, '')
        )
        plotGeneOnDiffMap(gene=input$geneDMdiff)
      }else{
        validate(
          need( input$geneDM %in% ann$gene, '')
        )
        plotGeneOnDiffMap(gene=input$geneDM)
      }
    }, height=380, width=950)

    output$downloadPlotDM <- downloadHandler(
      filename = function() { paste0(input$geneDM, '_expression_diffusionMap.pdf') },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plotGeneOnDiffMap(gene=input$geneDM), device = "pdf", width = 12, height = 5, units = "in")
        dev.off()
      }
    )
    output$downloadPlotDMdiff <- downloadHandler(
      filename = function() { paste0(input$geneDMdiff, '_expression_diffusionMap.pdf') },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plotGeneOnDiffMap(gene=input$geneDMdiff), device = "pdf", width = 12, height = 5, units = "in")
        dev.off()
      }
    )
    
    datasetDynamicGenes <- reactive({
      retrieveDynamicGenes()
    })
    output$downloadDynGenes <- downloadHandler(
      filename = function() { 
        "dynamicGenes.tsv"
      },
      content = function(file) {
        write.table(datasetDynamicGenes(), file, quote = FALSE, sep="\t")
      }
    )

  }
)

