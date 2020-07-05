## server

shinyServer(
  
  function(input, output, session) {
    ## populate gene lists
    updateSelectizeInput(session = session, inputId = 'gene', choices = selectList, 
                         selected = 'Nkx2-5', server = TRUE)
    updateSelectizeInput(session = session, inputId = 'gene2', choices = selectList, 
                         selected = 'Nkx2-5', server = TRUE)
    updateSelectizeInput(session = session, inputId = 'geneDM', choices = selectList, 
                         selected = 'Nkx2-5', server = TRUE)
    updateSelectizeInput(session = session, inputId = 'gene_ref', choices = selectList_ref, 
                         selected = 'Nkx2-5', server = TRUE)
    
    ## videos
    output$video_nkx <- renderUI({
      tags$video(src = "Nkx25_Stage1_Figure1B.webm", type = "video/webm", 
                 controls = "controls", width = 400, height = 230)
    })
    
    output$video_sox17 <- renderUI({
      tags$video(src = "Actinin2_Sox17_Stage1_Figure1D.webm", type = "video/webm", 
                 controls = "controls", width = 400, height = 230)
    })
    
    output$videos <- renderUI({
      switch(input$videoChoice,
             "fst0" = tags$video(src = "Fst_Mab21l2_Nkx25_Stage0_Figure2I.webm", type = "video/webm", 
                                controls = "controls", width = 400, height = 230),
             "fst2" = tags$video(src = "Fst_Mab21l2_Nkx25_Stage2_Figure5E.webm", type = "video/webm", 
                                 controls = "controls", width = 400, height = 230),
             "mab" = tags$video(src = "Mab21l2_Tbx18_Nkx25_Stage2_Figure2B.webm", 
                                type = "video/webm", controls = "controls", width = 400, height = 230),
             "tbx-1" = tags$video(src = "Tbx1_Asb2_Nkx25_Stage-1_Figure2G.webm", 
                                type = "video/webm", controls = "controls", width = 400, height = 230),
             "tbx2" = tags$video(src = "Tbx1_Asb2_Nkx25_Stage2_Figure2H.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "tbxlht" = tags$video(src = "Tbx1_Asb2_Nkx25_LHT_FigureSupp6C.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsn-1" = tags$video(src = "Vsnl1_Mab21l2_Nkx25_Stage-1_Figure4D.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsn1" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_Stage1_SuppFigure11A.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsn2" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_Stage2_Figure4C.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsnlht" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_LHT_SuppFigure11C.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230)
      )
    })
    
    output$caption <- renderUI({
      caption <- switch(input$videoChoice,
                  'fst0' = HTML("3D volume rendering movie of a <b>stage 0</b> cardiac crescent using 
                              multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Fst (green)</b>, <b>Mab21l2 
                              (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 2I."),
                  'fst2' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                              multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Fst (green)</b>, <b>Mab21l2 
                              (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 5E."),
                  'mab' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                             multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Mab21l2 (green)</b>, <b>Tbx18 
                             (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 2B."),
                  'tbx-1' = HTML("3D volume rendering movie of a <b>stage -1</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Tbx1 (green)</b>, <b>Asb2 
                               (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 2G."),
                  'tbx2' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Tbx1 (green)</b>, <b>Asb2 
                               (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 2H."),
                  'tbxlht' = HTML("3D volume rendering movie of a <b>stage LHT</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Tbx1 (green)</b>, <b>Asb2 
                               (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Supplementary Figure 6C."),
                  'vsn-1' = HTML("3D volume rendering movie of a <b>stage -1</b> cardiac crescent using 
                               whole mount immunofluorescence to label <b> Vsnl1 (green)</b>, <b>Mab21l2 
                               (red)</b> and <b>Nkx2-5 (magenta)</b>. Related to Figure 4D."),
                  'vsn1' = HTML("3D volume rendering movie of a <b>stage 1</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Vsnl1 (green)</b>, <b>Mab21l2 
                               (red)</b> and <b>Fsd2 (magenta)</b>. Related to Supplementary Figure 11A."),
                  'vsn2' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Vsnl1 (green)</b>, <b>Mab21l2 
                              (red)</b> and <b>Fsd2 (magenta)</b>. Related to Figure 4C."),
                  'vsnlht' = HTML("3D volume rendering movie of a <b>stage LHT</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b> Vsnl1 (green)</b>, <b>Mab21l2 
                              (red)</b> and <b>Fsd2 (magenta)</b>. Related to Supplementary Figure 11C.")
                  )
    })
    
    output$annotation <- renderUI({
      ann <- switch(input$videoChoice,
                    'fst0' = HTML("<p style='color:#00FF00';><b>Fst</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                  <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'fst2' = HTML("<p style='color:#00FF00';><b>Fst</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                  <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'mab' = HTML("<p style='color:#00FF00';><b>Mab21l2</b></bp>
                                  <p style='color:#FF0000';><b>Tbx18</b></bp>
                                 <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'tbx-1' = HTML("<p style='color:#00FF00';><b>Tbx1</b></bp>
                                  <p style='color:#FF0000';><b>Asb2</b></bp>
                                   <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'tbx2' = HTML("<p style='color:#00FF00';><b>Tbx1</b></bp>
                                  <p style='color:#FF0000';><b>Asb2</b></bp>
                                  <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'tbxlht' = HTML("<p style='color:#00FF00';><b>Tbx1</b></bp>
                                  <p style='color:#FF0000';><b>Asb2</b></bp>
                                    <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'vsn-1' = HTML("<p style='color:#00FF00';><b>Vsnl1</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                   <p style='color:#FF00FF';><b>Nkx2-5</b></bp>"),
                    'vsn1' = HTML("<p style='color:#00FF00';><b>Vsnl1</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                  <p style='color:#FF00FF';><b>Fsd2</b></bp>"),
                    'vsn2' = HTML("<p style='color:#00FF00';><b>Vsnl1</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                  <p style='color:#FF00FF';><b>Fsd2</b></bp>"),
                    'vsnlht' = HTML("<p style='color:#00FF00';><b>Vsnl1</b></bp>
                                  <p style='color:#FF0000';><b>Mab21l2</b></bp>
                                    <p style='color:#FF00FF';><b>Fsd2</b></bp>")
                    )
    })
    

    ## geneExpression
    output$expression_UMAPplot <- renderPlotly({
      validate(
        need( input$gene %in% ann$gene, '')
      )
      plotlyGeneOnUMAP(gene=input$gene)
    })
    
    output$expression_boxplot <- renderPlotly({
      validate(
        need( input$gene %in% ann$gene, '')
      )
      expr_boxplot(gene=input$gene)
    })
    
    output$numPositiveCells <- renderTable({
      validate(
        need( input$gene %in% ann$gene, '')
      )
      numberPositiveCells(gene = input$gene)
    }, rownames = TRUE, digits = 1)
    
    
    ## geneExpression per-stage
    output$UMAPplotPerStage_type <- renderPlot({ 
      plotUMAPperStage()
    }, height=175, width=950)
    
    output$UMAPplotPerStage_expr <- renderPlot({ 
      validate(
        need( input$gene2 %in% ann$gene, '')
      )
      plotGeneOnUMAPperStage(gene=input$gene2)
    }, height=225, width=950)
    
    output$downloadPlotStage <- downloadHandler(
      filename = function() { paste0(input$gene2, '_expression_UMAP_perStage.pdf') },
      content = function(file) {
        pdf(file = NULL)
        ggsave(file, plot = plotGeneOnUMAPperStage(gene=input$gene2), 
               device = "pdf", width = 13, height = 3, units = "in")
        dev.off()
      }
    )

    ## clusterMarkers
    output$clusterMarkerTable <- DT::renderDataTable(
      datatable( printClusterMarkerTable(pop=input$selectCluster, type=input$test), 
                 colnames = c("gene", "FDR", "log2 fold-change") 
      ),
      server=TRUE
    )

    datasetClusterMarkers <- reactive({
      retrieveClusterMarkers(pop=input$selectCluster, type=input$test)
    })
    output$downloadClusterMarkers <- downloadHandler(
      filename = function() { 
        if(input$test == 'all tests significant'){ 
          paste0('clusterMarkers_', 
                 names(markers.all)[as.numeric(input$selectCluster)], 
                 '_all.tsv')
        }else{
          paste0('clusterMarkers_', 
                 names(markers.some)[as.numeric(input$selectCluster)], 
                 '_over0.75.tsv')
        }
      },
      content = function(file) {
        write.table(datasetClusterMarkers(), file, quote = FALSE, sep="\t")
      }
    )

    output$umapClusterMarker <- renderPlotly({
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
    })


    ## trajectories
    output$diffMapPlot <- renderPlotly({
        validate( need( input$geneDM %in% ann$gene, '') )
        plotGeneOnDiffMap(gene=input$geneDM)
    })
    
    output$dynamicGenesTable <- DT::renderDataTable(
      datatable( printDynamicGenesTable(trajectory=input$trajectory, cluster=input$cluster), 
                 colnames = c("gene", "trajectory", "cluster_Me5", "cluster_Me7"),
                 options = list(lengthMenu = c(5, 10, 20, 50), pageLength = 5),
                 filter = "top"),
      server=TRUE
    )
    
    output$diffMapPlot_fromTable <- renderPlotly({
      sel <- input$dynamicGenesTable_row_last_clicked
      req(sel)
      df <- dynGenes[,c(1,4:6)]
      validate(need(sel <= nrow(df), ''))
      
      plotGeneOnDiffMap_fromTable(sel=input$dynamicGenesTable_row_last_clicked, 
                                  trajectory=input$trajectory,
                                  cluster=input$cluster)
    })
    
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
    
    ## geneExpression in reference dataset
    output$expression_UMAPplot_reference <- renderPlotly({
      validate(
        need( input$gene_ref %in% ann_ref$gene, '')
      )
      plotlyGeneOnUMAP_ref(gene=input$gene_ref)
    })
    
    output$expression_boxplot_reference <- renderPlotly({
      validate(
        need( input$gene_ref %in% ann_ref$gene, '')
      )
      expr_boxplot_ref(gene=input$gene_ref, col_by=input$colour_by)
    })
    
    output$numPositiveCells_reference <- renderTable({
      validate(
        need( input$gene_ref %in% ann_ref$gene, '')
      )
      numberPositiveCells_ref(gene = input$gene_ref, col_by=input$colour_by)
    }, rownames = TRUE, digits = 1)
    
    output$ImportantGenesTable <- DT::renderDataTable(
      datatable( printImportantGenesTable(), 
                 colnames = c("gene", "MeanDecreaseAccuracy", "MeanDecreaseGini"),
                 options = list(lengthMenu = c(5, 10, 20, 50), pageLength = 5)),
      server=TRUE
    )
    
    output$expression_UMAPplot_reference_fromTable <- renderPlotly({
      sel <- input$ImportantGenesTable_row_last_clicked
      req(sel)
      plotlyGeneOnUMAP_ref_fromTable(sel=sel)
    })
    
    output$expression_boxplot_reference_fromTable <- renderPlotly({
      sel <- input$ImportantGenesTable_row_last_clicked
      req(sel)
      expr_boxplot_ref_fromTable(sel=sel, col_by=input$colour_by)
    })

  }
)

