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
      tags$video(src = "Nkx25_Stage1.webm", type = "video/webm", 
                 controls = "controls", width = 400, height = 230)
    })
    
    output$video_sox17 <- renderUI({
      tags$video(src = "Actinin2_Sox17_Stage1.webm", type = "video/webm", 
                 controls = "controls", width = 400, height = 230)
    })
    
    output$videos <- renderUI({
      switch(input$videoChoice,
             "tbx2" = tags$video(src = "Tbx1_Asb2_Nkx25_Stage2.webm", 
                                 type = "video/webm", controls = "controls", width = 400, height = 230),
             "tbx-1" = tags$video(src = "Tbx1_Asb2_Nkx25_Stage-1.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "tbxlht" = tags$video(src = "Tbx1_Asb2_Nkx25_LHT.webm", 
                                   type = "video/webm", controls = "controls", width = 400, height = 230),
             "mab" = tags$video(src = "Mab21l2_Tbx18_Nkx25_Stage2.webm", 
                                type = "video/webm", controls = "controls", width = 400, height = 230),
             "fst0" = tags$video(src = "Fst_Mab21l2_Nkx25_Stage0.webm", type = "video/webm", 
                                 controls = "controls", width = 400, height = 230),
             "vsn1" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_Stage1.webm", 
                                 type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsn2" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_Stage2.webm", 
                                 type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsnlht" = tags$video(src = "Vsnl1_Mab21l2_Fsd2_LHT.webm", 
                                   type = "video/webm", controls = "controls", width = 400, height = 230),
             "vsn-1" = tags$video(src = "Vsnl1_Mab21l2_Nkx25_Stage-1.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "timeLapse1" = tags$video(src = "timeLapse_1.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "timeLapse2" = tags$video(src = "timeLapse_2.webm", 
                                  type = "video/webm", controls = "controls", width = 400, height = 230),
             "fst2" = tags$video(src = "Fst_Mab21l2_Nkx25_Stage2.webm", 
                                 type = "video/webm", controls = "controls", width = 400, height = 230)
      )
    })
    
    output$caption <- renderUI({
      caption <- switch(input$videoChoice,
                  'tbx2' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Tbx1</i> (green)</b>, 
                              <b><i>Asb2</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 2G."),
                  'tbx-1' = HTML("3D volume rendering movie of a <b>stage -1</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Tbx1</i> (green)</b>, 
                                 <b><i>Asb2</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 6B."),
                  'tbxlht' = HTML("3D volume rendering movie of a <b>stage LHT</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Tbx1</i> (green)</b>, 
                               <b><i>Asb2</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Supplementary Figure 6D."),
                  'mab' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                             multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Mab21l2</i> (green)</b>, 
                               <b><i>Tbx18</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 3E."),
                  'fst0' = HTML("3D volume rendering movie of a <b>stage 0</b> cardiac crescent using 
                              multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Fst</i> (green)</b>, 
                              <b><i>Mab21l2</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 3F."),
                  'vsn1' = HTML("3D volume rendering movie of a <b>stage 1</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Vsnl1</i> (green)</b>, 
                                <b><i>Mab21l2</i> (red)</b> and <b><i>Fsd2</i> (magenta)</b>. Related to Figure 4B."),
                  'vsn2' = HTML("3D volume rendering movie of a <b>stage 2</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Vsnl1</i> (green)</b>, 
                                <b><i>Mab21l2</i> (red)</b> and <b><i>Fsd2</i> (magenta)</b>. Related to Supplementary Figure 14B."),
                  'vsnlht' = HTML("3D volume rendering movie of a <b>stage LHT</b> cardiac crescent using 
                               multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Vsnl1</i> (green)</b>, 
                                  <b><i>Mab21l2</i> (red)</b> and <b><i>Fsd2</i> (magenta)</b>. Related to Supplementary Figure 14C."),
                  'vsn-1' = HTML("3D volume rendering movie of a <b>stage -1</b> cardiac crescent using 
                               whole mount immunofluorescence to label <b><i>Vsnl1</i> (green)</b>, <b><i>Mab21l2</i>  
                                 (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 4D."),
                  'timeLapse1' = HTML("Time-lapse video of a <b>Nkx2-5-Cre; R26-nTnG</b> embryo imaged over a 16-hour period 
                                      with adaptive light-sheet microscopy showing an Me5 cell (circled) migrating towards the 
                                      heart from a region rostral to the developing cardiac crescent. As cells migrated towards 
                                      the cardiac crescent, the cardiac progenitor marker <i>Nkx2-5</i> was upregulated, 
                                      as observed by an increase in the expression of nuclear GFP. Circles indicate the same cell 
                                      tracked across the time-series shown in different channels."),
                  'timeLapse2' = HTML("Time-lapse video of a <b>Nkx2-5-Cre; R26-nTnG</b> embryo imaged over a 16-hour period 
                                      with adaptive light-sheet microscopy showing an Me5 cell (circled) migrating towards the 
                                      heart from a region rostral to the developing cardiac crescent. As cells migrated towards 
                                      the cardiac crescent, the cardiac progenitor marker <i>Nkx2-5</i> was upregulated, 
                                      as observed by an increase in the expression of nuclear GFP. Circles indicate the same cell 
                                      tracked across the time-series shown in different channels."),
                  'fst2' = HTML("3D volume rendering movie of an <b>early headfold (EHF)</b> embryo using 
                              multiplexed in situ Hybridization Chain Reaction (HCR) to label <b><i>Fst</i> (green)</b>, 
                              <b><i>Mab21l2</i> (red)</b> and <b><i>Nkx2-5</i> (magenta)</b>. Related to Figure 5E.")
                  )
    })
    
    output$annotation <- renderUI({
      ann <- switch(input$videoChoice,
                    'fst0' = HTML("<p style='color:#00FF00';><b><i>Fst</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                  <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'fst2' = HTML("<p style='color:#00FF00';><b><i>Fst</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                  <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'mab' = HTML("<p style='color:#00FF00';><b><i>Mab21l2</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Tbx18</i></b></bp>
                                 <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'tbx-1' = HTML("<p style='color:#00FF00';><b><i>Tbx1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Asb2</i></b></bp>
                                   <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'tbx2' = HTML("<p style='color:#00FF00';><b><i>Tbx1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Asb2</i></b></bp>
                                  <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'tbxlht' = HTML("<p style='color:#00FF00';><b><i>Tbx1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Asb2</i></b></bp>
                                    <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'vsn-1' = HTML("<p style='color:#00FF00';><b><i>Vsnl1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                   <p style='color:#FF00FF';><b><i>Nkx2-5</i></b></bp>"),
                    'vsn1' = HTML("<p style='color:#00FF00';><b><i>Vsnl1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                  <p style='color:#FF00FF';><b><i>Fsd2</i></b></bp>"),
                    'vsn2' = HTML("<p style='color:#00FF00';><b><i>Vsnl1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                  <p style='color:#FF00FF';><b><i>Fsd2</i></b></bp>"),
                    'vsnlht' = HTML("<p style='color:#00FF00';><b><i>Vsnl1</i></b></bp>
                                  <p style='color:#FF0000';><b><i>Mab21l2</i></b></bp>
                                    <p style='color:#FF00FF';><b><i>Fsd2</i></b></bp>")
                    )
    })
    
    ## TIFF downloads
    output$figure6B <- downloadHandler(
      filename = "Figure6B_E80_40x_cTnT_ER_YFP_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Figure6B_E80_40x_cTnT_ER_YFP_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figure6C <- downloadHandler(
      filename = "Figure6C_E85_40x_SarcoActinin_YFP_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Figure6C_E85_40x_SarcoActinin_YFP_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figure6D <- downloadHandler(
      filename = "Figure6D_E95_5x_WT1_cTnT_YFP_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Figure6D_E95_5x_WT1_cTnT_YFP_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figure6E <- downloadHandler(
      filename = "Figure6E_E105_10x_SarcoActinin_YFP_Actin_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Figure6E_E105_10x_SarcoActinin_YFP_Actin_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figure6F <- downloadHandler(
      filename = "Figure6F_E105_10x_Wt1_YFP_Actin_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Figure6F_E105_10x_Wt1_YFP_Actin_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figureS20E <- downloadHandler(
      filename = "Supp20E_E95_10x_Tile5x4_WT1_YFP_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Supp20E_E95_10x_Tile5x4_WT1_YFP_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )
    output$figureS20F <- downloadHandler(
      filename = "Supp20F_E105_SarcoActinin_YFP_Actin_DAPI_DIC.tif.zip",
      content <- function(file) {
        file.copy("www/Supp20F_E105_SarcoActinin_YFP_Actin_DAPI_DIC.tif.zip", file)
      },
      contentType = "application/zip"
    )

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

