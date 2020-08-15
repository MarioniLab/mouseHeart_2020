## ui

library(shiny)
library(shinycssloaders)
library(DT)
library(plotly)
library(ggplot2)
library(ggpubr)
library(grid)
source("helper.R")

shinyUI(fluidPage(
  titlePanel(h4("Single-cell atlas of the mouse embryonic heart")),
  tabsetPanel(
    tabPanel("About",
      fluidRow(
        column(12,
          br(),
          HTML( "<p style='font-size:15px'><b>Website accompanying the paper </b>
                <i>Anatomical and Transcriptional identification of cardiac 
                progenitors during heart formation reveals a novel heart field,</i>
                by Tyser et al., xxx, 2020.</p><br>" ),
          hr()
        )
      ),
      fluidRow(
        column(7,
          HTML( "<p>Each of the tabs allows exploring the data presented in the manuscript.</p>" ),
          HTML("<p><li><b>Videos: </b>3D volume rendering videos are available for a selection of 
               marker genes. The one in this page shows the expression domain of the canonical cardiac 
               marker NKX2-5. Additional genes are available in the next tab.</p>"),
          HTML( "<p><li><b>Gene expression: </b>type in a gene of interest to plot its expression
                in the UMAP low dimensional space.</li></p>" ),
          HTML( "<p><li><b>Gene expression per stage: </b>type in a gene of interest to plot its expression
                in the UMAP low dimensional space, with cells from different stages plotted separately.</li></p>" ),
          HTML( "<p><li><b>Cluster marker genes: </b>select a cluster of interest to display the genes that show
                preferential expression in that cluster over others. Clicking on a row in the table will 
                show the corresponding expression plot. On the left, you can select the type of test, which 
                controls how many clusters (all or 75%) have lower expression of the gene, compared to the 
                selected cluster.</li></p>" ),
          HTML( "<p><li><b>Diffusion maps: </b>type in a gene of interest to plot its expression
                in the diffusion map low dimensional space. Checking the 'Show dynamic genes' box instead
                displays a table with the genes that change expression along pseudotime. Clicking on a row
                in the table will show the corresponding expression plot. The table can be filtered (using the boxes 
                at the top of each column) to subset to genes with dynamic expression in one or both of the trajectories. 
                Genes can also be filtered based on their expression profile, according to the clusters shown
                on the left. </li></p>" ),
          HTML( "<p><li><b>Reference dataset: </b>type in a gene of interest to plot its expression
                in the UMAP low dimensional space, for the reference dataset produced from microdissection
                of different mesoderm regions. Checking the 'Show most important genes for random forest'
                displays a table with the 500 most important genes for the random forest built from this dataset.
                Clicking on a row in the table will show the corresponding expression plot.</li></p><br>" )
        ),
        column(4,
          # div(style = "margin-top:-80px"),
          uiOutput("video_nkx"),
          HTML("<b>Movie S1:</b> 3D volume rendering movie of a <b>stage 1</b> cardiac crescent using whole mount 
               immunofluorescence to label the cardiac precursor marker <b>NKX2-5 (red)</b>. 
               DAPI labelling of nuclei shown in grey. Related to Figure 1B.")
        )
      ),
      fluidRow(
        column(11,
          hr(),
          HTML( "<p>Raw data is available from the European Nucleotide Archive database under study 
                <a href='https://www.ebi.ac.uk/ena/data/view/PRJEB14363'>PRJEB14363</a>, and ArrayExpress under accession 
                <a href='https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7403/'>E-MTAB-7403</a>. Identifiers 
                for each sample are provided in Supplementary Tables 5 and 7, along with relevant metadata.</p>" ),
          HTML( "<p>Processed data can be downloaded from 
                <a href='https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/'>here</a>.</p>" ),
          HTML( "<p>Code associated with this study is available in 
                <a href='https://github.com/MarioniLab/mouseHeart_2020/'>GitHub</a>.</p>" ),
          hr(),
          br(),
          br()
        )
      )
    ),
    tabPanel("Videos",
      fluidRow(
        column(4, offset=1,
          div(style = "height:20px"),
          h4("Cardiac crescent in relation to the endoderm"),
          uiOutput("video_sox17"),
          HTML("<b>Movie S2:</b> 3D volume rendering movie of a <b>stage 1</b> cardiac crescent using whole mount 
           immunofluorescence to label the contractile protein <b> Sarcomeric Alpha-Actinin (green)</b> 
           and the transcription factor <b>SOX17 (magenta)</b>. Sarcomeric Alpha-Actinin marks 
           cardiomyocytes and highlights the cardiac crescent. SOX17 marks the endoderm overlying 
           the heart, and endothelial cells. DAPI labelling of nuclei is shown in grey.
           Related to Figure 1D.")
        ),
        column(2,
          div(style = "height:70px"),
          HTML("<p style='color:#00FF00';><b>ACTN2</b></bp>
               <p style='color:#FF00FF';><b>SOX17</b></bp>")
        )
      ),
      fluidRow(
        hr()
      ),
      fluidRow(
        column(4, offset=1,
          h4("Additional markers"),
          div(style = "margin-top:-10px"),
          selectInput("videoChoice", label = "",
                      choices = list("Movie S3: Tbx1-Asb2-Nkx2.5 - Stage 2" = "tbx2",
                                     "Movie S4: Tbx1-Asb2-Nkx2.5 - Stage -1" = "tbx-1",
                                     "Movie S5: Tbx1-Asb2-Nkx2.5 - Stage LHT" = "tbxlht",
                                     "Movie S6: Mab21l2-Tbx18-Nkx2.5 - Stage 2" = "mab", 
                                     "Movie S7: Fst-Mab21l2-Nkx2.5 - Stage 0" = "fst0", 
                                     "Movie S8: Vsnl1-Mab21l2-Fsd2 - Stage 1" = "vsn1",
                                     "Movie S9: Vsnl1-Mab21l2-Fsd2 - Stage 2" = "vsn2",
                                     "Movie S10: Vsnl1-Mab21l2-Fsd2 - Stage LHT" = "vsnlht",
                                     "Movie S11: Vsnl1-Mab21l2-Nkx2.5 - Stage -1" = "vsn-1",
                                     "Movie S12:  Me5 cell 1 in Nkx2-5-Cre; R26-nTnG embryo" = "timeLapse1",
                                     "Movie S13:  Me5 cell 2 in Nkx2-5-Cre; R26-nTnG embryo" = "timeLapse2",
                                     "Movie S14: Fst-Mab21l2-Nkx2.5 - Stage 2" = "fst2"),
                      selected = 1,
                      width = "400px"),
          uiOutput("videos"),
          htmlOutput("caption")
        ),
        column(1,
          div(style = "height:110px"),
          htmlOutput("annotation")
        ),
        column(5,
          div(style = "height:50px"),
          img(src = "clusterMarkers.png", height = 420, width = 520),
          br(),
          br(),
          br(),
          br()
        )
      ),
      fluidRow(
        hr()
      ),
      fluidRow(
        column(10, offset = 1,
          h4("Additional data"),
          HTML("TIFF files containing individual Z-sections for 3D whole-mount immunostaining data 
               are available to download, for the figure panels indicated. Channel names and 
               magnification are specified in the file name. Actin was stained using Phalloidin. 
               Differential interference contrast (DIC) was also included to provide a brightfield 
               image of the sample. For more details please refer to the corresponding figure legend 
               and methods."),
          div(style = "height:20px"),
          downloadButton("figure6B", label = "Figure 6B"),
          downloadButton("figure6C", label = "Figure 6C"),
          downloadButton("figure6D", label = "Figure 6D"),
          downloadButton("figure6E", label = "Figure 6E"),
          downloadButton("figure6F", label = "Figure 6F"),
          downloadButton("figureS20E", label = "Supp Figure 20E"),
          downloadButton("figureS20F", label = "Supp Figure 20F"),
          br(),
          br(),
          br(),
          br()
        )
      )
    ),
    tabPanel("Gene expression",
      fluidRow(
        column(3,
          div(style = "height:20px"),
          selectizeInput("gene", label = "Type gene of interest", 
                         choices = NULL, size=10,
                         options = list( placeholder = 'Please type a gene name' )
          ),
          img(src = "UMAP.png", height = 220, width = 270)
        ),
        column(4,
          div(style = "height:70px"),
          withSpinner( plotlyOutput("expression_UMAPplot"), type = 6, size = 0.5)
        ),
        column(4,
          div(style = "height:70px"),
          withSpinner( plotlyOutput("expression_boxplot"), type = 6, size = 0.5)
        )
      ),
      fluidRow(
        column(4, offset=3,
          div(style = "margin-top:-50px"),
          h6("Number of cells positive for the selected gene"),
          tableOutput("numPositiveCells"),
          br(),
          br()
        )
      )
    ),
    tabPanel("Gene expression per stage",
      fluidRow(
        column(3,
          div(style = "height:20px"),
          selectizeInput("gene2", label = "Type gene of interest", 
                         choices = NULL, size=10,
                         options = list( placeholder = 'Please type a gene name' )
          ),
          img(src = "UMAP.png", height = 220, width = 270),
          div(style = "height:10px"),
          downloadButton('downloadPlotStage', 'Download expression plot')
        ),
        column(8,
          fluidRow(
            div(style = "height:80px"),
            plotOutput("UMAPplotPerStage_type")
          ),
          fluidRow(
            div(style = "margin-top:-225px"),
            withSpinner( plotOutput("UMAPplotPerStage_expr"), type = 6, size = 0.5)
          )
        )
      )
    ),
    tabPanel("Cluster marker genes",
      fluidRow(
        column(3,
          div(style = "height:5px"),
          selectInput("selectCluster", label = h5("Cluster"),
                      choices = list( "En1"=3,"En2"=5,"Ec1"=7,"Ec2"=9,
                                      "Me1"=12,"Me2"=11,"Me3"=1,"Me4"=6,"Me5"=4,
                                      "Me6"=10,"Me7"=2,"Me8"=8),
                      selected = 3
          ),
          img(src = "UMAP.png", height = 220, width = 270),
          radioButtons("test", label = h6("Type of test"), 
                       choices = c('all tests significant', 'over 75% of tests significant'), 
                       selected = 'all tests significant'),
          div(style = "height:10px"),
          downloadButton('downloadClusterMarkers', 'Download marker table')
        ),
        column(4, 
          div(style = "height:10px"),
          h6('Click on a row to plot gene expression'),
          withSpinner( DT::dataTableOutput("clusterMarkerTable"), type = 6, size = 0.5)
        ),
        column(1,
          h6('')
        ),
        column(4,
          fluidRow(
            div(style = "height:30px"),
            withSpinner( plotlyOutput("umapClusterMarker"), type = 6, size = 0.5)
          )
        )
      )
    ),
    tabPanel("Diffusion maps",
      fluidRow(
        column(3,
          div(style = "height:20px"),
          selectizeInput("geneDM", label = "Type gene of interest", 
                         choices = NULL, size=10,
                         options = list( placeholder = 'Please type a gene name' )
          ),
          img(src = "diffMap.png", height = 200, width = 200),
          div(style = "height:10px"),
          checkboxInput("diffOnly", label = "Show dynamic genes", value = FALSE, width = NULL),
          conditionalPanel(
            condition = "input.diffOnly == 1",
            img(src = "curves.png", height = 220, width = 220),
            div(style = "height:20px"),
            downloadButton('downloadDynGenes', 'Download dynamic gene table')
          )
        ),
        column(9,
          conditionalPanel(
            condition = "input.diffOnly == 0",
            div(style = "height:60px"),
            withSpinner( plotlyOutput("diffMapPlot"), type = 6, size = 0.5)
          ),
          conditionalPanel(
            condition = "input.diffOnly == 1",
              div(style = "height:15px"),
              h6('Click on the table below to plot gene expression'),
              div(style = "height:30px"),
              withSpinner( plotlyOutput("diffMapPlot_fromTable"), type = 6, size = 0.5),
              div(style = "margin-top:-100px"),
              DT::dataTableOutput("dynamicGenesTable"),
              div(style = "height:15px"),
              h6('Click on the boxes at the top of each column to filter the table by specific categories.'),
              br(),
              br()
          )
        )
      )
    ),
    tabPanel("Reference dataset",
      fluidRow(
        column(3,
          div(style = "height:20px"),
          selectizeInput("gene_ref", label = "Type gene of interest", 
                         choices = NULL, size=10,
                         options = list( placeholder = 'Please type a gene name' )
          ),
          div(style = "height:10px"),
          conditionalPanel(
            condition = "input.colour_by == 'mesodermType'",
            img(src = "UMAP_ref.png", height = 155, width = 290)
          ),
          conditionalPanel(
            condition = "input.colour_by == 'cluster'",
            img(src = "UMAP_ref_cluster.png", height = 155, width = 290)
          ),
          div(style = "height:20px"),
          radioButtons("colour_by", label = "Colour by:",
                       choices = c("mesodermType", "cluster"), 
                       inline = TRUE),
          div(style = "height:20px"),
          checkboxInput("importantGenes", label = "Show most important genes for random forest", 
                        value = FALSE, width = NULL)
        ),
        conditionalPanel(
          condition = "input.importantGenes == 0",
          column(4,
            div(style = "height:70px"),
            withSpinner( plotlyOutput("expression_UMAPplot_reference"), type = 6, size = 0.5)
          ),
          column(4,
            div(style = "height:70px"),
            withSpinner( plotlyOutput("expression_boxplot_reference"), type = 6, size = 0.5)
          )
        ),
        conditionalPanel(
          condition = "input.importantGenes == 1",
          column(4,
            div(style = "height:15px"),
            h6('Click on the table below to plot gene expression'),
            div(style = "height:25px"),
            withSpinner( plotlyOutput("expression_UMAPplot_reference_fromTable"), type = 6, size = 0.5)
          ),
          column(4,
            div(style = "height:70px"),
            withSpinner( plotlyOutput("expression_boxplot_reference_fromTable"), type = 6, size = 0.5)
          )
        )
      ),
      fluidRow(
        conditionalPanel(
          condition = "input.importantGenes == 0",
          column(4, offset=3,
            div(style = "margin-top:-50px"),
            h6("Number of cells positive for the selected gene"),
            tableOutput("numPositiveCells_reference"),
            br(),
            br()
          )
        ),
        conditionalPanel(
          condition = "input.importantGenes == 1",
          column(5, offset=3,
            div(style = "margin-top:-50px"),
            DT::dataTableOutput("ImportantGenesTable"),
            br(),
            br()
          )
        )
      )
    )
  )
))
