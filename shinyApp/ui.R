library(shiny)
shinyUI(fluidPage(
  titlePanel(h5("Single-cell RNA-seq data of mouse heart")),
  tabsetPanel(
    tabPanel("Gene expression",
             fluidRow(
                 column(3,
                        selectizeInput("gene", label = h4("Gene of interest"), choices = NULL,
                                       options = list(placeholder = 'type a gene name')),
                        imageOutput("myImage")
                 ),
                 column(4,
                        div(style = "height:65px;"),
                        plotOutput("UMAPplot"),
                        div(style = "height:5px;"),
                        downloadButton('downloadPlot', 'Download plot')
                 )
             )
    ),
    tabPanel("Gene expression per stage",
             fluidRow(
                 column(3,
                        selectizeInput("gene2", label = h4("Gene of interest"), choices = NULL,
                                       options = list(placeholder = 'type a gene name')),
                        imageOutput("myImage2")
                 ),
                 column(8,
                        fluidRow(
                          div(style = "height:50px;"),
                          imageOutput("myImageStages")
                        ),
                        fluidRow(
                          div(style = "margin-top:-225px"),
                          plotOutput("UMAPplotPerStage")
                        ),
                        fluidRow(
                          div(style = "margin-top:-200px"),
                          downloadButton('downloadPlotStage', 'Download plot')
                        )
                 )
             )
    ),
    tabPanel("Cluster marker genes",
               fluidRow(
                 column(3,
                        selectInput("selectCluster", label = h5("Cluster"),
                                    choices = list( "En1"=3,"En2"=5,"Ec1"=7,"Ec2"=9,
                                                    "Me1"=12,"Me2"=11,"Me3"=1,"Me4"=6,"Me5"=4,
                                                    "Me6"=10,"Me7"=2,"Me8"=8),
                                    selected = 3),
                        imageOutput("myImage3", height = "300px"),
                        radioButtons("test", label = h6("Type of test"), choices = c('all tests significant', 'over 75% of tests significant'), selected = 'all tests significant')
                 ),
                 column(4, offset=1,
                        h6('Click on a row to plot gene expression'),
                        DT::dataTableOutput("clusterMarkerTable"),
                        div(style = "height:20px;"),
                        downloadButton('downloadClusterMarkers', 'Download table')
                 ),
                 column(4, offset=0.5,
                        plotOutput("umapClusterMarker"),
                        div(style = "height:200px;")
                 )
             )
    ),
    tabPanel("Diffusion maps",
             fluidRow(
                 column(3,
                        checkboxInput("diffOnly", label = "Show DE genes only", value = FALSE, width = NULL),
                        conditionalPanel(
                          condition = "input.diffOnly",
                          selectizeInput("geneDMdiff", label = h4("Gene of interest"), choices = NULL, selected=NULL,
                                         options = list(placeholder = 'type a gene name'))
                        ),
                        conditionalPanel(
                          condition = "!input.diffOnly",
                          selectizeInput("geneDM", label = h4("Gene of interest"), choices = NULL, selected=NULL,
                                         options = list(placeholder = 'type a gene name'))
                        ),
                        imageOutput("myImageDM"),
                        div(style = "height:-500px;"),
                        downloadButton('downloadDynGenes', 'Download dynamic genes')
                 ),
                 column(8,
                        div(style = "height:115px;"),
                        plotOutput("diffMapPlot"),
                        div(style = "height:20px;"),
                        conditionalPanel(
                          condition = "input.diffOnly",
                          downloadButton('downloadPlotDMdiff', 'Download plot')
                        ),
                        conditionalPanel(
                          condition = "!input.diffOnly",
                          downloadButton('downloadPlotDM', 'Download plot')
                        )
                 )
             )
    )
  )
))
