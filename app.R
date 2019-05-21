library(shiny)
library(ggplot2)
library(pheatmap)
library(viridis)
library(plotly)
library(RColorBrewer)
library(shinycssloaders)
library(shinythemes)
library(data.table)

# Define UI ----
ui <- fluidPage(
  
  navbarPage("TCGA110 Web App", 
             
             tabPanel("About", 
                      tags$h2("A comprehensive transcriptomic analysis of cell lines as models of primary tumor samples across 22 tumor types"),
                      div(HTML("Katharine Yu<sup>1</sup>, Bin Chen<sup>1,3</sup>, Dvir Aran<sup>1</sup>, Theodore Goldstein<sup>1</sup>, Marina Sirota<sup>1,2</sup>")),
                      div(HTML('<font size="2">1. Institute for Computational Health Sciences, University of California, San Francisco, CA, USA</font>')),
                      div(HTML('<font size="2">2. Department of Pediatrics, University of California, San Francisco, CA</font>')),
                      div(HTML('<font size="2">3. Department of Pediatrics and Human Development, Department of Pharmacology and Toxicology, Michigan State University, Grand Rapids, MI, USA</font>')),                
                      tags$hr(),
                      tags$h3("Cancer cell lines are commonly used as models for cancer biology. While they are limited in their ability to capture complex interactions between tumors and their surrounding environment, they are a cornerstone of cancer research and many important findings have been discovered utilizing cell line models. Not all cell lines are appropriate models of primary tumors, however, which may contribute to the difficulty in translating in vitro findings to patients. We present here a comprehensive pan-cancer analysis utilizing approximately 9,000 transcriptomic profiles from The Cancer Genome Atlas and the Cancer Cell Line Encyclopedia to evaluate cell lines as models of primary tumors across 22 different tumor types. Our analysis of the 22 tumor types are available here in our web app as a resource to the cancer research community, and we hope it will allow researchers to select more appropriate cell line models and increase the translatability of in vitro findings.")
             ),
             
             tabPanel("Summary", 
                      tags$h2("Heatmap showing the median correlations between cell lines and primary tumor samples across all 22 tumor types"),
                      withSpinner(plotlyOutput("heat", width = "auto"), type = 5)
             ),
             
             tabPanel("Explore the data", 
                      sidebarPanel(
                        helpText("Correlation analysis of TCGA primary tumor samples and CCLE cell lines using the 5,000 most variable genes"),
                        selectInput("cancer", "Choose a tumor type", 
                                    choices = list("BLCA", "BRCA", "CHOL", "COADREAD", "DLBC", "ESCA", "GBM", "HNSC", "KIRC", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC"))
                      ),
                      mainPanel(
                        
                        tags$h4(textOutput("selected_var"), align = "center"),
                        downloadButton("downloadData_cancer", "Download cell line and primary tumor correlations"),
                        plotOutput("celllinePlot1", width = "100%", height = "100%"),
                        tags$h4(textOutput("selected_var2"), align = "center"),
                        plotOutput("primarytissuePlot1", height = "100%"),
                        tags$h4(textOutput("selected_var3"), align = "center"),
                        plotOutput("heatmapPlot1", height = "100%"),
                        tags$h4(textOutput("selected_var4"), align = "center"),
                        dataTableOutput("table_ccle"),
                        plotOutput("Subtype1", height = "100%"),
                        plotOutput("Subtype2", height = "100%"),
                        plotOutput("Subtype3", height = "100%"),
                        plotOutput("Subtype4", height = "100%"),
                        tags$style(type="text/css",
                                   ".shiny-output-error { visibility: hidden; }",
                                   ".shiny-output-error:before { visibility: hidden; }")
                  
                      )
             ),
             
             tabPanel("TCGA-110",
                      tags$h2("The TCGA-110 Cell Line Panel"),
                      tags$hr(),
                      downloadButton("downloadData", "Download"),
                      tags$h1(" "),
                      column(12,dataTableOutput('table'))
             ),
             tags$head(tags$script('
                                   var dimension = [0, 0];
                                   $(document).on("shiny:connected", function(e) {
                                   dimension[0] = window.innerWidth;
                                   dimension[1] = window.innerHeight;
                                   Shiny.onInputChange("dimension", dimension);
                                   });
                                   $(window).resize(function(e) {
                                   dimension[0] = window.innerWidth;
                                   dimension[1] = window.innerHeight;
                                   Shiny.onInputChange("dimension", dimension);
                                   });
                                   '))
             
             )
             )

# Define server logic ----
server <- function(input, output) {
  
  observeEvent(input$dimension,{
    output$heat <- renderPlotly({
      correlation = data.matrix(readRDS("data/summary_heatmap.rds"))
      nms <- colnames(correlation)
      plot_ly(width = (0.6*as.numeric(input$dimension[1])), height = .5*as.numeric(input$dimension[1]), colors = inferno(50), x = nms, y = nms, z = correlation, 
              key = correlation, type = "heatmap", source = "heatplot") %>%
        layout(xaxis = list(title = "Cell Lines"),
               yaxis = list(title = "Primary Tumors"))
    })
  })
  
  output$celllinePlot1 <- renderImage({
    filename <- normalizePath(file.path('www', paste(input$cancer, "_specific_cell_lines.png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$primarytissuePlot1 <- renderImage({
    filename <- normalizePath(file.path('www', paste(input$cancer, "_by_primary_sites.png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$heatmapPlot1 <- renderImage({
    filename <- normalizePath(file.path('www', paste(input$cancer, "_heatmap_TCGA_CCLE.png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$Subtype1 <- renderImage({
    subtype = unique(na.omit(readRDS(paste("data/subtypes_", input$cancer, ".rds", sep = ""))$published.subtype))[1]
    filename <- normalizePath(file.path('www', paste(input$cancer,"_specific_cell_lines_", subtype, ".png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)

  output$Subtype2 <- renderImage({
    subtype = unique(na.omit(readRDS(paste("data/subtypes_", input$cancer, ".rds", sep = ""))$published.subtype))[2]
    filename <- normalizePath(file.path('www', paste(input$cancer,"_specific_cell_lines_", subtype, ".png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$Subtype3 <- renderImage({
    subtype = unique(na.omit(readRDS(paste("data/subtypes_", input$cancer, ".rds", sep = ""))$published.subtype))[3]
    filename <- normalizePath(file.path('www', paste(input$cancer,"_specific_cell_lines_", subtype, ".png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$Subtype4 <- renderImage({
    subtype = unique(na.omit(readRDS(paste("data/subtypes_", input$cancer, ".rds", sep = ""))$published.subtype))[4]
    filename <- normalizePath(file.path('www', paste(input$cancer,"_specific_cell_lines_", subtype, ".png", sep='')))
    list(src = filename, width = 0.5*as.numeric(input$dimension[1]))
  }, deleteFile = FALSE)
  
  output$selected_var <- renderText({ 
    paste("Correlation analysis of ", (input$cancer), " primary tumors and cell lines using the 5000 most variable genes")
  })
  output$selected_var2 <- renderText({ 
    paste("Correlations between", (input$cancer), " primary tumors and all cell lines \n (separated by cell line tissue of origin)")
  })
  
  output$selected_var3 <- renderText({ 
    paste("Heatmap of correlations between", (input$cancer), " primary tumors and ",  (input$cancer), "cell lines")
  })
  
  output$selected_var4 <- renderText({ 
    if (input$cancer %in% c("BRCA", "COADREAD", "ESCA", "HNSC", "LUAD", "LUSC", "PAAD", "SKCM", "UCEC")){
    paste("Subtype predictions for ", (input$cancer), "cell lines")
    }else{
      invisible()
    }
  })
  
  output$table_ccle <- renderDataTable({
    
    cancer_subtypes = c("BRCA","COADREAD", "ESCA", "HNSC", "LUAD", "LUSC", "PAAD", "SKCM", "UCEC")
    if (input$cancer %in% cancer_subtypes){
      readRDS(paste("data/ccle_ntp_", input$cancer, ".rds", sep = ""))
    }else{
      invisible()
    }
  })
  
  output$table <- renderDataTable({
    readRDS("data/TCGA110.rds")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "TCGA110.csv"
    },
    content = function(file) {
      write.csv(readRDS("data/TCGA110.rds"), file, row.names = FALSE, quote = F)
    }
  )
  
  output$downloadData_cancer <- downloadHandler(
    filename = function() {
      paste(input$cancer, "_correlations.csv", sep = "")
    },
    content = function(file) {
      write.csv(fread(paste("data/Correlation_Matrices/Correlation_Matrix_", input$cancer, ".txt", sep = "")), file, quote = F, row.names = F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
