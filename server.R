#==============================================================================#
# Name: server.R
# Description: server of the ArrayAnalysis Shiny app
#==============================================================================#

# Increase connection size to allow for bigger data uploads
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

# Start server
server <- function(input, output, session){
  
  if (!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
  ##############################################################################
  
  # Preparations
  
  ##############################################################################

  # Set options for data upload
  options(shiny.maxRequestSize=125*1024^10)
  
  # Remove loading screen
  hideTab("navbar", target = "loading_panel")
  
  # Hide the RNA-seq (raw) tabs
  hideTab("navbar", target = "panel_upload_rnaseq_raw")
  hideTab("navbar", target = "panel_preprocessing_rnaseq_raw")
  hideTab("navbar", target = "panel_statistics_rnaseq_raw" )
  hideTab("navbar", target = "panel_ORA_rnaseq_raw")
  
  # Hide the RNA-seq (norm) tabs
  hideTab("navbar", target = "panel_upload_rnaseq_norm")
  hideTab("navbar", target = "panel_preprocessing_rnaseq_norm")
  hideTab("navbar", target = "panel_statistics_rnaseq_norm")
  hideTab("navbar", target = "panel_ORA_rnaseq_norm")
  
  # Hide the microarray (raw) tabs
  hideTab("navbar", target = "panel_upload_microarray_raw")
  hideTab("navbar", target = "panel_preprocessing_microarray_raw")
  hideTab("navbar", target = "panel_statistics_microarray_raw" )
  hideTab("navbar", target = "panel_ORA_microarray_raw")
  
  # Hide the microarray (norm) tabs
  hideTab("navbar", target = "panel_upload_microarray_norm")
  hideTab("navbar", target = "panel_preprocessing_microarray_norm")
  hideTab("navbar", target = "panel_statistics_microarray_norm")
  hideTab("navbar", target = "panel_ORA_microarray_norm")
  
  # Check version
  observe({
    versionMessage <- reactive({tryCatch({
      indexFile <- readLines("https://raw.githubusercontent.com/jarnokoetsier/ArrayAnalysis/refs/heads/main/docs/index.html")
      latest_version <- substr((stringr::str_remove(indexFile[stringr::str_detect(indexFile,"Version")][1],".*Version ")),1,5)
      this_version <- substr(ArrayAnalysis_version, nchar(ArrayAnalysis_version)-4, nchar(ArrayAnalysis_version))
      
      if ((!online) & (this_version < latest_version)){
        message <- paste0("<b style='color:red;'>There is a new ArrayAnalysis version available (v",
                          latest_version,
                          "). Click <a href = 'https://arrayanalysis.org/installation', target = '_blank'>here</a> to install the latest version.</b>")
      } else{
        message <- " "
      }
      return(message)
    }, error = function(cond){
      return(" ")
    })})
    
    output$versionMessage <- renderUI({HTML(versionMessage())})
  })
  
  
  observe({
    
    # Advanced settings
    observeEvent(input$advancedSettings, {
      
      # Show modal with options
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "l",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   h1(strong("About"), "ArrayAnalysis"),
                   h4(style="text-align: justify; line-height:25px;",
                      "ArrayAnalysis is designed to make transcriptomic data analysis 
                      accessible to everyone. By providing a user-friendly interface, ArrayAnalysis aims to empower
                      scientists, regardless of their expertise in high-throughput data analysis, 
                      to analyze and interpret transcriptomic experiments."),
                   hr(),
                   h4("You can run four analysis workflows in ArrayAnalysis: "),
                   br(),
                   h4(
                     HTML('<ol>
                                            <li><b>Raw RNA-seq count table (.tsv/.csv):</b>
                         Select <ins>RNA-Seq analysis</ins> and <ins>Raw counts</ins>.</li>
                         <br>
                         <li><b>Processed RNA-seq count table (.tsv/.csv):</b>
                         Select <ins>RNA-Seq analysis</ins> and <ins>Processed counts</ins>.</li>
                         <br>
                   <li><b>Microarray CEL files (.CEL.gz):</b> 
                   Select <ins>Microarray analysis</ins> and <ins>CEL files</ins>.</li>
                   <br>
                        <li><b>Microarray intensity table (.tsv/.csv):</b>
                        Select <ins>Microarray analysis</ins> and <ins>Processed intensities</ins>.</li>
                        </ol>')),
                   hr(),
                   h4(style="text-align: justify; line-height:25px;",
                      "Click on", HTML("<ins>Start Analysis</ins>"), "to upload your own dataset or to explore the app with an example dataset. Please visit our", a("help page",
                                                                                                                                                          href = "https://arrayanalysis.org/help",
                                                                                                                                                          target = "_blank"),"for more information."),
            )
            
          ) # EO fluidRow
        ) # EO tagList
      )) # EO showModal
      
    })
    
    ############################################################################
    
    # Home page
    
    ############################################################################
    
    # Save input values as a reactive value
    observeEvent(input$startAnalysis,{
      
      # Analyse microarray data or RNA-seq data (or go to documentation)
      microarray_or_rnaseq <- reactive({
        req(input$microarray_or_rnaseq)
        return(input$microarray_or_rnaseq)
      })
      
      # Analyse raw or pre-processed data
      raw_or_norm <- reactive({
        if (input$microarray_or_rnaseq == "Microarray"){
          return(input$raw_or_norm_microarray)
        }
        if (input$microarray_or_rnaseq == "RNA-Seq"){
          return(input$raw_or_norm_rnaseq)
        }
      })
      
      
      ############################################################################
      
      # RNA-seq Analysis
      
      ############################################################################
      
      if (microarray_or_rnaseq() == "RNA-Seq"){
        
        
        #************************************************************************#
        # Raw RNA-seq data
        #************************************************************************#
        if (raw_or_norm() == "Raw data"){
          source("ServerFiles/server_rnaseq_raw.R", local = TRUE)
        }
        
        #************************************************************************#
        # Preprocessed RNA-seq data
        #************************************************************************#
        if (raw_or_norm() == "Processed data"){
          source("ServerFiles/server_rnaseq_norm.R", local = TRUE)
        }
        
      } # EO Microarray or RNA-seq
      
      
      
      ############################################################################
      
      # Microarray Analysis
      
      ############################################################################
      
      if (microarray_or_rnaseq() == "Microarray"){
        
        
        #************************************************************************#
        # Raw microarray data
        #************************************************************************#
        if (raw_or_norm() == "Raw data"){
          source("ServerFiles/server_microarray_raw.R", local = TRUE)
        } # raw or norm
        
        
        #************************************************************************#
        # Preprocessed microarray data
        #************************************************************************#
        if (raw_or_norm() == "Processed data"){
          source("ServerFiles/server_microarray_norm.R", local = TRUE)
        } 
        
      } # EO Microarray or RNA-seq
    }) # EO observeEvent
  }) # EO Observe
} # server