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
  
  # Checck version
  observe({
    versionMessage <- reactive({tryCatch({
      indexFile <- readLines("https://raw.githubusercontent.com/jarnokoetsier/ArrayAnalysis/refs/heads/main/docs/index.html")
      latest_version <- substr((str_remove(indexFile[str_detect(indexFile,"Version")],".*Version ")),1,5)
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
                   h4(style="text-align: justify;",
                   "ArrayAnalysis is designed to make transcriptomic data analysis 
                      accessible to everyone. With its user-friendly interface, it empowers 
                      scientists, regardless of their expertise in high-throughput data analysis, 
                      to analyze and interpret transcriptomic experiments with ease."),
                   br(),
                   h4("You can run four analysis workflows in ArrayAnalysis: "),
                   h4(
                   HTML('<ol>
                   <li><b>Microarray CEL files:</b> 
                   Select <code>Microarray analysis</code> and <code>CEL files</code>.</li>
                        <li><b>Microarray intensity table:</b>
                        Select <code>Microarray analysis</code> and <code>Processed intensities</code>.</li>
                         <li><b>Raw RNA-seq count table:</b>
                         Select <code>RNA-Seq analysis</code> and <code>Raw counts</code>.</li>
                         <li><b>Processed RNA-seq count table:</b>
                         Select <code>RNA-Seq analysis</code> and <code>Processed counts</code>.</li>
                        </ol>')),
                   br(),
                   h4(style="text-align: justify;",
                     "Click on", code("Start Analysis"), "to upload your own dataset or run explore the app with an example dataset. Please visit our", a("help page",
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