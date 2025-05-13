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
  
  RODAF <- reactive({return(FALSE)})
  observe({
    
    # Advanced settings
    observeEvent(input$advancedSettings, {
      
      # Show modal with options
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   h1(strong("Advanced settings")),
                   h4("There are no advanced settings to select now. This option 
                      will become available in an upcoming update.")
                   
                   # # R-ODAF pipeline
                   # br(),
                   # awesomeCheckbox(inputId = "R_ODAF",
                   #                 label = "Run pipeline with R-ODAF parameters", 
                   #                 value = FALSE,
                   #                 status = "warning")
                   )
            
          ) # EO fluidRow
        ) # EO tagList
      )) # EO showModal
      
      
      # RODAF <- reactive({return(input$RODAF)})
      
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