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
  
  # Make list for reactive values
  rv <- reactiveValues()
  
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
  
  observe({
    
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
        req(input$raw_or_norm)
        return(input$raw_or_norm)
      })
      
      ############################################################################
      
      # RNA-seq Analysis
      
      ############################################################################
      
      if (microarray_or_rnaseq() == "RNA-Seq"){
        
        
        #************************************************************************#
        # Raw RNA-seq data
        #************************************************************************#
        if (raw_or_norm() == "Raw data"){
          
          #======================================================================#
          # TAB2: Data Upload
          #======================================================================#
          
          #----------------------------------------------------------------------#
          # Go to data upload tab
          #----------------------------------------------------------------------#
          
          observeEvent(input$startAnalysis,{
            
            # Show RNA-seq (raw) upload tab
            showTab("navbar", target = "panel_upload_rnaseq_raw")
            
            # Remove the other RNA-seq (raw) tabs
            hideTab("navbar", target = "panel_preprocessing_rnaseq_raw")
            hideTab("navbar", target = "panel_statistics_rnaseq_raw" )
            hideTab("navbar", target = "panel_ORA_rnaseq_raw")
            
            # Remove the other RNA-seq (norm) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_norm")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_norm")
            hideTab("navbar", target = "panel_statistics_rnaseq_norm")
            hideTab("navbar", target = "panel_ORA_rnaseq_norm")
            
            # Remove the other microarray (raw) tabs
            hideTab("navbar", target = "panel_upload_microarray_raw")
            hideTab("navbar", target = "panel_preprocessing_microarray_raw")
            hideTab("navbar", target = "panel_statistics_microarray_raw" )
            hideTab("navbar", target = "panel_ORA_microarray_raw")
            
            # Remove the other microarray (norm) tabs
            hideTab("navbar", target = "panel_upload_microarray_norm")
            hideTab("navbar", target = "panel_preprocessing_microarray_norm")
            hideTab("navbar", target = "panel_statistics_microarray_norm")
            hideTab("navbar", target = "panel_ORA_microarray_norm")
            
            
            # Go to microarray tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_upload_rnaseq_raw")
            
            # Example meta data file
            output$downloadmeta_example_rnaseq_raw <- downloadHandler(
              filename = "MetaData_example.csv",
              content = function(file){
                write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
              }
            )
          })
          
          #----------------------------------------------------------------------#
          # Upload files
          #----------------------------------------------------------------------#
          
          observeEvent(input$upload_rnaseq_raw,{
            
            # Show modal
            show_modal_spinner(text = "Reading data...",
                               color="#0dc5c1")
            
            # Read expression data
            if (!is.null(input$uploadExprData_rnaseq_raw)){
              rv$gxData <- readRNASeq(path=input$uploadExprData_rnaseq_raw$datapath)
            } else{
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
              
            }
            
            # Read metadata
            if (!is.null(input$uploadMeta_rnaseq_raw_tsv)){
              if (input$MetaFileType_rnaseq_raw ==".tsv/.csv file"){
                rv$metaData <- getMetaData(path = input$uploadMeta_rnaseq_raw_tsv$datapath,
                                           celfiles = colnames(rv$gxData),
                                           filetype = input$MetaFileType_rnaseq_raw)
              }
            }else if(!is.null(input$uploadMeta_rnaseq_raw_smf)){ 
              if (input$MetaFileType_rnaseq_raw =="Series Matrix File"){
                rv$metaData <- getMetaData(path = input$uploadMeta_rnaseq_raw_smf$datapath,
                                           celfiles =  colnames(rv$gxData),
                                           filetype = input$MetaFileType_rnaseq_raw)
              }
            }else{
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
              
            }
            
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(rv$gxData)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- rv$gxData[,rownames(rv$metaData)]
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_rnaseq_raw <- renderDataTable(NULL)
              output$metaTable_rnaseq_raw <- renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_rnaseq_raw <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(rv$gxData,6)
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              # output$metaTable_rnaseq_raw <- DT::renderDataTable({
              #   req(rv$metaData)
              #   return(rv$metaData)
              # }, options = list(pageLength = 6))
              
              output$metaTable_rnaseq_raw <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              observeEvent(input$metaTable_rnaseq_raw_cell_edit, {
                row  <- input$metaTable_rnaseq_raw_cell_edit$row
                clmn <- input$metaTable_rnaseq_raw_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_rnaseq_raw_cell_edit$value
              })
              
              
              # Render UI for main tab
              output$UI_upload_rnaseq_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_rnaseq_raw") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_rnaseq_raw") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_rnaseq_raw <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_rnaseq_raw")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(rv$gxData)){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_rnaseq_raw",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # EO observeEvent
          
          #----------------------------------------------------------------------#
          # Run Example
          #----------------------------------------------------------------------#
          
          observeEvent(input$example_rnaseq_raw,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Read expression data
            rv$gxData <- readRNASeq(path="Data/RNAseq/rawExpr_GSE128380.csv")
            
            # Get metadata
            rv$metaData <- getMetaData(path = "Data/RNAseq/sampleInfo_GSE128380.csv",
                                       celfiles = colnames(rv$gxData),
                                       filetype = ".tsv/.csv file")
            
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(rv$gxData)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- rv$gxData[,rownames(rv$metaData)]
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_rnaseq_raw <- DT::renderDataTable(NULL)
              output$metaTable_rnaseq_raw <- DT::renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_rnaseq_raw <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(rv$gxData,6)
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_rnaseq_raw <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              observeEvent(input$metaTable_rnaseq_raw_cell_edit, {
                row  <- input$metaTable_rnaseq_raw_cell_edit$row
                clmn <- input$metaTable_rnaseq_raw_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_rnaseq_raw_cell_edit$value
              })
              
              # Render UI for main tab
              output$UI_upload_rnaseq_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_rnaseq_raw") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Metadata table")),
                             h5("This is a preview of the metadata table. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_rnaseq_raw") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_rnaseq_raw <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_rnaseq_raw")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(rv$gxData)){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_rnaseq_raw",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # EO observeEvent
          
          
          #======================================================================#
          # TAB3: Data Pre-processing
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_upload_rnaseq_raw,{
            
            # Go to RNA-seq tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_preprocessing_rnaseq_raw")
          })
          
          
          # 1. Select outliers
          output$UI_outlier_rnaseq_raw <- renderUI({
            if(!input$outlier_rnaseq_raw){
              samples <- rownames(rv$metaData)
              
              shinyWidgets::pickerInput(inputId = "select_outliers_rnaseq_raw",
                                        label = tags$span(
                                          "Select samples to be removed", 
                                          tags$span(
                                            icon(
                                              name = "question-circle",
                                            ) 
                                          ) |>
                                            prompter::add_prompt(message = "Select one or more samples to exclude from the analysis.", 
                                                                 position = "right",
                                                                 size = "large")
                                        ),
                                        choices = samples,
                                        multiple = TRUE)
            } else{
              NULL
            }
            
          })
          
          # 2. Select experimental groups
          output$UI_groupselect_rnaseq_raw <- renderUI({
            shinyWidgets::pickerInput(inputId = "groupselect_rnaseq_raw",
                                      label = NULL,
                                      choices = colnames(rv$metaData),
                                      selected = autoGroup(rv$metaData),
                                      multiple = TRUE)
          })
          
          
          # print the experimental levels
          output$experimentallevels_rnaseq_raw <- renderText({
            req(input$groupselect_rnaseq_raw)
            
            if(length(input$groupselect_rnaseq_raw) > 1){
              experimentFactor <- make.names(apply(rv$metaData[,input$groupselect_rnaseq_raw], 1, paste, collapse = "_" ))
            } else{
              experimentFactor <- make.names(rv$metaData[,input$groupselect_rnaseq_raw])
            }
            levels <- unique(experimentFactor)
            if (length(levels) <= 6){
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(levels, collapse = ", "),
                             "</p>")
            } else {
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(head(levels), collapse = ", "),
                             "...</p>")
            }
            return(text)
            
          })
          
          
          # 3. pre-process the raw data
          observeEvent(input$start_preprocessing_rnaseq_raw,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Pre-processing data...",
                                          color="#0dc5c1")
            
            hideTab("navbar", target = "panel_statistics_rnaseq_raw")
            hideTab("navbar", target = "panel_ORA_rnaseq_raw")
            rv$top_table <- NULL
            
            # Select outlier
            if (!isTRUE(input$outier_rnaseq_raw)){
              rv$outlier <- input$select_outliers_rnaseq_raw
            } else{
              rv$outlier <- NULL
            }
            
            # Remove outlier
            if (!is.null(rv$outlier)){
              gxMatrix_temp <- rv$gxData
              rv$gxData_fil <- gxMatrix_temp[,setdiff(colnames(gxMatrix_temp),rv$outlier)]
              rm(gxMatrix_temp)
            } else
              rv$gxData_fil <- rv$gxData
            
            # Filter metadata and expression data (samples in correct order)
            rv$metaData_fil <- rv$metaData[colnames(rv$gxData_fil),]
            
            # Experiment factor
            if(length(input$groupselect_rnaseq_raw) > 1){
              rv$experimentFactor <- factor(make.names(apply(rv$metaData_fil[,input$groupselect_rnaseq_raw], 1, paste, collapse = "_" )))
              rv$experimentName <- input$groupselect_rnaseq_raw
            } else{
              rv$experimentFactor <- factor(make.names(rv$metaData_fil[,input$groupselect_rnaseq_raw]))
              rv$experimentName <- input$groupselect_rnaseq_raw
            }
            
            # Get filter
            rv$filterThreshold <- input$countfilter_rnaseq_raw
            
            # Normalization
            rv$normData <- RNASeqNormalization(gxData = rv$gxData_fil,
                                               metaData = rv$metaData_fil,
                                               filterThres = rv$filterThreshold,
                                               smallestGroupSize = min(table(rv$experimentFactor)))
            
            
            # Collect chosen pre-processing settings into a dataframe
            rv$processingSettings <- data.frame(
              Option = c("Removed samples",
                         "Experimental group",
                         "Experimental levels",
                         "Filter threshold"),
              Selected = c(paste(rv$outlier, collapse = "; "),
                           paste(input$groupselect_rnaseq_raw, collapse = "; "),
                           paste(unique(rv$experimentFactor), collapse = "; "),
                           rv$filterThreshold
              )
            )
            #======================#
            # Make QC plots/tables
            #======================#
            
            #********************************************************************#
            # Output 1: Expression values
            #********************************************************************#
            
            # Print expression table
            output$exprTable_rnaseq_raw <- DT::renderDataTable({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show message
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Info",
                text = "The data has been pre-processed. Please check the different 
              QC plots on this page to assess the pre-processing quality.",
                type = "info")
              
              # Show microarray statistics tab
              showTab("navbar", target = "panel_statistics_rnaseq_raw")
              
              output <- rv$normData
              return(output)
              
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$downloadNormalizedData_rnaseq_raw <- downloadHandler(
              filename = "normalizedData.csv",
              content = function(file){
                write.csv(rv$normData, file, quote = FALSE, row.names = TRUE)
              }
            )
            
            # Change color by click on button
            rv$colorOrder <- 1:length(levels(rv$experimentFactor))
            observeEvent(input$geneboxplot_changeOrder_rnaseq_raw,{
              all_orders <- permute(1:length(levels(rv$experimentFactor))) 
              sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
              
              if (sel < length(all_orders)){
                rv$colorOrder <- all_orders[[sel+1]]
              } else{
                rv$colorOrder <- all_orders[[1]]
              }
            })
            
            
            # Boxplot of single gene (based on selected row in the expression matrix)
            output$ExprBoxplot_rnaseq_raw <- renderPlot({
              req(input$exprTable_rnaseq_raw_rows_selected)
              
              if (length(levels(rv$experimentFactor)) > 4){
                legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
              } else{
                legendColors <- c(input$geneboxplot_col1_rnaseq_raw,
                                  input$geneboxplot_col2_rnaseq_raw,
                                  input$geneboxplot_col3_rnaseq_raw,
                                  input$geneboxplot_col4_rnaseq_raw,
                                  input$geneboxplot_col5_rnaseq_raw,
                                  input$geneboxplot_col6_rnaseq_raw)
              }
              # Make boxplot
              rv$temp1 <- geneBoxplot(experimentFactor = rv$experimentFactor, 
                          normMatrix = rv$normData, 
                          sel_row = input$exprTable_rnaseq_raw_rows_selected,
                          legendColors = legendColors[rv$colorOrder],
                          groupOrder = input$geneboxplot_order_rnaseq_raw,
                          rnaseq = TRUE,
                          jitter = input$jitter_geneboxplot_rnaseq_raw,
                          seed = sample(1:1000,1))
              return(rv$temp1)
            })
            
            
            # Get number of experimental groups
            output$length_geneboxplot_rnaseq_raw <- reactive({
              length(levels(rv$experimentFactor))
            })
            outputOptions(output, "length_geneboxplot_rnaseq_raw", suspendWhenHidden = FALSE) 
            
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_geneboxplot_rnaseq_raw <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp1, 
                       filename = file,
                       width = input$width_geneboxplot_rnaseq_raw,
                       height = input$height_geneboxplot_rnaseq_raw,
                       units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_geneboxplot_rnaseq_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_geneboxplot_rnaseq_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_geneboxplot_rnaseq_raw", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_geneboxplot_rnaseq_raw', 
                                          'Download')
                           )
                  )
                )
              ))
            })
            
            #********************************************************************#
            # Output 2: Boxplots
            #********************************************************************#
            
            # Boxplots of all genes together
            output$boxplots_rnaseq_raw <- renderImage({
              req(session$clientData$output_boxplots_rnaseq_raw_width)
              req(session$clientData$output_boxplots_rnaseq_raw_height)
              getBoxplots(experimentFactor = rv$experimentFactor,
                          normData = rv$normData,
                          RNASeq = TRUE,
                          width = session$clientData$output_boxplots_rnaseq_raw_width,
                          height = session$clientData$output_boxplots_rnaseq_raw_height)
            }, deleteFile = TRUE)
            
            
            #***************************#
            # Modal to download figure
            #***************************#
            
            # Download plot
            output$realdownload_boxplots_rnaseq_raw <- downloadHandler(
              filename = function(){"QC_Boxplots.png"},
              content = function(file){
                png(file,
                    width=input$width_boxplots_rnaseq_raw,
                    height=input$height_boxplots_rnaseq_raw,
                    pointsize=24)
                getBoxplots_download(experimentFactor = rv$experimentFactor,
                                     normData = rv$normData,
                                     RNASeq = TRUE)
                dev.off()
              }
            )
            
            
            # Make modal
            observeEvent(input$download_boxplots_rnaseq_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_boxplots_rnaseq_raw", 
                                       "Height",
                                       min = 1200, max = 1600,
                                       value = 1500, step = 1,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_boxplots_rnaseq_raw", 
                                       "Width",
                                       min = 800, max = 1200,
                                       value = 1000, step = 1,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_boxplots_rnaseq_raw', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Output 3: Density plots
            #********************************************************************#
            
            # Densityplots of all genes together
            output$densityplots_rnaseq_raw <- plotly::renderPlotly({
              getDensityplots(experimentFactor = rv$experimentFactor,
                              normMatrix = rv$normData,
                              RNASeq = TRUE)
              
            })
            
            #********************************************************************#
            # Output 4: Heatmap
            #********************************************************************#
            
            # Heatmap of sample-sample correlations
            output$heatmap_rnaseq_raw  <- plotly::renderPlotly({
              
              # Make colors
              if(length(input$heatmapFactor_rnaseq_raw) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$heatmapFactor_rnaseq_raw], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$heatmapFactor_rnaseq_raw])
              }
              
              getHeatmap(experimentFactor = colorFactor,
                         normMatrix = rv$normData,
                         clusterOption1 = input$clusteroption1_rnaseq_raw,
                         clusterOption2 = input$clusteroption2_rnaseq_raw,
                         theme = input$heatmaptheme_rnaseq_raw)
            })
            
            
            #********************************************************************#
            # Output 5: PCA
            #********************************************************************#
            
            #Perform PCA
            rv$PCA_data <- prcomp(t(rv$normData[apply(rv$normData, 1, var) != 0,]),
                                  retx = TRUE, 
                                  center = TRUE,
                                  scale.= TRUE)
            
            
            # Make PCA plot
            output$PCA_rnaseq_raw <- plotly::renderPlotly({
              
              # Make colors
              if(length(input$colorFactor_rnaseq_raw) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_rnaseq_raw], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$colorFactor_rnaseq_raw])
              }
              
              # Make PCA score plot
              plot_PCA(PC_data = rv$PCA_data, 
                       colorFactor = colorFactor, 
                       xpc = as.numeric(stringr::str_remove(input$xpca_rnaseq_raw,"PC")), 
                       ypc = as.numeric(stringr::str_remove(input$ypca_rnaseq_raw,"PC")), 
                       zpc = ifelse(input$xyz_rnaseq_raw,as.numeric(stringr::str_remove(input$zpca_rnaseq_raw,"PC")),3), 
                       xyz = input$xyz_rnaseq_raw)
              
            })
            
            #********************************************************************#
            # Output 6: Overview of pre-processing settings
            #********************************************************************#
            # Print table with settings
            output$processingSettings_rnaseq_raw <- DT::renderDataTable({
              return(rv$processingSettings)
            },options = list(pageLength = 10),
            selection = "none")
            
            # Download button
            output$downloadProcessingSettings_rnaseq_raw <- downloadHandler(
              filename = "preprocessingSettings.csv",
              content = function(file){
                write.csv(rv$processingSettings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
            #======================#
            # UI
            #======================#
            
            # Render UI for main tab
            output$UI_QC_rnaseq_raw <- renderUI({
              tagList(
                tabsetPanel(
                  
                  # TAB1: expression values
                  tabPanel("Expression values",
                           icon = icon("fas fa-mouse-pointer"),
                           
                           # Title + description
                           h3(strong("Normalized expression values")),
                           h5("Here you can view the normalized log intensity (expression) values."),
                           hr(),
                           
                           # Table
                           DT::dataTableOutput(outputId = "exprTable_rnaseq_raw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1"),
                           
                           # Download button
                           downloadButton("downloadNormalizedData_rnaseq_raw", 
                                          "Download table"),
                           
                           br(),
                           br(),
                           
                           # Dropdown Button to adjust the plot settings
                           shinyWidgets::dropdownButton(
                             tags$div(
                               style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                               
                               
                               # Change order of the boxplots:
                               tags$h4("Drag to change boxplot order"),
                               shinyjqui::orderInput(inputId = 'geneboxplot_order_rnaseq_raw', 
                                                     label = NULL, 
                                                     items = levels(rv$experimentFactor),
                                                     item_class = 'default'),
                               br(),
                               
                               # Change colour of the boxplots by button. 
                               # This is used when there are more than 6 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_rnaseq_raw > 6",
                                 tags$h4("Click to change boxplot colours"),
                                 shinyWidgets::actionBttn("geneboxplot_changeOrder_rnaseq_raw",
                                                          label = "Change color",
                                                          style = "simple",
                                                          color = "primary",
                                                          icon = icon("sync"))
                               ),
                               
                               # Change colour of the boxplots by colour picker
                               # This is used when there are less than 7 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_rnaseq_raw < 7",
                                 tags$h4("Click to select boxplot colours"),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 0",
                                   colourpicker::colourInput("geneboxplot_col1_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[1])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 1",
                                   colourpicker::colourInput("geneboxplot_col2_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[2])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 2",
                                   colourpicker::colourInput("geneboxplot_col3_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[3])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 3",
                                   colourpicker::colourInput("geneboxplot_col4_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[4])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 4",
                                   colourpicker::colourInput("geneboxplot_col5_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[5])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_raw > 5",
                                   colourpicker::colourInput("geneboxplot_col6_rnaseq_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[6])
                                 )
                               ),
                               br(),
                               tags$h4("Drag to change jitter"),
                               sliderInput("jitter_geneboxplot_rnaseq_raw", 
                                           NULL,
                                           min = 0, max = 0.3,
                                           value = 0.1, step = 0.01),
                               br()
                             ),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"),
                             
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to personalize plot!")
                             
                           ), # EO dropdownButton
                           
                           # Boxplot of the selected gene's expression values
                           plotOutput("ExprBoxplot_rnaseq_raw")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1"),
                           
                           actionButton("download_geneboxplot_rnaseq_raw", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           br(),
                           br()
                  ),
                  
                  # TAB2: Boxplots of all gene's expression values
                  tabPanel("Boxplots",
                           icon = icon("fas fa-file"),
                           br(),
                           actionButton("download_boxplots_rnaseq_raw", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           plotOutput(outputId = "boxplots_rnaseq_raw",
                                      width = "65vw", height = "80vw")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # TAB3: Density plots of all gene's expression values
                  tabPanel("Density plots",
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           h2(strong("Density plot of normalized counts"), align = "center"),
                           h4("Distributions should be comparable between samples", align = "center"),
                           plotly::plotlyOutput(outputId = "densityplots_rnaseq_raw",
                                                width = "65vw", height = "40vw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # TAB4: Heatmap of sample-sample correlations
                  tabPanel("Correlation plot", 
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           fluidRow(
                             # Select istance method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption1_rnaseq_raw",
                                      label = "Distance calculation 
                                              method",
                                      choices = c("Pearson","Spearman",
                                                  "Euclidean"),
                                      options = list(
                                        style = "btn-primary"))
                             ),
                             
                             #Clustering method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption2_rnaseq_raw",
                                      label = "Clustering method",
                                      choices = c("ward.D2","single",
                                                  "complete","average",
                                                  "mcquitty","median",
                                                  "centroid"),
                                      options = list(
                                        style = "btn-info"))
                             )),
                           
                           hr(),
                           
                           #Theme
                           shinyWidgets::dropdownButton(
                             selectInput(inputId = "heatmapFactor_rnaseq_raw",
                                                       label = "Side colors",
                                                       choices = colnames(rv$metaData_fil),
                                                       selected = rv$experimentName,
                                                       multiple = TRUE),
                             
                             selectInput(inputId = 'heatmaptheme_rnaseq_raw',
                                         label = "Heatmap theme",
                                         choices = c("Default", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"), width = "300px",
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to change colors!")
                           ),
                           
                           plotly::plotlyOutput("heatmap_rnaseq_raw", 
                                                width = "1000px", 
                                                height="600px") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1", 
                                                          proxy.height = "400px")
                           
                  ),
                  
                  # TAB5: PCA score plot
                  tabPanel("PCA",
                           icon = icon("fas fa-mouse-pointer"),
                           # Title + text
                           fluidRow(
                             h3(strong("Principal Component Analysis (PCA)")),
                             h5("Here you can view the PCA score plot."),
                             hr(),
                           ),
                           
                           # Set color + 3D/2D
                           fluidRow(
                             column(3,
                                    # Color by which factor?
                                    shinyWidgets::pickerInput(inputId = "colorFactor_rnaseq_raw",
                                                              label = "Color by:",
                                                              choices = colnames(rv$metaData_fil),
                                                              selected = rv$experimentName,
                                                              multiple = TRUE)
                             ),
                             column(3,
                                    br(),
                                    # 3D plot?
                                    shinyWidgets::materialSwitch(
                                      inputId = "xyz_rnaseq_raw",
                                      label = "3D",
                                      value = FALSE, 
                                      status = "danger")
                                    
                             )
                           ),
                           
                           # Set axes
                           fluidRow(
                             column(3,
                                    #X-axis
                                    selectInput(inputId = "xpca_rnaseq_raw", 
                                                label = "x-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC1")
                             ),
                             column(3,
                                    #Y-axis
                                    selectInput(inputId = "ypca_rnaseq_raw", 
                                                label = "y-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC2")
                             ),
                             column(3,
                                    #Z-axis
                                    conditionalPanel(
                                      condition = "input.xyz_rnaseq_raw==true",
                                      selectInput(inputId = "zpca_rnaseq_raw", 
                                                  label = "z-axis",
                                                  choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                              "PC6", "PC7", "PC8"),
                                                  selected = "PC3")
                                    )
                             )
                           ),
                           
                           # Print plot
                           fluidRow(
                             hr(),
                             plotly::plotlyOutput("PCA_rnaseq_raw")%>% 
                               shinycssloaders::withSpinner(color="#0dc5c1"),
                           ),

                           
                  ), # EO PCA tabpanel
                  
                  # TAB6: Settings table
                  tabPanel("Settings overview",
                           icon = icon("fas fa-file"),
                           h3(strong("Pre-processing settings")),
                           h5("Here you can see an overview of the chosen pre-processing settings."),
                           hr(),
                           DT::dataTableOutput(outputId = "processingSettings_rnaseq_raw") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadProcessingSettings_rnaseq_raw", 
                                          "Download"),
                           
                  ) # EO Settings tabPanel
                ) # EO tabsetPanel
              ) # EO tagList
            }) # EO renderUI
            
            # Allow user to go to next tab
            output$UI_next_preprocessing_rnaseq_raw <- renderUI({
              req(rv$normData)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                shinyWidgets::actionBttn(inputId = "next_preprocessing_rnaseq_raw",
                                         label = "Next",
                                         style = "jelly",
                                         color = "danger",
                                         icon = icon("arrow-right"))
              )
            })
            
          }) # EO observeEvent
          
          
          
          #======================================================================#
          # TAB4: Statistical analysis
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_preprocessing_rnaseq_raw,{
            
            # Go to statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_statistics_rnaseq_raw")
          })
          
          # SIDEPANEL:
          
          # Select continuous covariates
          output$UI_covGroups_num_rnaseq_raw <- renderUI({
            tagList(
              shinyWidgets::pickerInput(inputId = "covGroups_num_rnaseq_raw",
                                        label = "Continuous covariates (e.g., age):",
                                        choices = setdiff(colnames(rv$metaData_fil),
                                                          rv$experimentName),
                                        selected = NULL,
                                        multiple = TRUE)
            )
          })
          
          # Select discrete covariates
          output$UI_covGroups_char_rnaseq_raw <- renderUI({
            tagList(
              shinyWidgets::pickerInput(inputId = "covGroups_char_rnaseq_raw",
                                        label = "Discrete covariates (e.g., sex):",
                                        choices = setdiff(colnames(rv$metaData_fil),
                                                          rv$experimentName),
                                        selected = NULL,
                                        multiple = TRUE)
            )
          })
          
          # Select comparisons
          output$UI_comparisons_rnaseq_raw <- renderUI({
            
            tagList(
              shinyWidgets::multiInput(
                inputId = "comparisons_rnaseq_raw",
                label = "Comparisons:", 
                choices = makeComparisons(make.names(unique(rv$experimentFactor))),
                selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1]
              )
            )
          })
          
          
          # Select comparisons
          output$UI_biomart_dataset_rnaseq_raw <- renderUI({
            req(input$addAnnotation_rnaseq_raw)
            shinyWidgets::pickerInput(inputId = "biomart_dataset_rnaseq_raw",
                                      label = "Organism:",
                                      choices = c("Homo sapiens" = "hsapiens_gene_ensembl" ,
                                                  "Bos taurus" = "btaurus_gene_ensembl",
                                                  "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                                  "Mus musculus" = "mmusculus_gene_ensembl",
                                                  "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                                      selected = "hsapiens_gene_ensembl",
                                      multiple = FALSE)
          })
          
          output$UI_addAnnotations_rnaseq_raw <- renderUI({
            req(input$addAnnotation_rnaseq_raw)
            req(input$biomart_dataset_rnaseq_raw)
            
            tagList(
              
              shinyWidgets::pickerInput(inputId = "biomart_filter_rnaseq_raw",
                                        label = "Gene ID",
                                        choices = c("Ensembl Gene ID",
                                                    "Entrez Gene ID",
                                                    "Gene Symbol/Name"),
                                        selected = "Entrez Gene ID",
                                        multiple = FALSE),
              
              shinyWidgets::pickerInput(inputId = "biomart_attributes_rnaseq_raw",
                                        label = "Output",
                                        choices = c("Ensembl Gene ID",
                                                    "Entrez Gene ID",
                                                    "Gene Symbol/Name"),
                                        selected = "Gene Symbol/Name",
                                        multiple = TRUE)
            )
            
          })
          
          
          
          #=========================================#
          # Generate output of statistical analysis
          #=========================================#
          observeEvent(input$calculate_statistics_rnaseq_raw,{
            shinybusy::show_modal_spinner(text = "Statistical analysis...",
                                          color="#0dc5c1")
            
            # Calculate statistics
            if (isTRUE(input$addAnnotation_rnaseq_raw)){
              rv$top_table_list <- getStatistics_RNASeq(rawMatrix = rv$gxData_fil, 
                                                   metaData = rv$metaData_fil, 
                                                   expFactor = rv$experimentName,
                                                   covGroups_num = input$covGroups_num_rnaseq_raw,
                                                   covGroups_char = input$covGroups_char_rnaseq_raw,
                                                   comparisons = input$comparisons_rnaseq_raw,
                                                   filterThres = rv$filterThreshold,
                                                   smallestGroupSize = length(unique(rv$experimentFactor)),
                                                   addAnnotation = input$addAnnotation_rnaseq_raw,
                                                   biomart_dataset = input$biomart_dataset_rnaseq_raw,
                                                   biomart_attributes = unique(c(input$biomart_filter_rnaseq_raw,
                                                                                 input$biomart_attributes_rnaseq_raw)),
                                                   biomart_filters = input$biomart_filter_rnaseq_raw)
              
            } else{
              rv$top_table_list <- getStatistics_RNASeq(rawMatrix = rv$gxData_fil, 
                                                   metaData = rv$metaData_fil, 
                                                   expFactor = rv$experimentName, 
                                                   covGroups_num = input$covGroups_num_rnaseq_raw,
                                                   covGroups_char = input$covGroups_char_rnaseq_raw,
                                                   comparisons = input$comparisons_rnaseq_raw,
                                                   filterThres = rv$filterThreshold,
                                                   smallestGroupSize = length(unique(rv$experimentFactor)),
                                                   addAnnotation = input$addAnnotation_rnaseq_raw,
                                                   biomart_dataset = NULL,
                                                   biomart_attributes = NULL,
                                                   biomart_filters = NULL)
              
            }
            
            rv$newFactor <- rv$experimentFactor
            rv$newData <- rv$normData
            
            # Select comparison for output
            observe({
              
              if (!is.null(rv$top_table_list)){
                rv$top_table <- rv$top_table_list[[1]]
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show comparisons
                output$UI_comparisons_view_rnaseq_raw <- renderUI({
                  shinyWidgets::pickerInput(inputId = "comparisons_view_rnaseq_raw",
                                            label = "Select comparison:",
                                            choices = names(rv$top_table),
                                            selected = names(rv$top_table)[1],
                                            multiple = FALSE)
                })
                
                # Show message
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = rv$top_table_list[[2]],
                  type = "info")
                
                # Show microarray statistics tab
                showTab("navbar", target = "panel_ORA_rnaseq_raw")
              } else{
                shinybusy::remove_modal_spinner()
                # Show message
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = "Oops...something went wrong! Please try again!",
                  type = "error")
              }
              
            })
            
            
            #=========================================#
            # Plot output of statistical analysis
            #=========================================#
            
            #********************************************************************#
            # TAB1: top table
            #********************************************************************#
            
            observe({
              req(input$comparisons_view_rnaseq_raw)
              req(rv$top_table)
              
              # print top table
              output$top_table_rnaseq_raw <- DT::renderDataTable({
                
                if (is.null(input$comparisons_view_rnaseq_raw)){
                  output <- rv$top_table[[1]]
                } else{
                  output <- rv$top_table[[input$comparisons_view_rnaseq_raw]]
                }
                
                return(output)
              },options = list(pageLength = 6),
              selection = list(mode = "single", selected = 1), escape = FALSE)
              
              # Download button
              output$download_top_table_rnaseq_raw <- downloadHandler(
                filename = paste0("topTable_",input$comparisons_view_rnaseq_raw,".csv"),
                content = function(file){
                  write.csv(rv$top_table[[input$comparisons_view_rnaseq_raw]], file, quote = FALSE, row.names = FALSE)
                }
              )
            })
            
            # Change plotting data depending on whether all experimental groups will be plotted  
            observe({
              req(input$comparisons_view_rnaseq_raw)
              req(rv$top_table)
              
              if (!is.null(input$boxplotAll_rnaseq_raw)){
                if (input$boxplotAll_rnaseq_raw){
                  rv$newFactor <- rv$experimentFactor
                  rv$newData <- rv$normData
                }
                if (!input$boxplotAll_rnaseq_raw){
                  if(length(rv$experimentName) > 1){
                    t <- make.names(apply(rv$metaData_fil[,rv$experimentName], 1, paste, collapse = "_" ))
                  } else{
                    t <- make.names(rv$metaData_fil[,rv$experimentName])
                  }
                  rv$newData <- rv$normData[,t %in% (unlist(stringr::str_split(input$comparisons_view_rnaseq_raw, " - ")))]
                  rv$newFactor <- factor(as.character(rv$experimentFactor)[t %in% (unlist(stringr::str_split(input$comparisons_view_rnaseq_raw, " - ")))])
                }
              }

              # Change color by click on button
              rv$colorOrder <- 1:length(levels(rv$newFactor))
              observeEvent(input$statboxplot_changeOrder_rnaseq_raw,{
                all_orders <- permute(1:length(levels(rv$newFactor))) 
                sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
                
                if (sel < length(all_orders)){
                  rv$colorOrder <- all_orders[[sel+1]]
                } else{
                  rv$colorOrder <- all_orders[[1]]
                }
              })
              
              
              # Boxplot of single gene (based on selected row in the top table)
              output$ExprBoxplot_statistics_rnaseq_raw <- renderPlot({
                req(input$top_table_rnaseq_raw_rows_selected)
                req(rv$top_table)
                
                if (length(levels(rv$newFactor)) > 6){
                  legendColors <- colorsByFactor(rv$newFactor)$legendColors
                } else{
                  legendColors <- c(input$statboxplot_col1_rnaseq_raw,
                                    input$statboxplot_col2_rnaseq_raw,
                                    input$statboxplot_col3_rnaseq_raw,
                                    input$statboxplot_col4_rnaseq_raw,
                                    input$statboxplot_col5_rnaseq_raw,
                                    input$statboxplot_col6_rnaseq_raw)
                }
                
                gene <- rv$top_table[[input$comparisons_view_rnaseq_raw]]$GeneID[input$top_table_rnaseq_raw_rows_selected]
                sel_row <- which(as.character(rownames(rv$normData)) %in% as.character(gene))
                
                # Make boxplot
                rv$temp <- geneBoxplot(experimentFactor = rv$newFactor, 
                            normMatrix = rv$newData, 
                            sel_row = sel_row,
                            legendColors = legendColors[rv$colorOrder],
                            groupOrder = input$statboxplot_order_rnaseq_raw,
                            jitter = input$jitter_statboxplot_rnaseq_raw,
                            rnaseq=TRUE,
                            seed = sample(1:1000,1))
                
                return(rv$temp)
                
              })
              
              # Get number of experimental groups
              output$length_statboxplot_rnaseq_raw <- reactive({
                length(levels(rv$experimentFactor))
              })
              outputOptions(output, "length_statboxplot_rnaseq_raw", suspendWhenHidden = FALSE) 
              
            })
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_statboxplot_rnaseq_raw <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp, 
                                filename = file,
                                width = input$width_statboxplot_rnaseq_raw,
                                height = input$height_statboxplot_rnaseq_raw,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_statboxplot_rnaseq_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_statboxplot_rnaseq_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_statboxplot_rnaseq_raw", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_statboxplot_rnaseq_raw', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # TAB2: P value and logFC histograms
            #********************************************************************#
            
            observe({
              req(input$comparisons_view_rnaseq_raw)
              req(rv$top_table)
              
              if (input$comparisons_view_rnaseq_raw %in% names(rv$top_table)){
                # P value histogram
                output$Phistogram_rnaseq_raw <- plotly::renderPlotly({
                  req(rv$top_table)
                  p <- makePHistogram(rv$top_table[[input$comparisons_view_rnaseq_raw]][,"p-value"])
                  return(p)
                })
                
                # logFC histrogram
                output$logFChistogram_rnaseq_raw <- plotly::renderPlotly({
                  req(rv$top_table)
                  p <- makelogFCHistogram(rv$top_table[[input$comparisons_view_rnaseq_raw]][,"log2FC"])
                  return(p)
                })
              }
            })
            
            #********************************************************************#
            # TAB3: Volcano plot
            #********************************************************************#
            
            observeEvent(input$plot_volcano_rnaseq_raw, {
              req(rv$top_table)
              req(input$rawp_volcano_rnaseq_raw)
              req(input$p_thres_volcano_rnaseq_raw)
              req(input$logFC_thres_volcano_rnaseq_raw)
              req(input$comparisons_view_rnaseq_raw)
              
              if (input$comparisons_view_rnaseq_raw %in% names(rv$top_table)){
                p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_rnaseq_raw]], 
                                 p = input$rawp_volcano_rnaseq_raw, 
                                 p_threshold = input$p_thres_volcano_rnaseq_raw, 
                                 logFC_threshold = input$logFC_thres_volcano_rnaseq_raw)
                
                output$volcano_rnaseq_raw <- plotly::renderPlotly(p)
              }
            }, ignoreNULL = FALSE) 
            # ignoreNULL: generate plot even if action button is not pressed
            
            
            #=========================================#
            #  # UI: Output in different tabs
            #=========================================#
            
            observe({
              req(rv$top_table)
              
              output$UI_boxplotAll_rnaseq_raw <- renderUI({
                tagList(
                  # Change order of the boxplots:
                  tags$h4("Drag to change boxplot order"),
                  shinyjqui::orderInput(inputId = 'statboxplot_order_rnaseq_raw', 
                                        label = NULL, 
                                        items = levels(rv$newFactor),
                                        item_class = 'default'),
                  br(),
                  
                  # Change colour of the boxplots by button. 
                  # This is used when there are more than 6 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_rnaseq_raw > 6",
                    tags$h4("Click to change boxplot colours"),
                    shinyWidgets::actionBttn("statboxplot_changeOrder_rnaseq_raw",
                                             label = "Change color",
                                             style = "simple",
                                             color = "primary",
                                             icon = icon("sync"))
                  ),
                  
                  # Change colour of the boxplots by colour picker
                  # This is used when there are less than 7 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_rnaseq_raw < 7",
                    tags$h4("Click to select boxplot colours"),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 0",
                      colourpicker::colourInput("statboxplot_col1_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[1])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 1",
                      colourpicker::colourInput("statboxplot_col2_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[2])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 2",
                      colourpicker::colourInput("statboxplot_col3_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[3])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 3",
                      colourpicker::colourInput("statboxplot_col4_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[4])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 4",
                      colourpicker::colourInput("statboxplot_col5_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[5])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_raw > 5",
                      colourpicker::colourInput("statboxplot_col6_rnaseq_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[6])
                    )
                  ),
                  br(),
                  tags$h4("Drag to change jitter"),
                  sliderInput("jitter_statboxplot_rnaseq_raw", 
                              NULL,
                              min = 0, max = 0.3,
                              value = 0.1, step = 0.01),
                  br()
                )
              })
            })
            
            observe({
              if (is.null(rv$top_table)){
                output$UI_output_statistics_rnaseq_raw <- renderUI(NULL)
              } else{

              output$UI_output_statistics_rnaseq_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    
                    #********************************************************************#
                    # top table tab
                    #********************************************************************#
                    
                    tabPanel("Top table",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             h3(strong("Top Table")),
                             h5("The Top Table includes the output of the selected statistical analysis."),
                             hr(),
                             DT::dataTableOutput(outputId = "top_table_rnaseq_raw") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1"),
                             downloadButton("download_top_table_rnaseq_raw", 
                                            "Download table"),
                             br(),
                             br(),
                             
                             # Dropdown Button to adjust the plot settings
                             shinyWidgets::dropdownButton(
                               tags$div(
                                 style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                                 
                                 # Plot all experimental groups?
                                 conditionalPanel(
                                   condition = "output.length_statboxplot_rnaseq_raw > 2",
                                   tags$h4("Plot all experimental groups?"),
                                   shinyWidgets::materialSwitch(inputId = "boxplotAll_rnaseq_raw",
                                                                label = NULL, 
                                                                value = TRUE,
                                                                status = "primary"),
                                   
                                   br()
                                 ),
                                 uiOutput("UI_boxplotAll_rnaseq_raw")
                                 
                                 
                                 
                               ),
                               circle = TRUE, status = "info",
                               icon = icon("fas fa-cog"),
                               
                               tooltip = shinyWidgets::tooltipOptions(
                                 title = "Click to personalize plot!")
                               
                             ), # EO dropdownButton
                             plotOutput("ExprBoxplot_statistics_rnaseq_raw")%>% 
                               shinycssloaders::withSpinner(color="#0dc5c1"),
                             actionButton("download_statboxplot_rnaseq_raw", 
                                          "Download figure",
                                          icon = shiny::icon("download")),
                             br(),
                             br()
                    ),
                    #********************************************************************#
                    # histogram tab
                    #********************************************************************#
                    
                    tabPanel("Histograms",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             plotly::plotlyOutput("Phistogram_rnaseq_raw")%>% 
                               shinycssloaders::withSpinner(color="#0dc5c1"),
                             hr(),
                             plotly::plotlyOutput("logFChistogram_rnaseq_raw")%>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")
                             
                    ),
                    
                    #********************************************************************#
                    # volcano tab
                    #********************************************************************#
                    
                    tabPanel("Volcano plot",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             
                             fluidRow(
                               column(3,
                                      # Raw or adjusted P value?
                                      prettyRadioButtons(
                                        inputId = "rawp_volcano_rnaseq_raw",
                                        label = "P value", 
                                        choices = 
                                          c("Raw P value" = "raw", 
                                            "Adjusted P value" = "adj"))
                               ),
                               column(3,
                                      #P value threshold
                                      numericInput(
                                        inputId = "p_thres_volcano_rnaseq_raw",
                                        label = "P threshold",
                                        value = 0.05)
                               ),
                               column(3,
                                      #logFC threshold
                                      numericInput(
                                        inputId = "logFC_thres_volcano_rnaseq_raw",
                                        label = "logFC threshold",
                                        value = 1)
                               )
                             ),
                             hr(),
                             
                             # Button to reload the volcano plot
                             shinyWidgets::actionBttn(inputId = "plot_volcano_rnaseq_raw", 
                                                      label = "Plot",
                                                      style = "jelly",
                                                      color = "primary",
                                                      icon = icon("sync")),
                             br(),
                             br(),
                             # Volcano plot output
                             plotly::plotlyOutput("volcano_rnaseq_raw")%>% 
                               withSpinner(color="#0dc5c1")
                             
                    )
                    
                  ) # tabsetpanel
                ) # taglist
              }) # renderUI
              }
            }) # Observe
            
            # Allow user to go to next tab
            output$UI_next_statistics_rnaseq_raw <- renderUI({
              req(rv$top_table)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                shinyWidgets::actionBttn(inputId = "next_statistics_rnaseq_raw",
                                         label = "Next",
                                         style = "jelly",
                                         color = "danger",
                                         icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          
          #======================================================================#
          # TAB5: ORA
          #======================================================================#
          # Go to data upload tab
          observeEvent(input$next_statistics_rnaseq_raw,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_ORA_rnaseq_raw")
          })
          
          #********************************************************************#
          # Options for ORA
          #********************************************************************#
          
          # Comparisons for which ORA should be performed
          observe({
            req(rv$top_table)
            output$UI_comparisons_view_ORA_rnaseq_raw <- renderUI({
              shinyWidgets::pickerInput(inputId = "comparisons_view_ORA_rnaseq_raw",
                                        label = NULL,
                                        choices = names(rv$top_table),
                                        selected = names(rv$top_table)[1],
                                        multiple = FALSE)
            })
          })
          
          observe({
            
            # Get all columns of the top table that can possibly contain the gene IDs
            req(input$comparisons_view_ORA_rnaseq_raw)
            col_choice <- 1
            if (length(colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]])) > 6){
              col_choice <- c(1,7:ncol(rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]]))
            }
            
            # Several options need to be selected before ORA can be performed
            output$UI_geneID_ORA_rnaseq_raw <- renderUI({
              tagList(
                
                # Select organism
                selectInput(inputId = "organism_ORA_rnaseq_raw",
                            label = "Organism",
                            choices = c("Bos taurus",
                                        "Caenorhabditis elegans",
                                        "Homo sapiens",
                                        "Mus musculus", 
                                        "Rattus norvegicus"),
                            selected = "Homo sapiens"),
                
                # Which columns of the top table contains the gene ids?
                shinyWidgets::pickerInput(inputId = "geneID_ORA_rnaseq_raw",
                                          label = "Which column of the top table contains the gene IDs?",
                                          choices = colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]])[col_choice],
                                          selected = colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]])[1],
                                          multiple = FALSE),
                
                # Which gene IDs do they column contain?
                shinyWidgets::pickerInput(inputId = "selID_ORA_rnaseq_raw",
                                          label = "Which gene ID to use?",
                                          choices = c("Ensembl Gene ID" = "ENSEMBL", 
                                                      "Entrez Gene ID" = "ENTREZID", 
                                                      "Gene Symbol/Name" = "SYMBOL"),
                                          selected = "ENTREZID",
                                          multiple = FALSE),
              )
            })
          })
          
          #********************************************************************#
          # Perform ORA
          #********************************************************************#
          
          observeEvent(input$calculate_ORA_rnaseq_raw,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Overrepresentation analysis...",
                                          color="#0dc5c1")
            
            # Perform ORA:
            
            # Perform ORA based on logFC/P value threshold(s)
            if (input$topNorThres_rnaseq_raw == "Threshold"){
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]],
                                 geneset = input$geneset_ORA_rnaseq_raw,
                                 geneID_col = input$geneID_ORA_rnaseq_raw,
                                 geneID_type = input$selID_ORA_rnaseq_raw,
                                 organism = input$organism_ORA_rnaseq_raw,
                                 updown = input$updown_ORA_rnaseq_raw,
                                 topN = FALSE,
                                 N = NULL,
                                 rawadj = input$rawp_ORA_rnaseq_raw,
                                 p_thres = input$p_thres_ORA_rnaseq_raw,
                                 logFC_thres = input$logFC_thres_ORA_rnaseq_raw)
              
              # Perform ORA based on top N most significant genes
            } else{
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]],
                                 geneset = input$geneset_ORA_rnaseq_raw,
                                 geneID_col = input$geneID_ORA_rnaseq_raw,
                                 geneID_type = input$selID_ORA_rnaseq_raw,
                                 organism = input$organism_ORA_rnaseq_raw,
                                 updown = input$updown_ORA_rnaseq_raw,
                                 topN = TRUE,
                                 N = input$topN_rnaseq_raw,
                                 rawadj = NULL,
                                 p_thres = NULL,
                                 logFC_thres = NULL)
            }
            
            
            #********************************************************************#
            # Generate ORA output
            #********************************************************************#
            
            observe({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show error message
              if (is.null(rv$ORA_data)){
                sendSweetAlert(
                  session = session,
                  title = "Error!",
                  text = "No significant genes!",
                  type = "error")
                
                # Show success message
              }else{
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = "Gene overrepresentation analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
                  type = "info")
                
                #--------------------------------------------------------------#
                # ORA statistics table
                #--------------------------------------------------------------#
                
                output$ORA_table_rnaseq_raw <- DT::renderDataTable({
                  req(input$geneset_ORA_rnaseq_raw)
                  output <- rv$ORA_data@result
                  
                  # Link to wikipathways website if geneset ID is from WikiPathways
                  if (input$geneset_ORA_rnaseq_raw == "WikiPathways"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.wikipathways.org/pathways/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  # Link to QuickGO website if geneset ID is a GO term
                  if (input$geneset_ORA_rnaseq_raw == "KEGG"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.genome.jp/pathway/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_rnaseq_raw %in% c("GO-BP", "GO-MF", "GO-CC")){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.ebi.ac.uk/QuickGO/term/",
                        output$ID
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  
                  return(output)
                },options = list(pageLength = 6),
                selection = list(mode = "single", selected = 1), escape = FALSE)
                
                # Download button
                output$download_ORA_table_rnaseq_raw <- downloadHandler(
                  filename = paste0("ORATable_",input$comparisons_view_ORA_rnaseq_raw,"_",input$geneset_ORA_rnaseq_raw,".csv"),
                  content = function(file){
                    write.csv(rv$ORA_data@result, file, quote = FALSE, row.names = FALSE)
                  }
                )
                
                # Print statistics of genes in selected Term
                output$ORAgene_table_rnaseq_raw <- DT::renderDataTable({
                  req(input$ORA_table_rnaseq_raw_rows_selected)
                  
                  # Make ORA gene table
                  output <- make_ORAgene_table(ORA_data = rv$ORA_data,
                                               top_table = rv$top_table[[input$comparisons_view_ORA_rnaseq_raw]],
                                               geneID_col = input$geneID_ORA_rnaseq_raw,
                                               sel_row_ORA = input$ORA_table_rnaseq_raw_rows_selected)
                  
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  output$meanExpr <- round(output$meanExpr,3)
                  output$log2FC <- round(output$log2FC,3)
                  output$`log2FC SE` <- round(output$`log2FC SE`,3)
                  
                  return(output)
                }, options = list(pageLength = 6), escape = FALSE)
                
                # Text for gene table
                output$text_ORAgene_table_rnaseq_raw <- renderText({
                  text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_rnaseq_raw_rows_selected,"ID"],
                                 "</b></h3>")
                  return(text)
                })
                
                #--------------------------------------------------------------#
                # ORA barchart
                #--------------------------------------------------------------#
                
                observeEvent(input$plot_ORAplot_rnaseq_raw, {
                  req(input$nSets_ORAplot_rnaseq_raw)
                  p <- makeORAplot(rv$ORA_data,
                                   nSets = input$nSets_ORAplot_rnaseq_raw)
                  
                  output$ORAplot_rnaseq_raw <- plotly::renderPlotly(p)
                  
                }, ignoreNULL = FALSE)
                
                #--------------------------------------------------------------#
                # ORA network diagram
                #--------------------------------------------------------------#
                
                observeEvent(input$plot_ORAnetwork_rnaseq_raw, {
                  req(input$layout_network_rnaseq_raw)
                  req(input$nSets_network_rnaseq_raw)
                  p <- makeORAnetwork(ORA_data = rv$ORA_data,
                                      layout = input$layout_network_rnaseq_raw,
                                      nSets = input$nSets_network_rnaseq_raw)
                  
                  output$ORAnetwork_rnaseq_raw <- renderPlot(p,
                                                             height = 500, 
                                                             width = 800)
                  
                }, ignoreNULL = FALSE)
                
              } # EO ifelse
            }) # EO observe
            
            
            
            #********************************************************************#
            # UI: ORA output
            #********************************************************************#
            observe({
              if (!is.null(rv$ORA_data)){
                output$UI_output_ORA_rnaseq_raw <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #--------------------------------------------------------------#
                      # ORA statistics table
                      #--------------------------------------------------------------#
                      tabPanel("ORA table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               # Title + description of statistics table
                               h3(strong("Statistics table")),
                               h5("The ORA statistics table encompasses the output of the gene 
                              overrepresentation analysis."),
                               hr(),
                               
                               # Statistics table
                               DT::dataTableOutput(outputId = "ORA_table_rnaseq_raw") %>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1"),
                               
                               # Download button
                               downloadButton("download_ORA_table_rnaseq_raw", 
                                              "Download"),
                               br(),
                               
                               # Title + description of gene table
                               htmlOutput("text_ORAgene_table_rnaseq_raw"),
                               h5(paste0("The gene table encompasses the statistics of all genes 
                              from the selected geneset.")),
                               hr(),
                               
                               # Gene table
                               DT::dataTableOutput(outputId = "ORAgene_table_rnaseq_raw") %>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1")
                      ),
                      
                      
                      #--------------------------------------------------------------#
                      # ORA bar chart
                      #--------------------------------------------------------------#
                      tabPanel("Bar chart",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               # Title + description of bar chart
                               h3(strong("ORA bar chart")),
                               h5("The ORA bar chart visualizes the results from the overrepresentation analysis."),
                               hr(),
                               
                               # Select number of genesets shown in graph
                               fluidRow(
                                 column(3,
                                        # Number of genesets
                                        numericInput(
                                          inputId = "nSets_ORAplot_rnaseq_raw",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 )
                               ),
                               hr(),
                               
                               # Actionbutton: press to reload plot
                               shinyWidgets::actionBttn(inputId = "plot_ORAplot_rnaseq_raw", 
                                                        label = "Plot",
                                                        style = "jelly",
                                                        color = "primary",
                                                        icon = icon("sync")),
                               br(),
                               br(),
                               
                               # Interactive plot output
                               plotly::plotlyOutput("ORAplot_rnaseq_raw") %>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1")
                      ),
                      
                      #--------------------------------------------------------------#
                      # ORA network diagram
                      #--------------------------------------------------------------#
                      tabPanel("Network diagram",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               # Title + description of the network diagram
                               h3(strong("ORA network diagram")),
                               h5("The ORA network diagram visualize the similarity between the most significant genesets."),
                               hr(),
                               
                               # Select plotting options
                               fluidRow(
                                 column(3,
                                        # Select number of genesets in network
                                        numericInput(
                                          inputId = "nSets_network_rnaseq_raw",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 ),
                                 column(3,
                                        # Select network layout
                                        shinyWidgets::pickerInput(inputId = "layout_network_rnaseq_raw",
                                                                  label = "Network layout",
                                                                  choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                                              'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                                  selected = 'graphopt',
                                                                  multiple = FALSE)
                                 )
                               ),
                               hr(),
                               
                               # Actionbutton: press to reload plot
                               shinyWidgets::actionBttn(inputId = "plot_ORAnetwork_rnaseq_raw", 
                                                        label = "Plot",
                                                        style = "jelly",
                                                        color = "primary",
                                                        icon = icon("sync")),
                               br(),
                               br(),
                               
                               # Make plot
                               plotOutput("ORAnetwork_rnaseq_raw") %>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1")
                      ),
                      
                    ) # EO tabSetPanel
                  ) # EO tagList
                }) # EO renderUI
              }# EO if !is.null(rv$ORA_data)
            }) # EO observe
            
          }) # EO observeEvent
          
          
        } # EO raw or norm
        
        
        
        #************************************************************************#
        # Preprocessed RNA-seq data
        #************************************************************************#
        if (raw_or_norm() == "Processed data"){
          
          #======================================================================#
          # Data Upload
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$startAnalysis,{
            
            # Show RNA-seq (norm) upload tab
            showTab("navbar", target = "panel_upload_rnaseq_norm")
            
            # Remove the other RNA-seq (raw) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_raw")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_raw")
            hideTab("navbar", target = "panel_statistics_rnaseq_raw" )
            hideTab("navbar", target = "panel_ORA_rnaseq_raw")
            
            # Remove the other RNA-seq (norm) tabs
            hideTab("navbar", target = "panel_preprocessing_rnaseq_norm")
            hideTab("navbar", target = "panel_statistics_rnaseq_norm")
            hideTab("navbar", target = "panel_ORA_rnaseq_norm")
            
            # Remove the other microarray (raw) tabs
            hideTab("navbar", target = "panel_upload_microarray_raw")
            hideTab("navbar", target = "panel_preprocessing_microarray_raw")
            hideTab("navbar", target = "panel_statistics_microarray_raw" )
            hideTab("navbar", target = "panel_ORA_microarray_raw")
            
            # Remove the other microarray (norm) tabs
            hideTab("navbar", target = "panel_upload_microarray_norm")
            hideTab("navbar", target = "panel_preprocessing_microarray_norm")
            hideTab("navbar", target = "panel_statistics_microarray_norm")
            hideTab("navbar", target = "panel_ORA_microarray_norm")
            
            # Go to RNA-seq tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_upload_rnaseq_norm")
            
            # Example meta data file
            output$downloadmeta_example_rnaseq_norm <- downloadHandler(
              filename = "MetaData_example.csv",
              content = function(file){
                write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
              }
            )
          })
          
          #----------------------------------------------------------------------#
          # Upload files
          #----------------------------------------------------------------------#
          
          observeEvent(input$upload_rnaseq_norm,{
            
            # Show modal
            show_modal_spinner(text = "Reading data...",
                               color="#0dc5c1")
            
            # Read expression data
            if (!is.null(input$uploadExprData_rnaseq_norm)){
              rv$gxData <- readRNASeq(path=input$uploadExprData_rnaseq_norm$datapath)
            }else{
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
            }
            
            
            # Read metadata
            if (!is.null(input$uploadMeta_rnaseq_norm_tsv)){
              if (input$MetaFileType_rnaseq_norm ==".tsv/.csv file"){
                rv$metaData <- getMetaData(path = input$uploadMeta_rnaseq_norm_tsv$datapath,
                                           celfiles = colnames(rv$gxData),
                                           filetype = input$MetaFileType_rnaseq_norm)
              }
            }else if (!is.null(input$uploadMeta_rnaseq_norm_smf)){
              if (input$MetaFileType_rnaseq_norm =="Series Matrix File"){
                rv$metaData <- getMetaData(path = input$uploadMeta_rnaseq_norm_smf$datapath,
                                           celfiles =  colnames(rv$gxData),
                                           filetype = input$MetaFileType_rnaseq_norm)
              }
            }else{
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
            }
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(rv$gxData)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- rv$gxData[,rownames(rv$metaData)]
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_rnaseq_norm <- renderDataTable(NULL)
              output$metaTable_rnaseq_norm <- renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_rnaseq_norm <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(rv$gxData,6)
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_rnaseq_norm <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              observeEvent(input$metaTable_rnaseq_norm_cell_edit, {
                row  <- input$metaTable_rnaseq_norm_cell_edit$row
                clmn <- input$metaTable_rnaseq_norm_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_rnaseq_norm_cell_edit$value
              })
              
              # Render UI for main tab
              output$UI_upload_rnaseq_norm <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_rnaseq_norm") %>% 
                               withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_rnaseq_norm") %>% 
                               withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_rnaseq_norm <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_rnaseq_norm")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(rv$gxData)){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_rnaseq_norm",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # Observe event
          
          #----------------------------------------------------------------------#
          # Run Example
          #----------------------------------------------------------------------#
          
          observeEvent(input$example_rnaseq_norm,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Read expression data
            rv$gxData <- readRNASeq(path="Data/RNAseq/normExpr_GSE128380.csv")
            
            # Get metadata
            rv$metaData <- getMetaData(path = "Data/RNAseq/sampleInfo_GSE128380.csv",
                                       celfiles = colnames(rv$gxData),
                                       filetype = ".tsv/.csv file")
            
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(rv$gxData)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- rv$gxData[,rownames(rv$metaData)]
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_rnaseq_norm <- DT::renderDataTable(NULL)
              output$metaTable_rnaseq_norm <- DT::renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_rnaseq_norm <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(rv$gxData,6)
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_rnaseq_norm <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              observeEvent(input$metaTable_rnaseq_norm_cell_edit, {
                row  <- input$metaTable_rnaseq_norm_cell_edit$row
                clmn <- input$metaTable_rnaseq_norm_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_rnaseq_norm_cell_edit$value
              })
              
              
              # Render UI for main tab
              output$UI_upload_rnaseq_norm <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_rnaseq_norm") %>% 
                               withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_rnaseq_norm") %>% 
                               withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_rnaseq_norm <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_rnaseq_norm")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(rv$gxData)){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_rnaseq_norm",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # Observe event
          
          
          #======================================================================#
          # Data Pre-processing
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_upload_rnaseq_norm,{
            
            # Go to RNA-seq tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_preprocessing_rnaseq_norm")
          })
          
          
          # 1. Select outliers
          output$UI_outlier_rnaseq_norm <- renderUI({
            if(!input$outlier_rnaseq_norm){
              samples <- rownames(rv$metaData)
              
              pickerInput(inputId = "select_outliers_rnaseq_norm",
                          label = tags$span(
                            "Select samples to be removed", 
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select one or more samples to exclude from the analysis.", 
                                         position = "right",
                                         size = "large")
                          ),
                          choices = samples,
                          multiple = TRUE)
            } else{
              NULL
            }
            
          })
          
          # 2. Select experimental groups
          output$UI_groupselect_rnaseq_norm <- renderUI({
            pickerInput(inputId = "groupselect_rnaseq_norm",
                        label = NULL,
                        choices = colnames(rv$metaData),
                        selected = autoGroup(rv$metaData),
                        multiple = TRUE)
          })
          
          # print the experimental levels
          output$experimentallevels_rnaseq_norm <- renderText({
            req(input$groupselect_rnaseq_norm)
            
            if(length(input$groupselect_rnaseq_norm) > 1){
              experimentFactor <- make.names(apply(rv$metaData[,input$groupselect_rnaseq_norm], 1, paste, collapse = "_" ))
            } else{
              experimentFactor <- make.names(rv$metaData[,input$groupselect_rnaseq_norm])
            }
            levels <- unique(experimentFactor)
            if (length(levels) <= 6){
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(levels, collapse = ", "),
                             "</p>")
            } else {
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(head(levels), collapse = ", "),
                             "...</p>")
            }
            return(text)
            
          })
          
          # Automatically select transformation
          output$UI_transformation_rnaseq_norm <- renderUI({
            
            if (checkTransformation(gxMatrix = rv$gxData) == "Count matrix"){
              tagList(
                prettyRadioButtons(inputId = "transformation_rnaseq_norm", 
                                   label = NULL, 
                                   choices = c("Log2-transformation",
                                               "None"),
                                   inline = TRUE, 
                                   status = "danger",
                                   fill = TRUE,
                                   selected = "Log2-transformation"),
                h5(strong("NOTE: "),"Data seems not to be log-transformed. So, 
                 log2-transformation is needed.")
                
              )
              
            }else{
              tagList(
                prettyRadioButtons(inputId = "transformation_rnaseq_norm", 
                                   label = NULL, 
                                   choices = c("Log2-transformation",
                                               "None"),
                                   inline = TRUE, 
                                   status = "danger",
                                   fill = TRUE,
                                   selected = "None"),
                h5(strong("NOTE: "),"Data seems to be already log-transformed. So, 
                   no additional transformation is needed.")
              )
            }
          })
          
          # Automatically select filtering
          output$UI_filtering_rnaseq_norm <- renderUI({
            
            if (checkFiltering(gxMatrix = rv$gxData) == "No filtering"){
              tagList(
                materialSwitch(
                  inputId = "filtering_rnaseq_norm",
                  label = "Gene filtering?",
                  value = TRUE, 
                  status = "danger"),
                
                conditionalPanel(
                  condition = "input.filtering_rnaseq_norm==true",
                  numericInput(
                    inputId = "countfilter_rnaseq_norm",
                    label = "Minimum number of counts in smallest group size",
                    value = 10)
                  
                ),
                
                h5(strong("NOTE: "),"Data seems not to be filtered. So, 
                 (additional) filtering might be needed.")
                
              )
              
            }else{
              tagList(
                materialSwitch(
                  inputId = "filtering_rnaseq_norm",
                  label = "Gene filtering?",
                  value = FALSE, 
                  status = "danger"),
                
                conditionalPanel(
                  condition = "input.filtering_rnaseq_norm==true",
                  numericInput(
                    inputId = "countfilter_rnaseq_norm",
                    label = "Minimum number of counts in smallest group size",
                    value = 10)
                  
                ),
                
                h5(strong("NOTE: "),"Data seems to be already filtered. So, 
                   additional filtering might not be needed.")
              )
            }
          })
          
          # 3. pre-process the raw data
          observeEvent(input$start_preprocessing_rnaseq_norm,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Pre-processing data...",
                                          color="#0dc5c1")
            
            # Select outlier
            if (!isTRUE(input$outier_rnaseq_norm)){
              rv$outlier <- input$select_outliers_rnaseq_norm
            } else{
              rv$outlier <- NULL
            }
            
            # Remove outlier
            if (!is.null(rv$outlier)){
              gxMatrix_temp <- rv$gxData
              rv$gxData_fil <- gxMatrix_temp[,setdiff(colnames(gxMatrix_temp),rv$outlier)]
              rm(gxMatrix_temp)
            } else
              rv$gxData_fil <- rv$gxData
            
            # Filter metadata and expression data (samples in correct order)
            rv$metaData_fil <- rv$metaData[colnames(rv$gxData_fil),]
            
            # Experiment factor
            if(length(input$groupselect_rnaseq_norm) > 1){
              rv$experimentFactor <- factor(make.names(apply(rv$metaData_fil[,input$groupselect_rnaseq_norm], 1, paste, collapse = "_" )))
              rv$experimentName <- input$groupselect_rnaseq_norm
            } else{
              rv$experimentFactor <- factor(make.names(rv$metaData_fil[,input$groupselect_rnaseq_norm]))
              rv$experimentName <- input$groupselect_rnaseq_norm
            }
            
            # Normalization
            if (isTRUE(input$filtering_rnaseq_norm)){
              rv$normData <- RNASeqNormalization_processed(gxMatrix = rv$gxData_fil,
                                                           experimentFactor = rv$experimentFactor,
                                                           transMeth = input$transformation_rnaseq_norm,
                                                           normMeth = input$normMeth_rnaseq_norm,
                                                           perGroup_name = input$perGroup_rnaseq_norm,
                                                           filterThres = input$countfilter_rnaseq_norm)
            } else{
              rv$normData <- RNASeqNormalization_processed(gxMatrix = rv$gxData_fil,
                                                           experimentFactor = rv$experimentFactor,
                                                           transMeth = input$transformation_rnaseq_norm,
                                                           normMeth = input$normMeth_rnaseq_norm,
                                                           perGroup_name = input$perGroup_rnaseq_norm,
                                                           filterThres = 0)
            }
            
            # Collect chosen pre-processing settings into a dataframe
            rv$processingSettings <- data.frame(
              Option = c("Removed samples",
                         "Experimental group",
                         "Experimental levels",
                         "Normalization method",
                         "Transformation method",
                         "Filter threshold",
                         "Smallest group size"),
              Selected = c(paste(rv$outlier, collapse = "; "),
                           paste(input$groupselect_rnaseq_norm, collapse = "; "),
                           paste(unique(rv$experimentFactor), collapse = "; "),
                           paste(input$normMeth_rnaseq_norm, input$perGroup_rnaseq_norm, collapse = "; "),
                           input$transformation_rnaseq_norm,
                           ifelse(input$filtering_rnaseq_norm, input$countfilter_rnaseq_norm,0),
                           min(table(rv$experimentFactor))
              )
            )
            
            #======================================================================#
            # QC
            #======================================================================#
            
            #********************************************************************#
            # Expression values
            #********************************************************************#
            
            # Print expression table
            output$exprTable_rnaseq_norm <- DT::renderDataTable({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show message
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Info",
                text = "The data has been pre-processed. Please check the different 
              QC plots on this page to assess the pre-processing quality.",
                type = "info")
              
              # Show microarray statistics tab
              showTab("navbar", target = "panel_statistics_rnaseq_norm")
              
              output <- rv$normData
              return(output)
              
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$downloadNormalizedData_rnaseq_norm <- downloadHandler(
              filename = "normalizedData.csv",
              content = function(file){
                write.csv(rv$normData, file, quote = FALSE, row.names = TRUE)
              }
            )
            
            # Change color by click on button
            rv$colorOrder <- 1:length(levels(rv$experimentFactor))
            observeEvent(input$geneboxplot_changeOrder_rnaseq_norm,{
              all_orders <- permute(1:length(levels(rv$experimentFactor))) 
              sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
              
              if (sel < length(all_orders)){
                rv$colorOrder <- all_orders[[sel+1]]
              } else{
                rv$colorOrder <- all_orders[[1]]
              }
            })
            
            # Boxplot of single gene (based on selected row in the expression matrix)
            output$ExprBoxplot_rnaseq_norm <- renderPlot({
              req(input$exprTable_rnaseq_norm_rows_selected)
              
              if (length(levels(rv$experimentFactor)) > 4){
                legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
              } else{
                legendColors <- c(input$geneboxplot_col1_rnaseq_norm,
                                  input$geneboxplot_col2_rnaseq_norm,
                                  input$geneboxplot_col3_rnaseq_norm,
                                  input$geneboxplot_col4_rnaseq_norm,
                                  input$geneboxplot_col5_rnaseq_norm,
                                  input$geneboxplot_col6_rnaseq_norm)
              }
              # Make boxplot
              rv$temp1 <- geneBoxplot(experimentFactor = rv$experimentFactor, 
                                     normMatrix = rv$normData, 
                                     sel_row = input$exprTable_rnaseq_norm_rows_selected,
                                     legendColors = legendColors[rv$colorOrder],
                                     groupOrder = input$geneboxplot_order_rnaseq_norm,
                                     rnaseq = TRUE,
                                     jitter = input$jitter_geneboxplot_rnaseq_norm,
                                     seed = sample(1:1000,1))
              return(rv$temp1)
            })
            
            
            # Get number of experimental groups
            output$length_geneboxplot_rnaseq_norm <- reactive({
              length(levels(rv$experimentFactor))
            })
            outputOptions(output, "length_geneboxplot_rnaseq_norm", suspendWhenHidden = FALSE) 
            
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_geneboxplot_rnaseq_norm <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp1, 
                                filename = file,
                                width = input$width_geneboxplot_rnaseq_norm,
                                height = input$height_geneboxplot_rnaseq_norm,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_geneboxplot_rnaseq_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_geneboxplot_rnaseq_norm", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_geneboxplot_rnaseq_norm", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_geneboxplot_rnaseq_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Boxplots
            #********************************************************************#
            
            # Boxplots of all genes together
            output$boxplots_rnaseq_norm <- renderImage({
              req(session$clientData$output_boxplots_rnaseq_norm_width)
              req(session$clientData$output_boxplots_rnaseq_norm_height)
              getBoxplots(experimentFactor = rv$experimentFactor,
                          normData = rv$normData,
                          RNASeq = TRUE,
                          width = session$clientData$output_boxplots_rnaseq_norm_width,
                          height = session$clientData$output_boxplots_rnaseq_norm_height)
            }, deleteFile = TRUE)
            
            #***************************#
            # Modal to download figure
            #***************************#
            
            # Download plot
            output$realdownload_boxplots_rnaseq_norm <- downloadHandler(
              filename = function(){"QC_Boxplots.png"},
              content = function(file){
                png(file,
                    width=input$width_boxplots_rnaseq_norm,
                    height=input$height_boxplots_rnaseq_norm,
                    pointsize=24)
                getBoxplots_download(experimentFactor = rv$experimentFactor,
                                     normData = rv$normData,
                                     RNASeq = TRUE)
                dev.off()
              }
            )
            
            
            # Make modal
            observeEvent(input$download_boxplots_rnaseq_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_boxplots_rnaseq_norm", 
                                       "Height",
                                       min = 1200, max = 1600,
                                       value = 1500, step = 1,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_boxplots_rnaseq_norm", 
                                       "Width",
                                       min = 800, max = 1200,
                                       value = 1000, step = 1,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_boxplots_rnaseq_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Densityplots
            #********************************************************************#
            
            # Densityplots of all genes together
            output$densityplots_rnaseq_norm <- plotly::renderPlotly({
              getDensityplots(experimentFactor = rv$experimentFactor,
                              normMatrix = rv$normData,
                              RNASeq = TRUE)
              
            })
            
            #********************************************************************#
            # Heatmap
            #********************************************************************#
            
            # Heatmap of sample-sample correlations
            output$heatmap_rnaseq_norm  <- plotly::renderPlotly({
              
              # Make colors
              if(length(input$heatmapFactor_rnaseq_norm) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$heatmapFactor_rnaseq_norm], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$heatmapFactor_rnaseq_norm])
              }
              
              getHeatmap(experimentFactor = colorFactor, 
                         normMatrix = rv$normData,
                         clusterOption1 = input$clusteroption1_rnaseq_norm,
                         clusterOption2 = input$clusteroption2_rnaseq_norm,
                         theme = input$heatmaptheme_rnaseq_norm)
            })
            
            
            #********************************************************************#
            # PCA
            #********************************************************************#
            
            #Perform PCA
            rv$PCA_data <- prcomp(t(rv$normData[apply(rv$normData, 1, var) != 0,]),
                                  retx = TRUE, 
                                  center = TRUE,
                                  scale.= TRUE)
            
            
            # Make PCA plot
            output$PCA_rnaseq_norm <- plotly::renderPlotly({
              
              if(length(input$colorFactor_rnaseq_norm) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_rnaseq_norm], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$colorFactor_rnaseq_norm])
              }
              
              plot_PCA(PC_data = rv$PCA_data, 
                       colorFactor = colorFactor, 
                       xpc = as.numeric(stringr::str_remove(input$xpca_rnaseq_norm,"PC")), 
                       ypc = as.numeric(stringr::str_remove(input$ypca_rnaseq_norm,"PC")), 
                       zpc = ifelse(input$xyz_rnaseq_norm,as.numeric(stringr::str_remove(input$zpca_rnaseq_norm,"PC")),3), 
                       xyz = input$xyz_rnaseq_norm)
              
            })
            
            #********************************************************************#
            # Output 6: Overview of pre-processing settings
            #********************************************************************#
            # Print table with settings
            output$processingSettings_rnaseq_norm <- DT::renderDataTable({
              return(rv$processingSettings)
            },options = list(pageLength = 10),
            selection = "none")
            
            # Download button
            output$downloadProcessingSettings_rnaseq_norm <- downloadHandler(
              filename = "preprocessingSettings.csv",
              content = function(file){
                write.csv(rv$processingSettings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
            #********************************************************************#
            # UI
            #********************************************************************#
            
            # Render UI for main tab
            output$UI_QC_rnaseq_norm <- renderUI({
              tagList(
                tabsetPanel(
                  tabPanel("Expression values",
                           icon = icon("fas fa-mouse-pointer"),
                           h3(strong("Normalized expression values")),
                           h5("Here you can view the normalized log intensity (expression) values."),
                           hr(),
                           dataTableOutput(outputId = "exprTable_rnaseq_norm") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadNormalizedData_rnaseq_norm", 
                                          "Download"),
                           br(),
                           br(),
                           
                           # Dropdown Button to adjust the plot settings
                           shinyWidgets::dropdownButton(
                             tags$div(
                               style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                               
                               
                               # Change order of the boxplots:
                               tags$h4("Drag to change boxplot order"),
                               shinyjqui::orderInput(inputId = 'geneboxplot_order_rnaseq_norm', 
                                                     label = NULL, 
                                                     items = levels(rv$experimentFactor),
                                                     item_class = 'default'),
                               br(),
                               
                               # Change colour of the boxplots by button. 
                               # This is used when there are more than 6 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_rnaseq_norm > 6",
                                 tags$h4("Click to change boxplot colours"),
                                 shinyWidgets::actionBttn("geneboxplot_changeOrder_rnaseq_norm",
                                                          label = "Change color",
                                                          style = "simple",
                                                          color = "primary",
                                                          icon = icon("sync"))
                               ),
                               
                               # Change colour of the boxplots by colour picker
                               # This is used when there are less than 7 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_rnaseq_norm < 7",
                                 tags$h4("Click to select boxplot colours"),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 0",
                                   colourpicker::colourInput("geneboxplot_col1_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[1])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 1",
                                   colourpicker::colourInput("geneboxplot_col2_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[2])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 2",
                                   colourpicker::colourInput("geneboxplot_col3_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[3])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 3",
                                   colourpicker::colourInput("geneboxplot_col4_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[4])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 4",
                                   colourpicker::colourInput("geneboxplot_col5_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[5])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_rnaseq_norm > 5",
                                   colourpicker::colourInput("geneboxplot_col6_rnaseq_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[6])
                                 )
                               ),
                               br(),
                               tags$h4("Drag to change jitter"),
                               sliderInput("jitter_geneboxplot_rnaseq_norm", 
                                           NULL,
                                           min = 0, max = 0.3,
                                           value = 0.1, step = 0.01),
                               br()
                             ),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"),
                             
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to personalize plot!")
                             
                           ), # EO dropdownButton
                           
                           # Boxplot of the selected gene's expression values
                           plotOutput("ExprBoxplot_rnaseq_norm")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1"),
                           
                           actionButton("download_geneboxplot_rnaseq_norm", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           br(),
                           br()
                  ),
                  tabPanel("Boxplots",
                           icon = icon("fas fa-file"),
                           br(),
                           actionButton("download_boxplots_rnaseq_norm", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           plotOutput(outputId = "boxplots_rnaseq_norm",
                                      width = "65vw", height = "80vw")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  tabPanel("Density plots",
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           h2(strong("Density plot of normalized counts"), align = "center"),
                           h4("Distributions should be comparable between samples", align = "center"),
                           plotly::plotlyOutput(outputId = "densityplots_rnaseq_norm",
                                                width = "65vw", height = "40vw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  tabPanel("Correlation plot", 
                           icon = icon("fas fa-mouse-pointer"),
                           
                           br(),
                           
                           fluidRow(
                             
                             #Distance method
                             column(4,
                                    pickerInput(
                                      inputId = "clusteroption1_rnaseq_norm",
                                      label = "Distance calculation 
                                              method",
                                      choices = c("Pearson","Spearman",
                                                  "Euclidean"),
                                      options = list(
                                        style = "btn-primary"))
                             ),
                             
                             #Ckustering method
                             column(4,
                                    pickerInput(
                                      inputId = "clusteroption2_rnaseq_norm",
                                      label = "Clustering method",
                                      choices = c("ward.D2","single",
                                                  "complete","average",
                                                  "mcquitty","median",
                                                  "centroid"),
                                      options = list(
                                        style = "btn-info"))
                             )),
                           
                           hr(),
                           
                           #Theme
                           dropdownButton(
                             tags$h3("Theme"),
                             
                             selectInput(inputId = "heatmapFactor_rnaseq_norm",
                                         label = "Side colors",
                                         choices = colnames(rv$metaData_fil),
                                         selected = rv$experimentName,
                                         multiple = TRUE),
                             
                             selectInput(inputId = 'heatmaptheme_rnaseq_norm',
                                         label = "Heatmap theme",
                                         choices = c("Default", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"), width = "300px",
                             
                             tooltip = tooltipOptions(
                               title = "Click to change colors!")
                             
                           ),
                           
                           plotlyOutput("heatmap_rnaseq_norm", 
                                        width = "1000px", 
                                        height="600px") %>% 
                             withSpinner(color="#0dc5c1", 
                                         proxy.height = "400px")
                           
                  ),
                  tabPanel("PCA",
                           icon = icon("fas fa-mouse-pointer"),
                           # Title + text
                           fluidRow(
                             h3(strong("Principal Component Analysis (PCA)")),
                             h5("Here you can view the PCA score plot."),
                             hr(),
                           ),
                           
                           # Set color + 3D/2D
                           fluidRow(
                             column(3,
                                    # Color by which factor?
                                    pickerInput(inputId = "colorFactor_rnaseq_norm",
                                                label = "Color by:",
                                                choices = colnames(rv$metaData_fil),
                                                selected = rv$experimentName,
                                                multiple = TRUE)
                             ),
                             column(3,
                                    br(),
                                    # 3D plot?
                                    materialSwitch(
                                      inputId = "xyz_rnaseq_norm",
                                      label = "3D",
                                      value = FALSE, 
                                      status = "danger")
                                    
                             )
                           ),
                           
                           # Set axes
                           fluidRow(
                             column(3,
                                    #X-axis
                                    selectInput(inputId = "xpca_rnaseq_norm", 
                                                label = "x-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC1")
                             ),
                             column(3,
                                    #Y-axis
                                    selectInput(inputId = "ypca_rnaseq_norm", 
                                                label = "y-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC2")
                             ),
                             column(3,
                                    #Z-axis
                                    conditionalPanel(
                                      condition = "input.xyz_rnaseq_norm==true",
                                      selectInput(inputId = "zpca_rnaseq_norm", 
                                                  label = "z-axis",
                                                  choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                              "PC6", "PC7", "PC8"),
                                                  selected = "PC3")
                                    )
                             )
                           ),
                           
                           # Print plot
                           fluidRow(
                             hr(),
                             plotlyOutput("PCA_rnaseq_norm")%>% 
                               withSpinner(color="#0dc5c1")
                           ) 
                           
                  ), # EO PCA tabPanel
                  
                  # TAB6: Settings table
                  tabPanel("Settings overview",
                           icon = icon("fas fa-file"),
                           h3(strong("Pre-processing settings")),
                           h5("Here you can see an overview of the chosen pre-processing settings."),
                           hr(),
                           DT::dataTableOutput(outputId = "processingSettings_rnaseq_norm") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadProcessingSettings_rnaseq_norm", 
                                          "Download"),
                           
                  ) # EO Settings tabPanel
                ) # EO tabsetPanel
              ) # EO tagList
            }) # EO renderUI
            
            # Allow user to go to next tab
            output$UI_next_preprocessing_rnaseq_norm <- renderUI({
              req(rv$normData)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                actionBttn(inputId = "next_preprocessing_rnaseq_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          
          
          #======================================================================#
          # Statistical analysis
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_preprocessing_rnaseq_norm,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_statistics_rnaseq_norm")
          })
          
          # SIDEPANEL:
          
          # Select continuous covariates
          output$UI_covGroups_num_rnaseq_norm <- renderUI({
            tagList(
              pickerInput(inputId = "covGroups_num_rnaseq_norm",
                          label = "Continuous covariates (e.g., age):",
                          choices = setdiff(colnames(rv$metaData_fil),
                                            rv$experimentName),
                          selected = NULL,
                          multiple = TRUE)
            )
          })
          
          # Select discrete covariates
          output$UI_covGroups_char_rnaseq_norm <- renderUI({
            tagList(
              pickerInput(inputId = "covGroups_char_rnaseq_norm",
                          label = "Discrete covariates (e.g., sex):",
                          choices = setdiff(colnames(rv$metaData_fil),
                                            rv$experimentName),
                          selected = NULL,
                          multiple = TRUE)
            )
          })
          
          # Select comparisons
          output$UI_comparisons_rnaseq_norm <- renderUI({
            
            tagList(
              multiInput(
                inputId = "comparisons_rnaseq_norm",
                label = "Comparisons:", 
                choices = makeComparisons(make.names(unique(rv$experimentFactor))),
                selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1]
              )
            )
          })
          
          
          # Select comparisons
          output$UI_biomart_dataset_rnaseq_norm <- renderUI({
            req(input$addAnnotation_rnaseq_norm)
            pickerInput(inputId = "biomart_dataset_rnaseq_norm",
                        label = "Organism:",
                        choices = c("Homo sapiens" = "hsapiens_gene_ensembl" ,
                                    "Bos taurus" = "btaurus_gene_ensembl",
                                    "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                    "Mus musculus" = "mmusculus_gene_ensembl",
                                    "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                        selected = "hsapiens_gene_ensembl",
                        multiple = FALSE)
          })
          
          output$UI_addAnnotations_rnaseq_norm <- renderUI({
            req(input$addAnnotation_rnaseq_norm)
            req(input$biomart_dataset_rnaseq_norm)
            
            tagList(
              
              pickerInput(inputId = "biomart_filter_rnaseq_norm",
                          label = "Gene ID",
                          choices = c("Ensembl Gene ID",
                                      "Entrez Gene ID",
                                      "Gene Symbol/Name"),
                          selected = "Entrez Gene ID",
                          multiple = FALSE),
              
              pickerInput(inputId = "biomart_attributes_rnaseq_norm",
                          label = "Output",
                          choices = c("Ensembl Gene ID",
                                      "Entrez Gene ID",
                                      "Gene Symbol/Name"),
                          selected = "Gene Symbol/Name",
                          multiple = TRUE)
            )
            
          })
          
          
          
          #======================================================================#
          # Output of statistical analysis
          #======================================================================#
          observeEvent(input$calculate_statistics_rnaseq_norm,{
            show_modal_spinner(text = "Statistical analysis...",
                               color="#0dc5c1")
            
            # Calculate statistics
            if (isTRUE(input$addAnnotation_rnaseq_norm)){
              
              counts <- 2^rv$normData - 1 
              rv$top_table_list <- getStatistics_RNASeq_processed(normMatrix = counts, 
                                                                  metaData = rv$metaData_fil, 
                                                                  expFactor = rv$experimentName, 
                                                                  covGroups_num = input$covGroups_num_rnaseq_norm, 
                                                                  covGroups_char = input$covGroups_char_rnaseq_norm, 
                                                                  comparisons = input$comparisons_rnaseq_norm,
                                                                  addAnnotation = input$addAnnotation_rnaseq_norm,
                                                                  biomart_dataset = input$biomart_dataset_rnaseq_norm,
                                                                  biomart_attributes = unique(c(input$biomart_filter_rnaseq_norm,
                                                                                                input$biomart_attributes_rnaseq_norm)),
                                                                  biomart_filters = input$biomart_filter_rnaseq_norm)
              
            } else{
              counts <- 2^rv$normData - 1 
              rv$top_table_list <- getStatistics_RNASeq_processed(normMatrix = counts, 
                                                                  metaData = rv$metaData_fil, 
                                                                  expFactor = rv$experimentName, 
                                                                  covGroups_num = input$covGroups_num_rnaseq_norm, 
                                                                  covGroups_char = input$covGroups_char_rnaseq_norm, 
                                                                  comparisons = input$comparisons_rnaseq_norm,
                                                                  addAnnotation = input$addAnnotation_rnaseq_norm,
                                                                  biomart_dataset = NULL,
                                                                  biomart_attributes = NULL,
                                                                  biomart_filters = NULL)
              
            }
            
            # Select comparison for output
            observe({
              
              
              if (!is.null(rv$top_table_list)){
                rv$top_table <- rv$top_table_list[[1]]
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show comparisons
                output$UI_comparisons_view_rnaseq_norm <- renderUI({
                  pickerInput(inputId = "comparisons_view_rnaseq_norm",
                              label = "Select comparison:",
                              choices = names(rv$top_table),
                              selected = names(rv$top_table)[1],
                              multiple = FALSE)
                })
                
                # Show message
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = rv$top_table_list[[2]],
                  type = "info")
                
                # Show microarray statistics tab
                showTab("navbar", target = "panel_ORA_rnaseq_norm")
              } else{
                shinybusy::remove_modal_spinner()
                # Show message
                sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = "Oops...something went wrong! Please try again!",
                  type = "error")
              }
              
            })
            
            
            # Make outputs:
            
            #********************************************************************#
            # top table
            #********************************************************************#
            
            observe({
              req(input$comparisons_view_rnaseq_norm)
              # print top table
              output$top_table_rnaseq_norm <- DT::renderDataTable({
                
                output <- rv$top_table[[input$comparisons_view_rnaseq_norm]]
                
                return(output)
              },options = list(pageLength = 6),
              selection = list(mode = "single", selected = 1), escape = FALSE)
              
              # Download button
              output$download_top_table_rnaseq_norm <- downloadHandler(
                filename = paste0("topTable_",input$comparisons_view_rnaseq_norm,".csv"),
                content = function(file){
                  write.csv(rv$top_table[[input$comparisons_view_rnaseq_norm]], file, quote = FALSE, row.names = FALSE)
                }
              )
            })
              
            # Change plotting data depending on whether all experimental groups will be plotted  
            observe({
              req(input$comparisons_view_rnaseq_norm)
              req(rv$top_table)
              
              if (!is.null(input$boxplotAll_rnaseq_norm)){
                if (input$boxplotAll_rnaseq_norm){
                  rv$newFactor <- rv$experimentFactor
                  rv$newData <- rv$normData
                }
                if (!input$boxplotAll_rnaseq_norm){
                  if(length(rv$experimentName) > 1){
                    t <- make.names(apply(rv$metaData_fil[,rv$experimentName], 1, paste, collapse = "_" ))
                  } else{
                    t <- make.names(rv$metaData_fil[,rv$experimentName])
                  }
                  rv$newData <- rv$normData[,t %in% (unlist(stringr::str_split(input$comparisons_view_rnaseq_raw, " - ")))]
                  rv$newFactor <- factor(as.character(rv$experimentFactor)[t %in% (unlist(stringr::str_split(input$comparisons_view_rnaseq_raw, " - ")))])
                }
              }
              
              # Change color by click on button
              rv$colorOrder <- 1:length(levels(rv$newFactor))
              observeEvent(input$statboxplot_changeOrder_rnaseq_norm,{
                all_orders <- permute(1:length(levels(rv$newFactor))) 
                sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
                
                if (sel < length(all_orders)){
                  rv$colorOrder <- all_orders[[sel+1]]
                } else{
                  rv$colorOrder <- all_orders[[1]]
                }
              })
              
              
              # Boxplot of single gene (based on selected row in the top table)
              output$ExprBoxplot_statistics_rnaseq_norm <- renderPlot({
                req(input$top_table_rnaseq_norm_rows_selected)
                req(rv$top_table)
                
                if (length(levels(rv$newFactor)) > 6){
                  legendColors <- colorsByFactor(rv$newFactor)$legendColors
                } else{
                  legendColors <- c(input$statboxplot_col1_rnaseq_norm,
                                    input$statboxplot_col2_rnaseq_norm,
                                    input$statboxplot_col3_rnaseq_norm,
                                    input$statboxplot_col4_rnaseq_norm,
                                    input$statboxplot_col5_rnaseq_norm,
                                    input$statboxplot_col6_rnaseq_norm)
                }
                
                gene <- rv$top_table[[input$comparisons_view_rnaseq_norm]]$GeneID[input$top_table_rnaseq_norm_rows_selected]
                sel_row <- which(as.character(rownames(rv$normData)) %in% as.character(gene))
                
                # Make boxplot
                rv$temp <- geneBoxplot(experimentFactor = rv$newFactor, 
                                       normMatrix = rv$newData, 
                                       sel_row = sel_row,
                                       legendColors = legendColors[rv$colorOrder],
                                       groupOrder = input$statboxplot_order_rnaseq_norm,
                                       jitter = input$jitter_statboxplot_rnaseq_norm,
                                       rnaseq=TRUE,
                                       seed = sample(1:1000,1))
                
                return(rv$temp)
                
              })
              
              # Get number of experimental groups
              output$length_statboxplot_rnaseq_norm <- reactive({
                length(levels(rv$experimentFactor))
              })
              outputOptions(output, "length_statboxplot_rnaseq_norm", suspendWhenHidden = FALSE) 
              
            })
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_statboxplot_rnaseq_norm <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp, 
                                filename = file,
                                width = input$width_statboxplot_rnaseq_norm,
                                height = input$height_statboxplot_rnaseq_norm,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_statboxplot_rnaseq_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_statboxplot_rnaseq_norm", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_statboxplot_rnaseq_norm", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_statboxplot_rnaseq_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            
            #********************************************************************#
            # Histograms
            #********************************************************************#
            observe({
              req(input$comparisons_view_rnaseq_norm)
              
              if (input$comparisons_view_rnaseq_norm %in% names(rv$top_table)){
                output$Phistogram_rnaseq_norm <- plotly::renderPlotly({
                  req(rv$top_table)
                  p <- makePHistogram(rv$top_table[[input$comparisons_view_rnaseq_norm]][,"p-value"])
                  return(p)
                })
                
                output$logFChistogram_rnaseq_norm <- renderPlotly({
                  req(rv$top_table)
                  p <- makelogFCHistogram(rv$top_table[[input$comparisons_view_rnaseq_norm]][,"log2FC"])
                  return(p)
                })
              }
            })
            
            #********************************************************************#
            # Volcano plot
            #********************************************************************#
            observeEvent(input$plot_volcano_rnaseq_norm, {
              req(rv$top_table)
              req(input$rawp_volcano_rnaseq_norm)
              req(input$p_thres_volcano_rnaseq_norm)
              req(input$logFC_thres_volcano_rnaseq_norm)
              req(input$comparisons_view_rnaseq_norm)
              
              if (input$comparisons_view_rnaseq_norm %in% names(rv$top_table)){
                p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_rnaseq_norm]], 
                                 p = input$rawp_volcano_rnaseq_norm, 
                                 p_threshold = input$p_thres_volcano_rnaseq_norm, 
                                 logFC_threshold = input$logFC_thres_volcano_rnaseq_norm)
                
                output$volcano_rnaseq_norm <- plotly::renderPlotly(p)
              }
            }, ignoreNULL = FALSE)
            
            
            #=========================================#
            #  # UI: Output in different tabs
            #=========================================#
            
            observe({
              req(rv$top_table)
              
              output$UI_boxplotAll_rnaseq_norm <- renderUI({
                tagList(
                  # Change order of the boxplots:
                  tags$h4("Drag to change boxplot order"),
                  shinyjqui::orderInput(inputId = 'statboxplot_order_rnaseq_norm', 
                                        label = NULL, 
                                        items = levels(rv$newFactor),
                                        item_class = 'default'),
                  br(),
                  
                  # Change colour of the boxplots by button. 
                  # This is used when there are more than 6 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_rnaseq_raw > 6",
                    tags$h4("Click to change boxplot colours"),
                    shinyWidgets::actionBttn("statboxplot_changeOrder_rnaseq_norm",
                                             label = "Change color",
                                             style = "simple",
                                             color = "primary",
                                             icon = icon("sync"))
                  ),
                  
                  # Change colour of the boxplots by colour picker
                  # This is used when there are less than 7 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_rnaseq_norm < 7",
                    tags$h4("Click to select boxplot colours"),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_norm > 0",
                      colourpicker::colourInput("statboxplot_col1_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[1])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_norm > 1",
                      colourpicker::colourInput("statboxplot_col2_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[2])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_norm > 2",
                      colourpicker::colourInput("statboxplot_col3_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[3])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_norm > 3",
                      colourpicker::colourInput("statboxplot_col4_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[4])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_norm > 4",
                      colourpicker::colourInput("statboxplot_col5_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[5])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_rnaseq_rorm > 5",
                      colourpicker::colourInput("statboxplot_col6_rnaseq_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[6])
                    )
                  ),
                  br(),
                  tags$h4("Drag to change jitter"),
                  sliderInput("jitter_statboxplot_rnaseq_norm", 
                              NULL,
                              min = 0, max = 0.3,
                              value = 0.1, step = 0.01),
                  br()
                )
              })
            })
            
            observe({
              if (is.null(rv$top_table)){
                output$UI_output_statistics_rnaseq_norm <- renderUI(NULL)
              } else{
                output$UI_output_statistics_rnaseq_norm <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #********************************************************************#
                      # top table tab
                      #********************************************************************#
                      
                      tabPanel("Top table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("Top Table")),
                               h5("The Top Table includes the output of the selected statistical analysis."),
                               hr(),
                               dataTableOutput(outputId = "top_table_rnaseq_norm") %>% 
                                 withSpinner(color="#0dc5c1"),
                               downloadButton("download_top_table_rnaseq_norm", 
                                              "Download table"),
                               br(),
                               br(),
                               
                               # Dropdown Button to adjust the plot settings
                               shinyWidgets::dropdownButton(
                                 tags$div(
                                   style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                                   
                                   # Plot all experimental groups?
                                   conditionalPanel(
                                     condition = "output.length_statboxplot_rnaseq_norm > 2",
                                     tags$h4("Plot all experimental groups?"),
                                     shinyWidgets::materialSwitch(inputId = "boxplotAll_rnaseq_norm",
                                                                  label = NULL, 
                                                                  value = TRUE,
                                                                  status = "primary"),
                                     
                                     br()
                                   ),
                                   uiOutput("UI_boxplotAll_rnaseq_norm")
                                 ),
                                 circle = TRUE, status = "info",
                                 icon = icon("fas fa-cog"),
                                 
                                 tooltip = shinyWidgets::tooltipOptions(
                                   title = "Click to personalize plot!")
                                 
                               ), # EO dropdownButton
                               plotOutput("ExprBoxplot_statistics_rnaseq_norm")%>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1"),
                               actionButton("download_statboxplot_rnaseq_norm", 
                                            "Download figure",
                                            icon = shiny::icon("download")),
                               br(),
                               br()
                      ),
                      #********************************************************************#
                      # histogram tab
                      #********************************************************************#
                      
                      tabPanel("Histograms",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               plotlyOutput("Phistogram_rnaseq_norm")%>% 
                                 withSpinner(color="#0dc5c1"),
                               hr(),
                               plotlyOutput("logFChistogram_rnaseq_norm")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      ),
                      
                      #********************************************************************#
                      # volcano tab
                      #********************************************************************#
                      
                      tabPanel("Volcano plot",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               fluidRow(
                                 column(3,
                                        # Raw or adjusted P value?
                                        prettyRadioButtons(
                                          inputId = "rawp_volcano_rnaseq_norm",
                                          label = "P value", 
                                          choices = 
                                            c("Raw P value" = "raw", 
                                              "Adjusted P value" = "adj"))
                                 ),
                                 column(3,
                                        #P value threshold
                                        numericInput(
                                          inputId = "p_thres_volcano_rnaseq_norm",
                                          label = "P threshold",
                                          value = 0.05)
                                 ),
                                 column(3,
                                        #logFC threshold
                                        numericInput(
                                          inputId = "logFC_thres_volcano_rnaseq_norm",
                                          label = "logFC threshold",
                                          value = 1)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_volcano_rnaseq_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               # Volcano plot
                               plotlyOutput("volcano_rnaseq_norm")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      )
                      
                    ) # tabsetpanel
                  ) # taglist
                }) # renderUI
              }
            }) # Observe
            
            # Allow user to go to next tab
            output$UI_next_statistics_rnaseq_norm <- renderUI({
              req(rv$top_table)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                actionBttn(inputId = "next_statistics_rnaseq_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          
          #======================================================================#
          # ORA
          #======================================================================#
          # Go to data upload tab
          observeEvent(input$next_statistics_rnaseq_norm,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_ORA_rnaseq_norm")
          })
          
          # Comparisons for which ORA should be performed
          observe({
            req(rv$top_table)
            output$UI_comparisons_view_ORA_rnaseq_norm <- renderUI({
              pickerInput(inputId = "comparisons_view_ORA_rnaseq_norm",
                          label = NULL,
                          choices = names(rv$top_table),
                          selected = names(rv$top_table)[1],
                          multiple = FALSE)
            })
          })
          
          # Which columns contain gene ids
          observe({
            req(input$comparisons_view_ORA_rnaseq_norm)
            col_choice <- 1
            if (length(colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]])) > 6){
              col_choice <- c(1,7:ncol(rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]]))
            }
            output$UI_geneID_ORA_rnaseq_norm <- renderUI({
              tagList(
                selectInput(inputId = "organism_ORA_rnaseq_norm",
                            label = "Organism",
                            choices = c("Bos taurus",
                                        "Caenorhabditis elegans",
                                        "Homo sapiens",
                                        "Mus musculus", 
                                        "Rattus norvegicus"),
                            selected = "Homo sapiens"),
                
                pickerInput(inputId = "geneID_ORA_rnaseq_norm",
                            label = "Which column of the top table contains the gene IDs?",
                            choices = colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]])[col_choice],
                            selected = colnames(rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]])[1],
                            multiple = FALSE),
                
                pickerInput(inputId = "selID_ORA_rnaseq_norm",
                            label = "Which gene ID to use?",
                            choices = c("Ensembl Gene ID" = "ENSEMBL", 
                                        "Entrez Gene ID" = "ENTREZID", 
                                        "Gene Symbol/Name" = "SYMBOL"),
                            selected = "ENTREZID",
                            multiple = FALSE),
              )
            })
          })
          
          # Show modal
          observeEvent(input$calculate_ORA_rnaseq_norm,{
            show_modal_spinner(text = "Overrepresentation analysis...",
                               color="#0dc5c1")
            
            # Perform ORA
            if (input$topNorThres_rnaseq_norm == "Threshold"){
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_rnaseq_norm]],
                                 geneset = input$geneset_ORA_rnaseq_norm,
                                 geneID_col = input$geneID_ORA_rnaseq_norm,
                                 geneID_type = input$selID_ORA_rnaseq_norm,
                                 organism = input$organism_ORA_rnaseq_norm,
                                 updown = input$updown_ORA_rnaseq_norm,
                                 topN = FALSE,
                                 N = NULL,
                                 rawadj = input$rawp_ORA_rnaseq_norm,
                                 p_thres = input$p_thres_ORA_rnaseq_norm,
                                 logFC_thres = input$logFC_thres_ORA_rnaseq_norm)
            } else{
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]],
                                 geneset = input$geneset_ORA_rnaseq_norm,
                                 geneID_col = input$geneID_ORA_rnaseq_norm,
                                 geneID_type = input$selID_ORA_rnaseq_norm,
                                 organism = input$organism_ORA_rnaseq_norm,
                                 updown = input$updown_ORA_rnaseq_norm,
                                 topN = TRUE,
                                 N = input$topN_rnaseq_norm,
                                 rawadj = NULL,
                                 p_thres = NULL,
                                 logFC_thres = NULL)
            }
            
            observe({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show message
              if (is.null(rv$ORA_data)){
                sendSweetAlert(
                  session = session,
                  title = "Error!",
                  text = "No significant genes!",
                  type = "error")
                
              }else{
                
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = "Gene overrepresentation analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
                  type = "info")
                
                # print GO table
                output$ORA_table_rnaseq_norm <- DT::renderDataTable({
                  req(input$geneset_ORA_rnaseq_norm)
                  output <- rv$ORA_data@result
                  
                  if (input$geneset_ORA_rnaseq_norm == "WikiPathways"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.wikipathways.org/pathways/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_rnaseq_norm == "KEGG"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.genome.jp/pathway/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_rnaseq_norm %in% c("GO-BP", "GO-MF", "GO-CC")){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.ebi.ac.uk/QuickGO/term/",
                        output$ID
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  
                  return(output)
                },options = list(pageLength = 6),
                selection = list(mode = "single", selected = 1), escape = FALSE)
                
                # Download button
                output$download_ORA_table_rnaseq_norm <- downloadHandler(
                  filename = paste0("ORATable_",input$comparisons_view_ORA_rnaseq_norm,"_",input$geneset_ORA_rnaseq_norm,".csv"),
                  content = function(file){
                    write.csv(rv$ORA_data@result, file, quote = FALSE, row.names = FALSE)
                  }
                )
                
                # Print statistics of genes in selected Term
                output$ORAgene_table_rnaseq_norm <- DT::renderDataTable({
                  req(input$ORA_table_rnaseq_norm_rows_selected)
                  
                  # Make ORA gene table
                  output <- make_ORAgene_table(ORA_data = rv$ORA_data,
                                               top_table = rv$top_table[[input$comparisons_view_ORA_rnaseq_norm]],
                                               geneID_col = input$geneID_ORA_rnaseq_norm,
                                               sel_row_ORA = input$ORA_table_rnaseq_norm_rows_selected)
                  
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  output$meanExpr <- round(output$meanExpr,3)
                  output$log2FC <- round(output$log2FC,3)
                  output$`log2FC SE` <- round(output$`log2FC SE`,3)
                  
                  return(output)
                }, options = list(pageLength = 6), escape = FALSE)
                
                # Text for gene table
                output$text_ORAgene_table_rnaseq_norm <- renderText({
                  text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_rnaseq_norm_rows_selected,"ID"],
                                 "</b></h3>")
                  return(text)
                })
                
                # ORA barchart
                observeEvent(input$plot_ORAplot_rnaseq_norm, {
                  req(input$nSets_ORAplot_rnaseq_norm)
                  p <- makeORAplot(rv$ORA_data,
                                   nSets = input$nSets_ORAplot_rnaseq_norm)
                  
                  output$ORAplot_rnaseq_norm <- plotly::renderPlotly(p)
                  
                }, ignoreNULL = FALSE)
                
                # ORA network diagram
                observeEvent(input$plot_ORAnetwork_rnaseq_norm, {
                  req(input$layout_network_rnaseq_norm)
                  req(input$nSets_network_rnaseq_norm)
                  p <- makeORAnetwork(ORA_data = rv$ORA_data,
                                      layout = input$layout_network_rnaseq_norm,
                                      nSets = input$nSets_network_rnaseq_norm)
                  
                  output$ORAnetwork_rnaseq_norm <- renderPlot(p,
                                                              height = 500, 
                                                              width = 800)
                  
                }, ignoreNULL = FALSE)
                
              }
            }) # Observe
            
            
            
            # Render output table
            observe({
              if (!is.null(rv$ORA_data)){
                output$UI_output_ORA_rnaseq_norm <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #********************************************************************#
                      # top table tab
                      #********************************************************************#
                      
                      tabPanel("ORA table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("Statistics table")),
                               h5("The ORA statistics table encompasses the output of the gene 
                              overrepresentation analysis."),
                               hr(),
                               dataTableOutput(outputId = "ORA_table_rnaseq_norm") %>% 
                                 withSpinner(color="#0dc5c1"),
                               downloadButton("download_ORA_table_rnaseq_norm", 
                                              "Download"),
                               br(),
                               htmlOutput("text_ORAgene_table_rnaseq_norm"),
                               h5(paste0("The gene table encompasses the statistics of all genes 
                              from the selected geneset.")),
                               hr(),
                               dataTableOutput(outputId = "ORAgene_table_rnaseq_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                      tabPanel("Bar chart",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("ORA bar chart")),
                               h5("The ORA bar chart visualizes the results from the overrepresentation analysis."),
                               hr(),
                               fluidRow(
                                 column(3,
                                        # Number of genesets
                                        numericInput(
                                          inputId = "nSets_ORAplot_rnaseq_norm",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_ORAplot_rnaseq_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               plotlyOutput("ORAplot_rnaseq_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                      tabPanel("Network diagram",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("ORA network diagram")),
                               h5("The ORA network diagram visualize the similarity between the most significant genesets."),
                               hr(),
                               fluidRow(
                                 column(3,
                                        # Number of genesets
                                        numericInput(
                                          inputId = "nSets_network_rnaseq_norm",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 ),
                                 column(3,
                                        # Network layout
                                        pickerInput(inputId = "layout_network_rnaseq_norm",
                                                    label = "Network layout",
                                                    choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                                'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                    selected = 'graphopt',
                                                    multiple = FALSE)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_ORAnetwork_rnaseq_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               plotOutput("ORAnetwork_rnaseq_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                    ) # tabsetpanel
                  ) # taglist
                })
              }# renderUI
            }) # Observe
            
          }) #observe event
          
        } # raw or norm
        
        
      } # EO Microarray or RNA-seq
      
      
      
      ############################################################################
      
      # Microarray Analysis
      
      ############################################################################
      
      if (microarray_or_rnaseq() == "Microarray"){
        
        
        #************************************************************************#
        # Raw microarray data
        #************************************************************************#
        if (raw_or_norm() == "Raw data"){
          
          #======================================================================#
          # TAB 1: Data Upload
          #======================================================================#
          
          #----------------------------------------------------------------------#
          # Go to data upload tab
          #----------------------------------------------------------------------#
          # Observe "startAnalysis" input
          observeEvent(input$startAnalysis,{
            
            # Show microarray (raw) upload tab
            showTab("navbar", target = "panel_upload_microarray_raw")
            
            # Remove the other RNA-seq (raw) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_raw")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_raw")
            hideTab("navbar", target = "panel_statistics_rnaseq_raw" )
            hideTab("navbar", target = "panel_ORA_rnaseq_raw")
            
            # Remove the other RNA-seq (norm) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_norm")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_norm")
            hideTab("navbar", target = "panel_statistics_rnaseq_norm")
            hideTab("navbar", target = "panel_ORA_rnaseq_norm")
            
            # Remove the other microarray (raw) tabs
            hideTab("navbar", target = "panel_preprocessing_microarray_raw")
            hideTab("navbar", target = "panel_statistics_microarray_raw" )
            hideTab("navbar", target = "panel_ORA_microarray_raw")
            
            # Remove the other microarray (norm) tabs
            hideTab("navbar", target = "panel_upload_microarray_norm")
            hideTab("navbar", target = "panel_preprocessing_microarray_norm")
            hideTab("navbar", target = "panel_statistics_microarray_norm")
            hideTab("navbar", target = "panel_ORA_microarray_norm")
            
            # Go to microarray (raw) upload tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_upload_microarray_raw")
            
            # Example meta data file
            output$downloadmeta_example_microarray_raw <- downloadHandler(
              filename = "MetaData_example.csv",
              content = function(file){
                write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
              }
            )
          })
          
          #----------------------------------------------------------------------#
          # Upload files
          #----------------------------------------------------------------------#
          
          # Observe "upload" input
          observeEvent(input$upload_microarray_raw,{
            
            # Show loading modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Select (raw) CEL files
            rv$zippath <- input$uploadCEL_microarray_raw$datapath
            rv$celfiles <- getCELs(zippath = input$uploadCEL_microarray_raw)
            
            # Get metadata
            #req(rv$celfiles) # CEL files are required for the metadata table
            if(length(rv$celfiles)>0){
              
              # Metadata in tsv or csv format
              if (input$MetaFileType_microarray_raw==".tsv/.csv file"){
                rv$metaData <- getMetaData(path = input$uploadMeta_microarray_raw_tsv$datapath,
                                           celfiles = rv$celfiles,
                                           filetype = input$MetaFileType_microarray_raw)
              }
              
              # Metadata in Series Matrix File format
              if (input$MetaFileType_microarray_raw=="Series Matrix File"){
                rv$metaData <- getMetaData(path = input$uploadMeta_microarray_raw_smf$datapath,
                                           celfiles = rv$celfiles,
                                           filetype = input$MetaFileType_microarray_raw)
              }
            } else {
              
              # If there are no CEL files: give error message
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
            }
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # Read data
              rv$celfiles_fil <- rv$celfiles[stringr::str_remove(basename(rv$celfiles),"\\.CEL.*") %in% rownames(rv$metaData)]
              rv$gxData <- readCELs(rv$celfiles_fil, 
                                    zippath = rv$zippath, 
                                    rm = FALSE)
              
              # Check whether all expression samples have meta data available
              if (nrow(rv$metaData) < length(rv$celfiles)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              
              
              #------------------------------------------------------------------#
              # Outputs of uploaded files
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_microarray_raw <- renderDataTable(NULL)
              output$metaTable_microarray_raw <- renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_microarray_raw <- renderDataTable({
                req(rv$gxData)
                output <- head(exprs(rv$gxData),6)
                colnames(output) <- stringr::str_remove(colnames(output),"\\.CEL.*")
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_microarray_raw <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              # Allow the user to adjust meta data
              observeEvent(input$metaTable_microarray_raw_cell_edit, {
                row  <- input$metaTable_microarray_raw_cell_edit$row
                clmn <- input$metaTable_microarray_raw_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_microarray_raw_cell_edit$value
              })
              
              # Render UI for main tab
              output$UI_upload_microarray_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_microarray_raw <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_microarray_raw")
                
                # Show message
                if (nrow(rv$metaData) >= length(rv$celfiles)){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_microarray_raw",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              # Remove loading modal
              shinybusy::remove_modal_spinner()
            }
            
          }) # EO observeEvent
          
          #----------------------------------------------------------------------#
          # Run Example
          #----------------------------------------------------------------------#
          
          # Observe "upload" input
          observeEvent(input$example_microarray_raw,{
            
            # Show loading modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Select (raw) CEL files
            rv$zippath <- "Data/Microarray/GSE6955_RAW.zip"
            rv$celfiles <- getCELs(zippath = rv$zippath,
                                   shiny_upload = FALSE)
            
            # Get metadata
            req(rv$celfiles) # CEL files are required for the metadata table
            if(length(rv$celfiles)>0){
              
              # Metadata in tsv or csv format
              rv$metaData <- getMetaData(path = "Data/Microarray/metaData_GSE6955.csv",
                                         celfiles = rv$celfiles,
                                         filetype = ".tsv/.csv file")
            } else {
              # Send error message if there is no expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No expression data",
                type = "error")
            }
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # Read data
              rv$celfiles_fil <- rv$celfiles[stringr::str_remove(basename(rv$celfiles),"\\.CEL.*") %in% rownames(rv$metaData)]
              rv$gxData <- readCELs(rv$celfiles_fil, 
                                    zippath = rv$zippath,
                                    rm = FALSE)
              
              # Check whether all expression samples have meta data available
              if (nrow(rv$metaData) < length(rv$celfiles)){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              
              
              #------------------------------------------------------------------#
              # Outputs of example
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_microarray_raw <- DT::renderDataTable(NULL)
              output$metaTable_microarray_raw <- DT::renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_microarray_raw <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(exprs(rv$gxData),6)
                colnames(output) <- stringr::str_remove(colnames(output),"\\.CEL.*")
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_microarray_raw <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              # Allow the user to adjust the meta data
              observeEvent(input$metaTable_microarray_raw_cell_edit, {
                row  <- input$metaTable_microarray_raw_cell_edit$row
                clmn <- input$metaTable_microarray_raw_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_microarray_raw_cell_edit$value
              })
              
              # Render UI for main tab
              output$UI_upload_microarray_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    
                    # Tab with the expression matrix
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")),
                    
                    # Tab with the meta data
                    tabPanel("Meta data",                  
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_microarray_raw <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_microarray_raw")
                
                # Show success message (only if the meta data is available for 
                # each expression file)
                if (nrow(rv$metaData) >= length(rv$celfiles)){
                  
                  # Show information tab
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been succesfully uploaded. Please check 
                  the expression matrix and the meta data tables on this page 
                  to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_microarray_raw",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else { # if nrow(rv$metaData) == 0
              
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              # Remove loading modal
              shinybusy::remove_modal_spinner()
            }
            
          }) # EO observeEvent
          
          
          #======================================================================#
          # TAB 2: Data Pre-processing
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_upload_microarray_raw,{
            
            # Go to RNA-seq tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_preprocessing_microarray_raw")
            
            # Guess the organism
            rv$Organism <- getOrganism(gxData = rv$gxData)
          })
          
          #------------------------------------------------------------------#
          # Dynamic UI
          #------------------------------------------------------------------#
          
          # 1. Select outliers
          output$UI_outlier_microarray_raw <- renderUI({
            
            # If outliers are checked, show possible sample to select as outliers
            if(!input$outlier_microarray_raw){
              samples <- rownames(rv$metaData)
              shinyWidgets::pickerInput(inputId = "select_outliers_microarray_raw",
                                        label = tags$span(
                                          "Select samples to be removed", 
                                          tags$span(
                                            icon(
                                              name = "question-circle",
                                            ) 
                                          ) |>
                                            prompter::add_prompt(message = "Select one or more samples to exclude from the analysis.", 
                                                                 position = "right",
                                                                 size = "large")
                                        ),
                                        choices = samples,
                                        multiple = TRUE)
            } else{
              NULL
            }
            
          })
          
          # 2. Select experimental groups
          output$UI_groupselect_microarray_raw <- renderUI({
            shinyWidgets::pickerInput(inputId = "groupselect_microarray_raw",
                                      label = NULL,
                                      choices = colnames(rv$metaData),
                                      selected = autoGroup(rv$metaData),
                                      multiple = TRUE)
          })
          
          # print the experimental levels
          output$experimentallevels <- renderText({
            req(input$groupselect_microarray_raw)
            
            if(length(input$groupselect_microarray_raw) > 1){
              # If more than one metadata column is selected: concatenate  the vectors
              experimentFactor <- make.names(apply(rv$metaData[,input$groupselect_microarray_raw], 1, paste, collapse = "_" ))
            } else{
              # If only one metadata column is selected: nothing to do!
              experimentFactor <- make.names(rv$metaData[,input$groupselect_microarray_raw])
            }
            
            # Show the levels of the selected experimental factor
            levels <- unique(experimentFactor)
            if (length(levels) <= 6){
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(levels, collapse = ", "),
                             "</p>")
            } else {
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(head(levels), collapse = ", "),
                             "...</p>")
            }
            return(text)
            
          })
          
          # 3. Normalization
          # -
          
          # 4. Select Annotation
          
          # Select Organism
          output$UI_species_microarray_raw <- renderUI({
            selectInput(inputId = "species_microarray_raw",
                        label = "Species",
                        choices = c("Anopheles gambiae",
                                    "Arabidopsis thaliana",
                                    "Bos taurus",
                                    "Caenorhabditis elegans",
                                    "Canis familiaris", 
                                    "Danio rerio",
                                    "Drosophila melanogaster",
                                    "Gallus gallus",
                                    "Homo sapiens",
                                    "Macaca mulatta",
                                    "Mus musculus", 
                                    "Oryza sativa",
                                    "Rattus norvegicus",
                                    "Saccharomyces cerevisiae",
                                    "Schizosaccharomyces pombe",
                                    "Sus scrofa"),
                        selected = rv$Organism)
          })
          
          #------------------------------------------------------------------#
          # Perform pre-processing
          #------------------------------------------------------------------#
          
          # 5. Pre-process the raw data
          observeEvent(input$start_preprocessing_microarray_raw,{
            
            # Get Probeset annotation
            rv$annotations <- input$annotations_microarray_raw 
            if (rv$annotations == "Custom annotations"){
              rv$ProbeAnnotation <- input$CDFtype_microarray_raw 
              rv$Organism <- input$species_microarray_raw 
            } else{
              rv$ProbeAnnotation <- NULL
            }
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Pre-processing data...",
                                          color="#0dc5c1")
            
            # Select outlier
            if (!isTRUE(input$outier_microarray_raw)){
              rv$outlier <- input$select_outliers_microarray_raw
            } else{
              rv$outlier <- NULL
            }
            
            # Remove outlier
            rv$celfiles_sel <- rv$celfiles_fil
            if (!is.null(rv$outlier)){
              for (i in 1:length(rv$outlier)){
                rv$celfiles_sel <- rv$celfiles_sel[!stringr::str_detect(rv$celfiles_sel, rv$outlier[i])]
              }
            }
            
            # Filter metadata and expression data (and put samples in correct order)
            rv$gxData_fil <- readCELs(celfiles= rv$celfiles_sel, 
                                      zippath = rv$zippath, 
                                      rm = TRUE)
            rv$metaData_fil <- rv$metaData[stringr::str_remove(colnames(Biobase::exprs(rv$gxData_fil)),"\\.CEL.*"),]
            
            # Experiment factor
            if(length(input$groupselect_microarray_raw) > 1){
              rv$experimentFactor <- factor(make.names(apply(rv$metaData_fil[,input$groupselect_microarray_raw], 1, paste, collapse = "_" )))
              rv$experimentName <- input$groupselect_microarray_raw
            } else{
              rv$experimentFactor <- factor(make.names(rv$metaData_fil[,input$groupselect_microarray_raw]))
              rv$experimentName <- input$groupselect_microarray_raw
            }
            
            # If the data object is of a "geneFeatureSet" class, only RMA
            # normalization will be performed --> warn user about this!
            if (class(rv$gxData_fil) == "GeneFeatureSet"){
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Warning",
                text = "The data is from the Gene ST 2.x series of arrays. 
              Currently only RMA normalization without custom probeset annotation 
              is supported. So, the pre-processing will be performed accordingly.",
                type = "info")
              
              rv$annotations <- "No annotations"
              rv$ProbeAnnotation <- NULL 
            }
            
            # Collect chosen pre-processing settings into a dataframe
            rv$processingSettings <- data.frame(
              Option = c("Removed samples",
                         "Experimental group",
                         "Experimental levels",
                         "Normalization method",
                         "Annotation"),
              Selected = c(paste(rv$outlier, collapse = "; "),
                           paste(input$groupselect_microarray_raw, collapse = "; "),
                           paste(unique(rv$experimentFactor), collapse = "; "),
                           paste(input$normMeth_microarray_raw,"; ",input$perGroup_microarray_raw),
                           paste0(rv$annotations,
                                  ifelse(rv$annotations == "Custom annotations", paste0("; ",input$CDFtype_microarray_raw),""), 
                                  ifelse(rv$annotations == "Custom annotations", paste0("; ",input$species_microarray_raw), ""))
              )
            )
            
            # Perform normalization
            rv$normData <- microarrayNormalization(gxData = rv$gxData_fil,
                                                   experimentFactor = rv$experimentFactor,
                                                   normMeth = input$normMeth_microarray_raw,
                                                   CDFtype = input$CDFtype_microarray_raw,
                                                   species = input$species_microarray_raw,
                                                   annotations = rv$annotations,
                                                   perGroup_name = input$perGroup_microarray_raw,
                                                   annot_file_datapath = input$annot_file_microarray_raw$datapath)
            
            # Get expression values
            normMatrix <- Biobase::exprs(rv$normData)
            
            # If custom annotations were selected we need to remove the _at/_*
            # from the probe ids.
            if (rv$annotations == "Custom annotations"){
              id_names <- as.character(stringr::str_remove(rownames(normMatrix),"_.*"))
              normMatrix <- normMatrix[!duplicated(id_names),]
              rownames(normMatrix) <- id_names[!duplicated(id_names)]
            }
            rv$normMatrix <- normMatrix
            rm(normMatrix)
            
            
            #------------------------------------------------------------------#
            #Generate output: QC tables/plots
            #------------------------------------------------------------------#
            
            #******************************************************************#
            # Output 1: Table with normalized expression values
            #******************************************************************#
            
            # Print expression table
            output$exprTable_microarray_raw <- DT::renderDataTable({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show message
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Info",
                text = "The data has been pre-processed. Please check the different 
              QC plots on this page to assess the pre-processing quality.",
                type = "info")
              
              # Show microarray statistics tab
              showTab("navbar", target = "panel_statistics_microarray_raw")
              
              output <- rv$normMatrix
              
              # If custom probe annotations are selected, we can provide a link
              # to the NCBI or ENSEMBL websites
              if (!is.null(rv$ProbeAnnotation)){
                
                # ENTREZG: link to NCBI website
                if(rv$ProbeAnnotation == "ENTREZG"){
                  rownames(output) <- paste0(
                    '<a ',
                    'href=',
                    paste(
                      "https://www.ncbi.nlm.nih.gov/gene/",
                      rownames(output),
                      sep = ''
                    ),
                    ' target="_blank"',
                    '>',
                    rownames(output),
                    '</a>'
                  )
                }
                
                # ENSG: link to ENSEMBL website
                if(rv$ProbeAnnotation == "ENSG"){
                  rownames(output) <- paste0(
                    '<a ',
                    'href=',
                    paste(
                      "http://www.ensembl.org/id/",
                      rownames(output),
                      sep = ''
                    ),
                    ' target="_blank"',
                    '>',
                    rownames(output),
                    '</a>'
                  )
                }
              }
              return(output)
              
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$downloadNormalizedData_microarray_raw <- downloadHandler(
              filename = "normalizedData.csv",
              content = function(file){
                write.csv(rv$normMatrix, file, quote = FALSE, row.names = TRUE)
              }
            )
            
            # Change color by click on button
            rv$colorOrder <- 1:length(levels(rv$experimentFactor))
            observeEvent(input$geneboxplot_changeOrder_microarray_raw,{
              all_orders <- permute(1:length(levels(rv$experimentFactor))) 
              sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
              
              if (sel < length(all_orders)){
                rv$colorOrder <- all_orders[[sel+1]]
              } else{
                rv$colorOrder <- all_orders[[1]]
              }
            })
              
            
            # Boxplot of single gene (based on selected row in the expression matrix)
            output$ExprBoxplot_microarray_raw <- renderPlot({
              req(input$exprTable_microarray_raw_rows_selected)
              
              if (length(levels(rv$experimentFactor)) > 4){
                legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
              } else{
                legendColors <- c(input$geneboxplot_col1_microarray_raw,
                                  input$geneboxplot_col2_microarray_raw,
                                  input$geneboxplot_col3_microarray_raw,
                                  input$geneboxplot_col4_microarray_raw,
                                  input$geneboxplot_col5_microarray_raw,
                                  input$geneboxplot_col6_microarray_raw)
              }
              # Make boxplot
              rv$temp1 <- geneBoxplot(experimentFactor = rv$experimentFactor, 
                                     normMatrix = rv$normMatrix, 
                                     sel_row = input$exprTable_microarray_raw_rows_selected,
                                     legendColors = legendColors[rv$colorOrder],
                                     groupOrder = input$geneboxplot_order_microarray_raw,
                                     rnaseq = FALSE,
                                     jitter = input$jitter_geneboxplot_microarray_raw,
                                     seed = sample(1:1000,1))
              return(rv$temp1)
            })
            
            
            # Get number of experimental groups
            output$length_geneboxplot_microarray_raw <- reactive({
              length(levels(rv$experimentFactor))
            })
            outputOptions(output, "length_geneboxplot_microarray_raw", suspendWhenHidden = FALSE) 
            
            
            #***************************#
            # Modal to download boxplot
            #***************************#

            # Download plot
            output$realdownload_geneboxplot_microarray_raw <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp1, 
                                filename = file,
                                width = input$width_geneboxplot_microarray_raw,
                                height = input$height_geneboxplot_microarray_raw,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_geneboxplot_microarray_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_geneboxplot_microarray_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_geneboxplot_microarray_raw", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_geneboxplot_microarray_raw', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Output 2: Boxplots
            #********************************************************************#
            
            # Boxplots of all genes together
            output$boxplots_microarray_raw <- renderImage({
              req(session$clientData$output_boxplots_microarray_raw_width)
              req(session$clientData$output_boxplots_microarray_raw_height)
              getBoxplots(experimentFactor = rv$experimentFactor,
                          normData = rv$normData,
                          RNASeq = FALSE,
                          width = session$clientData$output_boxplots_microarray_raw_width,
                          height = session$clientData$output_boxplots_microarray_raw_height)
            }, deleteFile = TRUE)
            
            #***************************#
            # Modal to download figure
            #***************************#
            
            # Download plot
            output$realdownload_boxplots_microarray_raw <- downloadHandler(
              filename = function(){"QC_Boxplots.png"},
              content = function(file){
                png(file,
                    width=input$width_boxplots_microarray_raw,
                    height=input$height_boxplots_microarray_raw,
                    pointsize=24)
                getBoxplots_download(experimentFactor = rv$experimentFactor,
                                     normData = rv$normData,
                                     RNASeq = FALSE)
                dev.off()
              }
            )
            
            
            # Make modal
            observeEvent(input$download_boxplots_microarray_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_boxplots_microarray_raw", 
                                       "Height",
                                       min = 1200, max = 1600,
                                       value = 1500, step = 1,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_boxplots_microarray_raw", 
                                       "Width",
                                       min = 800, max = 1200,
                                       value = 1000, step = 1,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_boxplots_microarray_raw', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            
            #********************************************************************#
            # Output 3: Densityplots
            #********************************************************************#
            
            # Densityplots of all genes together
            output$densityplots_microarray_raw <- plotly::renderPlotly({
              getDensityplots(experimentFactor = rv$experimentFactor,
                              normMatrix = rv$normMatrix)
              
            })
            
            #********************************************************************#
            # Output 4: Sample-sample correlation heatmap
            #********************************************************************#
            
            # Heatmap of sample-sample correlations
            output$heatmap_microarray_raw  <- plotly::renderPlotly({
              
              # Make colors
              if(length(input$heatmapFactor_microarray_raw) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$heatmapFactor_microarray_raw], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$heatmapFactor_microarray_raw])
              }
              
              getHeatmap(experimentFactor = colorFactor, 
                         normMatrix = rv$normMatrix,
                         clusterOption1 = input$clusteroption1_microarray_raw,
                         clusterOption2 = input$clusteroption2_microarray_raw,
                         theme = input$heatmaptheme_microarray_raw)
            })
            
            #********************************************************************#
            # Output 5: PCA score plot
            #********************************************************************#
            
            #Perform PCA
            rv$PCA_data <- prcomp(t(rv$normMatrix[apply(rv$normMatrix, 1, var) != 0,]),
                                  retx = TRUE, 
                                  center = TRUE,
                                  scale.=TRUE)
            
            
            # Make PCA plot
            output$PCA_microarray_raw <- renderPlotly({
              
              # Get vector for colouring the samples
              if(length(input$colorFactor_microarray_raw) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_microarray_raw], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$colorFactor_microarray_raw])
              }
              
              # Make PCA plot
              plot_PCA(PC_data = rv$PCA_data, 
                       colorFactor = colorFactor, 
                       xpc = as.numeric(str_remove(input$xpca_microarray_raw,"PC")), 
                       ypc = as.numeric(str_remove(input$ypca_microarray_raw,"PC")), 
                       zpc = ifelse(input$xyz_microarray_raw,as.numeric(str_remove(input$zpca_microarray_raw,"PC")),3), 
                       xyz = input$xyz_microarray_raw)
              
            })
            
            #********************************************************************#
            # Output 6: Overview of pre-processing settings
            #********************************************************************#
            # Print table with settings
            output$processingSettings_microarray_raw <- DT::renderDataTable({
              return(rv$processingSettings)
            },options = list(pageLength = 10),
            selection = "none")
            
            # Download button
            output$downloadProcessingSettings_microarray_raw <- downloadHandler(
              filename = "preprocessingSettings.csv",
              content = function(file){
                write.csv(rv$processingSettings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
            
            #------------------------------------------------------------------#
            # Display output: QC tables/plots
            #------------------------------------------------------------------#
            
            
            # Render UI for main tab
            output$UI_QC_microarray_raw <- renderUI({
              tagList(
                tabsetPanel(
                  
                  # Expression table
                  tabPanel("Expression values",
                           icon = icon("fas fa-mouse-pointer"),
                           h3(strong("Normalized expression values")),
                           h5("Here you can view the normalized log intensity (expression) values."),
                           hr(),
                           DT::dataTableOutput(outputId = "exprTable_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadNormalizedData_microarray_raw", 
                                          "Download"),
                           br(),
                           br(),
                           
                           # Dropdown Button to adjust the plot settings
                           shinyWidgets::dropdownButton(
                             tags$div(
                               style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                               
                               
                               # Change order of the boxplots:
                               tags$h4("Drag to change boxplot order"),
                               shinyjqui::orderInput(inputId = 'geneboxplot_order_microarray_raw', 
                                                     label = NULL, 
                                                     items = levels(rv$experimentFactor),
                                                     item_class = 'default'),
                               br(),
                               
                               # Change colour of the boxplots by button. 
                               # This is used when there are more than 6 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_microarray_raw > 6",
                                 tags$h4("Click to change boxplot colours"),
                                 shinyWidgets::actionBttn("geneboxplot_changeOrder_microarray_raw",
                                                          label = "Change color",
                                                          style = "simple",
                                                          color = "primary",
                                                          icon = icon("sync"))
                               ),
                               
                               # Change colour of the boxplots by colour picker
                               # This is used when there are less than 7 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_microarray_raw < 7",
                                 tags$h4("Click to select boxplot colours"),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 0",
                                   colourpicker::colourInput("geneboxplot_col1_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[1])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 1",
                                   colourpicker::colourInput("geneboxplot_col2_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[2])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 2",
                                   colourpicker::colourInput("geneboxplot_col3_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[3])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 3",
                                   colourpicker::colourInput("geneboxplot_col4_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[4])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 4",
                                   colourpicker::colourInput("geneboxplot_col5_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[5])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_raw > 5",
                                   colourpicker::colourInput("geneboxplot_col6_microarray_raw", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[6])
                                 )
                               ),
                               br(),
                               tags$h4("Drag to change jitter"),
                               sliderInput("jitter_geneboxplot_microarray_raw", 
                                           NULL,
                                           min = 0, max = 0.3,
                                           value = 0.1, step = 0.01),
                               br()
                             ),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"),
                             
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to personalize plot!")
                             
                           ), # EO dropdownButton
                           
                           # Boxplot of the selected gene's expression values
                           plotOutput("ExprBoxplot_microarray_raw")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1"),
                           
                           actionButton("download_geneboxplot_microarray_raw", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           br(),
                           br()
                  ),
                  
                  # Boxplots
                  tabPanel("Boxplots",
                           icon = icon("fas fa-file"),
                           br(),
                           actionButton("download_boxplots_microarray_raw", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           plotOutput(outputId = "boxplots_microarray_raw",
                                      width = "65vw", height = "80vw")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # Density plot
                  tabPanel("Density plots",
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           h2(strong("Density plot of normalized intensities"), align = "center"),
                           h4("Distributions should be comparable between arrays", align = "center"),
                           plotly::plotlyOutput(outputId = "densityplots_microarray_raw",
                                                width = "65vw", height = "40vw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # Sample-sample correlation heatmap
                  tabPanel("Correlation plot", 
                           icon = icon("fas fa-mouse-pointer"),
                           
                           br(),
                           
                           fluidRow(
                             
                             #Distance method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption1_microarray_raw",
                                      label = "Distance calculation 
                                              method",
                                      choices = c("Pearson","Spearman",
                                                  "Euclidean"),
                                      options = list(
                                        style = "btn-primary"))
                             ),
                             
                             #Ckustering method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption2_microarray_raw",
                                      label = "Clustering method",
                                      choices = c("ward.D2","single",
                                                  "complete","average",
                                                  "mcquitty","median",
                                                  "centroid"),
                                      options = list(
                                        style = "btn-info"))
                             )),
                           
                           hr(),
                           
                           #Theme
                           dropdownButton(
                             tags$h3("Theme"),
                             
                             selectInput(inputId = "heatmapFactor_microarray_raw",
                                         label = "Side colors",
                                         choices = colnames(rv$metaData_fil),
                                         selected = rv$experimentName,
                                         multiple = TRUE),
                             
                             selectInput(inputId = 'heatmaptheme_microarray_raw',
                                         label = "Heatmap theme",
                                         choices = c("Default", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"), width = "300px",
                             
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to change colors!")
                             
                           ),
                           
                           # Heatmap
                           plotly::plotlyOutput("heatmap_microarray_raw", 
                                                width = "1000px", 
                                                height="600px") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1", 
                                                          proxy.height = "400px")
                           
                  ),
                  
                  # PCA score plot
                  tabPanel("PCA",
                           icon = icon("fas fa-mouse-pointer"),
                           # Title + text
                           fluidRow(
                             h3(strong("Principal Component Analysis (PCA)")),
                             h5("Here you can view the PCA score plot."),
                             hr(),
                           ),
                           
                           # Set color + 3D/2D
                           fluidRow(
                             column(3,
                                    # Color by which factor?
                                    shinyWidgets::pickerInput(inputId = "colorFactor_microarray_raw",
                                                              label = "Color by:",
                                                              choices = colnames(rv$metaData_fil),
                                                              selected = rv$experimentName,
                                                              multiple = TRUE)
                             ),
                             column(3,
                                    br(),
                                    # 3D plot?
                                    shinyWidgets::materialSwitch(
                                      inputId = "xyz_microarray_raw",
                                      label = "3D",
                                      value = FALSE, 
                                      status = "danger")
                                    
                             )
                           ),
                           
                           # Set axes
                           fluidRow(
                             column(3,
                                    #X-axis
                                    selectInput(inputId = "xpca_microarray_raw", 
                                                label = "x-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC1")
                             ),
                             column(3,
                                    #Y-axis
                                    selectInput(inputId = "ypca_microarray_raw", 
                                                label = "y-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC2")
                             ),
                             column(3,
                                    #Z-axis
                                    conditionalPanel(
                                      condition = "input.xyz_microarray_raw==true",
                                      selectInput(inputId = "zpca_microarray_raw", 
                                                  label = "z-axis",
                                                  choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                              "PC6", "PC7", "PC8"),
                                                  selected = "PC3")
                                    )
                             )
                           ),
                           
                           # Print plot
                           fluidRow(
                             hr(),
                             plotly::plotlyOutput("PCA_microarray_raw")%>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")
                           ) 
                           
                  ), 
                  # Settings table
                  tabPanel("Settings overview",
                           icon = icon("fas fa-file"),
                           h3(strong("Pre-processing settings")),
                           h5("Here you can see an overview of the chosen pre-processing settings."),
                           hr(),
                           DT::dataTableOutput(outputId = "processingSettings_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadProcessingSettings_microarray_raw", 
                                          "Download"),
                           
                  ) # EO Settings tabPanel
                ) # EO tabsetPanel
              ) # EO tagList
            }) # EO renderUI
            
            # Allow user to go to next tab
            output$UI_next_preprocessing_microarray_raw <- renderUI({
              req(rv$normData)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                shinyWidgets::actionBttn(inputId = "next_preprocessing_microarray_raw",
                                         label = "Next",
                                         style = "jelly",
                                         color = "danger",
                                         icon = icon("arrow-right"))
              ) # EO tagList
            }) # EO renderUI
          }) # EO observeEvent
          
          #====================================================================#
          # TAB 3: Statistical analysis
          #====================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_preprocessing_microarray_raw,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_statistics_microarray_raw")
          })
          
          # SIDEPANEL:
          
          # Select continuous covariates
          output$UI_covGroups_num_microarray_raw <- renderUI({
            tagList(
              shinyWidgets::pickerInput(inputId = "covGroups_num_microarray_raw",
                                        label = "Continuous covariates (e.g., age):",
                                        choices = setdiff(colnames(rv$metaData_fil),
                                                          rv$experimentName),
                                        selected = NULL,
                                        multiple = TRUE)
            )
          })
          
          # Select discrete covariates
          output$UI_covGroups_char_microarray_raw <- renderUI({
            tagList(
              shinyWidgets::pickerInput(inputId = "covGroups_char_microarray_raw",
                                        label = "Discrete covariates (e.g., sex):",
                                        choices = setdiff(colnames(rv$metaData_fil),
                                                          rv$experimentName),
                                        selected = NULL,
                                        multiple = TRUE)
            )
          })
          
          # Select comparisons
          output$UI_comparisons_microarray_raw <- renderUI({
            tagList(
              shinyWidgets::multiInput(
                inputId = "comparisons_microarray_raw",
                label = "Comparisons:", 
                choices = makeComparisons(make.names(unique(rv$experimentFactor))),
                selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1]
              )
            )
          })
          
          
          
          # Select comparisons
          output$UI_biomart_dataset_microarray_raw <- renderUI({
            req(input$addAnnotation_microarray_raw)
            shinyWidgets::pickerInput(inputId = "biomart_dataset_microarray_raw",
                                      label = "Organism:",
                                      choices = c("Homo sapiens" = "hsapiens_gene_ensembl" ,
                                                  "Bos taurus" = "btaurus_gene_ensembl",
                                                  "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                                  "Mus musculus" = "mmusculus_gene_ensembl",
                                                  "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                                      selected = switch(rv$Organism,
                                                        "Homo sapiens" = "hsapiens_gene_ensembl" ,
                                                        "Bos taurus" = "btaurus_gene_ensembl",
                                                        "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                                        "Mus musculus" = "mmusculus_gene_ensembl",
                                                        "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                                      multiple = FALSE)
          })
          
          output$UI_addAnnotations_microarray_raw <- renderUI({
            req(input$addAnnotation_microarray_raw)
            req(input$biomart_dataset_microarray_raw)
            
            tagList(
              
              shinyWidgets::pickerInput(inputId = "biomart_filter_microarray_raw",
                                        label = "Probeset ID",
                                        choices = filterList[[input$biomart_dataset_microarray_raw]],
                                        selected = selFilter(rv$ProbeAnnotation),
                                        multiple = FALSE),
              
              shinyWidgets::pickerInput(inputId = "biomart_attributes_microarray_raw",
                                        label = "Output",
                                        choices = c("Ensembl Gene ID",
                                                    "Entrez Gene ID",
                                                    "Gene Symbol/Name"),
                                        selected = "Gene Symbol/Name",
                                        multiple = TRUE)
            )
            
          })
          
          
          
          #======================================================================#
          # Output of statistical analysis
          #======================================================================#
          observeEvent(input$calculate_statistics_microarray_raw,{
            shinybusy::show_modal_spinner(text = "Statistical analysis...",
                                          color="#0dc5c1")
            
            # Calculate statistics
            if (isTRUE(input$addAnnotation_microarray_raw)){
              rv$top_table_list <- getStatistics(normMatrix = rv$normMatrix, 
                                            metaData = rv$metaData_fil, 
                                            expFactor = rv$experimentName, 
                                            covGroups_num = input$covGroups_num_microarray_raw,
                                            covGroups_char = input$covGroups_char_microarray_raw,
                                            comparisons = input$comparisons_microarray_raw,
                                            addAnnotation = input$addAnnotation_microarray_raw,
                                            biomart_dataset = input$biomart_dataset_microarray_raw,
                                            biomart_attributes = unique(c(input$biomart_filter_microarray_raw,
                                                                          input$biomart_attributes_microarray_raw)),
                                            biomart_filters = input$biomart_filter_microarray_raw)
              
            } else{
              rv$top_table_list <- getStatistics(normMatrix = rv$normMatrix, 
                                            metaData = rv$metaData_fil, 
                                            expFactor = rv$experimentName, 
                                            covGroups_num = input$covGroups_num_microarray_raw,
                                            covGroups_char = input$covGroups_char_microarray_raw,
                                            comparisons = input$comparisons_microarray_raw,
                                            addAnnotation = input$addAnnotation_microarray_raw,
                                            biomart_dataset = NULL,
                                            biomart_attributes = NULL,
                                            biomart_filters = NULL)
              
            }
            
            # Select comparison for output
            observe({
              req(rv$normMatrix)
              req(rv$experimentFactor)
              
              if (!is.null(rv$top_table_list)){
                rv$top_table <- rv$top_table_list[[1]]
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show comparisons
                output$UI_comparisons_view_microarray_raw <- renderUI({
                  shinyWidgets::pickerInput(inputId = "comparisons_view_microarray_raw",
                                            label = "Select comparison:",
                                            choices = names(rv$top_table),
                                            selected = names(rv$top_table)[1],
                                            multiple = FALSE)
                })
                
                # Show message
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = rv$top_table_list[[2]],
                  type = "info")
                
                # Show microarray statistics tab
                showTab("navbar", target = "panel_ORA_microarray_raw")
                
              } else{
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show message
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = "Oops...something went wrong! Please try again!",
                  type = "error")
              }
              
            })
            
            #--------------------------#
            # Generate output figures:
            #--------------------------#
            
            #********************************************************************#
            # top table
            #********************************************************************#
            
            # print top table
            output$top_table_microarray_raw <- DT::renderDataTable({
              
              # Required: selected comparison
              req(input$comparisons_view_microarray_raw)
              
              # Get statistics of selected comparison
              output <- rv$top_table[[input$comparisons_view_microarray_raw]]
              
              # Add link to ENSEMBL or NCBI website if custom annotations is selected
              if (!is.null(rv$ProbeAnnotation)){
                
                # If custom annotation = ENTREZ, link to NCBI website
                if(rv$ProbeAnnotation == "ENTREZG"){
                  output$GeneID <- paste0(
                    '<a ',
                    'href=',
                    paste(
                      "https://www.ncbi.nlm.nih.gov/gene/",
                      output$GeneID,
                      sep = ''
                    ),
                    ' target="_blank"',
                    '>',
                    output$GeneID,
                    '</a>'
                  )
                }
                
                # If custom annotation = ENSG, link to ENSEMBL website
                if(rv$ProbeAnnotation == "ENSG"){
                  output$GeneID <- paste0(
                    '<a ',
                    'href=',
                    paste(
                      "http://www.ensembl.org/id/",
                      output$GeneID,
                      sep = ''
                    ),
                    ' target="_blank"',
                    '>',
                    output$GeneID,
                    '</a>'
                  )
                }
              }
              return(output)
              
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$download_top_table_microarray_raw <- downloadHandler(
              filename = paste0("topTable_",input$comparisons_view_microarray_raw,".csv"),
              content = function(file){
                write.csv(rv$top_table[[input$comparisons_view_microarray_raw]], file, quote = FALSE, row.names = FALSE)
              }
            )
            
            # Change plotting data depending on whether all experimental groups will be plotted  
            observe({
              req(input$comparisons_view_microarray_raw)
              req(rv$top_table)
              
              if (!is.null(input$boxplotAll_microarray_raw)){
                if (input$boxplotAll_microarray_raw){
                  rv$newFactor <- rv$experimentFactor
                  rv$newData <- rv$normMatrix
                }
                if (!input$boxplotAll_microarray_raw){
                  if(length(rv$experimentName) > 1){
                    t <- make.names(apply(rv$metaData_fil[,rv$experimentName], 1, paste, collapse = "_" ))
                  } else{
                    t <- make.names(rv$metaData_fil[,rv$experimentName])
                  }
                  rv$newData <- rv$normMatrix[,t %in% (unlist(stringr::str_split(input$comparisons_view_microarray_raw, " - ")))]
                  rv$newFactor <- factor(as.character(rv$experimentFactor)[t %in% (unlist(stringr::str_split(input$comparisons_view_microarray_raw, " - ")))])
                }
              }
              
              # Change color by click on button
              rv$colorOrder <- 1:length(levels(rv$newFactor))
              observeEvent(input$statboxplot_changeOrder_microarray_raw,{
                all_orders <- permute(1:length(levels(rv$newFactor))) 
                sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
                
                if (sel < length(all_orders)){
                  rv$colorOrder <- all_orders[[sel+1]]
                } else{
                  rv$colorOrder <- all_orders[[1]]
                }
              })
              
              
              # Boxplot of single gene (based on selected row in the top table)
              output$ExprBoxplot_statistics_microarray_raw <- renderPlot({
                req(input$top_table_microarray_raw_rows_selected)
                req(rv$top_table)
                
                if (length(levels(rv$newFactor)) > 6){
                  legendColors <- colorsByFactor(rv$newFactor)$legendColors
                } else{
                  legendColors <- c(input$statboxplot_col1_microarray_raw,
                                    input$statboxplot_col2_microarray_raw,
                                    input$statboxplot_col3_microarray_raw,
                                    input$statboxplot_col4_microarray_raw,
                                    input$statboxplot_col5_microarray_raw,
                                    input$statboxplot_col6_microarray_raw)
                }
                
                gene <- rv$top_table[[input$comparisons_view_microarray_raw]]$GeneID[input$top_table_microarray_raw_rows_selected]
                sel_row <- which(as.character(rownames(rv$normMatrix)) %in% as.character(gene))
                
                # Make boxplot
                rv$temp <- geneBoxplot(experimentFactor = rv$newFactor, 
                                       normMatrix = rv$newData, 
                                       sel_row = sel_row,
                                       legendColors = legendColors[rv$colorOrder],
                                       groupOrder = input$statboxplot_order_microarray_raw,
                                       jitter = input$jitter_statboxplot_microarray_raw,
                                       rnaseq=FALSE,
                                       seed = sample(1:1000,1))
                
                return(rv$temp)
                
              })
              
              # Get number of experimental groups
              output$length_statboxplot_microarray_raw <- reactive({
                length(levels(rv$experimentFactor))
              })
              outputOptions(output, "length_statboxplot_microarray_raw", suspendWhenHidden = FALSE) 
              
            })
            
            #***************************#
            # Modal to download boxplot
            #***************************#

            # Download plot
            output$realdownload_statboxplot_microarray_raw <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp, 
                                filename = file,
                                width = input$width_statboxplot_microarray_raw,
                                height = input$height_statboxplot_microarray_raw,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_statboxplot_microarray_raw, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_statboxplot_microarray_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_statboxplot_microarray_raw", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_statboxplot_microarray_raw', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Histograms
            #********************************************************************#
            observe({
              req(input$comparisons_view_microarray_raw)
              
              if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
                # Make P value histrogram
                output$Phistogram_microarray_raw <- renderPlotly({
                  req(rv$top_table)
                  p <- makePHistogram(rv$top_table[[input$comparisons_view_microarray_raw]][,"p-value"])
                  return(p)
                })
                
                # Make log2FC histogram
                output$logFChistogram_microarray_raw <- renderPlotly({
                  req(rv$top_table)
                  p <- makelogFCHistogram(rv$top_table[[input$comparisons_view_microarray_raw]][,"log2FC"])
                  return(p)
                })
                
              }
            })
            
            #********************************************************************#
            # Volcano plot
            #********************************************************************#
            
            # Re-render the volcano plot after pressing the "plot" button
            # ignoreNULL = FALSE, so it will render the plot the first time 
            # without needing to press the "plot" button.
            observeEvent(input$plot_volcano_microarray_raw, {
              
              # Requirements
              req(rv$top_table) # Top table
              req(input$rawp_volcano_microarray_raw) # raw or adj. P value?
              req(input$p_thres_volcano_microarray_raw) # P value threshold?
              req(input$logFC_thres_volcano_microarray_raw) # log2FC threshold?
              req(input$comparisons_view_microarray_raw)
              
              if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
                
                # Make plot
                p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                                 p = input$rawp_volcano_microarray_raw, 
                                 p_threshold = input$p_thres_volcano_microarray_raw, 
                                 logFC_threshold = input$logFC_thres_volcano_microarray_raw)
                
                # Render plot
                output$volcano_microarray_raw <- renderPlotly(p)
              }
              
            }, ignoreNULL = FALSE) 
            
            
            #--------------------------------------#
            # dynamic UI: Output in different tabs
            #--------------------------------------#
            observe({
              req(rv$top_table)
              
              output$UI_boxplotAll_microarray_raw <- renderUI({
                tagList(
                  # Change order of the boxplots:
                  tags$h4("Drag to change boxplot order"),
                  shinyjqui::orderInput(inputId = 'statboxplot_order_microarray_raw', 
                                        label = NULL, 
                                        items = levels(rv$newFactor),
                                        item_class = 'default'),
                  br(),
                  
                  # Change colour of the boxplots by button. 
                  # This is used when there are more than 6 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_microarray_raw > 6",
                    tags$h4("Click to change boxplot colours"),
                    shinyWidgets::actionBttn("statboxplot_changeOrder_microarray_raw",
                                             label = "Change color",
                                             style = "simple",
                                             color = "primary",
                                             icon = icon("sync"))
                  ),
                  
                  # Change colour of the boxplots by colour picker
                  # This is used when there are less than 7 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_microarray_raw < 7",
                    tags$h4("Click to select boxplot colours"),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 0",
                      colourpicker::colourInput("statboxplot_col1_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[1])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 1",
                      colourpicker::colourInput("statboxplot_col2_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[2])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 2",
                      colourpicker::colourInput("statboxplot_col3_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[3])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 3",
                      colourpicker::colourInput("statboxplot_col4_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[4])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 4",
                      colourpicker::colourInput("statboxplot_col5_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[5])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_raw > 5",
                      colourpicker::colourInput("statboxplot_col6_microarray_raw", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[6])
                    )
                  ),
                  br(),
                  tags$h4("Drag to change jitter"),
                  sliderInput("jitter_statboxplot_microarray_raw", 
                              NULL,
                              min = 0, max = 0.3,
                              value = 0.1, step = 0.01),
                  br()
                )
              })
            })
            
            observe({
              if (is.null(rv$top_table)){
                output$UI_output_statistics_microarray_raw <- renderUI(NULL)
              } else{
                output$UI_output_statistics_microarray_raw <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #********************************************************************#
                      # top table tab
                      #********************************************************************#
                      
                      tabPanel("Top table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               # Title + description
                               h3(strong("Top Table")),
                               h5("The Top Table includes the output of the selected statistical analysis."),
                               hr(),
                               
                               # Top table
                               dataTableOutput(outputId = "top_table_microarray_raw") %>% 
                                 withSpinner(color="#0dc5c1"),
                               
                               # Download button
                               downloadButton("download_top_table_microarray_raw", 
                                              "Download table"),
                               br(),
                               br(),
                               
                               # Dropdown Button to adjust the plot settings
                               shinyWidgets::dropdownButton(
                                 tags$div(
                                   style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                                   
                                   # Plot all experimental groups?
                                   conditionalPanel(
                                     condition = "output.length_statboxplot_microarray_raw > 2",
                                     tags$h4("Plot all experimental groups?"),
                                     shinyWidgets::materialSwitch(inputId = "boxplotAll_microarray_raw",
                                                                  label = NULL, 
                                                                  value = TRUE,
                                                                  status = "primary"),
                                     
                                     br()
                                   ),
                                   uiOutput("UI_boxplotAll_microarray_raw")
                                   
                                   
                                   
                                 ),
                                 circle = TRUE, status = "info",
                                 icon = icon("fas fa-cog"),
                                 
                                 tooltip = shinyWidgets::tooltipOptions(
                                   title = "Click to personalize plot!")
                                 
                               ), # EO dropdownButton
                               plotOutput("ExprBoxplot_statistics_microarray_raw")%>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1"),
                               actionButton("download_statboxplot_microarray_raw", 
                                            "Download figure",
                                            icon = shiny::icon("download")),
                               br(),
                               br()
                      ),
                      
                      #********************************************************************#
                      # histogram tab
                      #********************************************************************#
                      
                      tabPanel("Histograms",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               # P value histogram
                               plotlyOutput("Phistogram_microarray_raw")%>% 
                                 withSpinner(color="#0dc5c1"),
                               hr(),
                               
                               # log2FC histogram
                               plotlyOutput("logFChistogram_microarray_raw")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      ),
                      
                      #********************************************************************#
                      # volcano tab
                      #********************************************************************#
                      
                      tabPanel("Volcano plot",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               fluidRow(
                                 column(3,
                                        # Raw or adjusted P value?
                                        prettyRadioButtons(
                                          inputId = "rawp_volcano_microarray_raw",
                                          label = "P value", 
                                          choices = 
                                            c("Raw P value" = "raw", 
                                              "Adjusted P value" = "adj"))
                                 ),
                                 column(3,
                                        # P value threshold
                                        numericInput(
                                          inputId = "p_thres_volcano_microarray_raw",
                                          label = "P threshold",
                                          value = 0.05)
                                 ),
                                 column(3,
                                        # log2FC threshold
                                        numericInput(
                                          inputId = "logFC_thres_volcano_microarray_raw",
                                          label = "logFC threshold",
                                          value = 1)
                                 )
                               ),
                               hr(),
                               
                               # Action button: press to (re-)render plot
                               actionBttn(inputId = "plot_volcano_microarray_raw", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               
                               # Volcano plot
                               plotlyOutput("volcano_microarray_raw")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      )
                      
                    ) # tabsetpanel
                  ) # taglist
                }) # renderUI
              }
            }) # Observe
            
            # Allow user to go to next tab
            output$UI_next_statistics_microarray_raw <- renderUI({
              req(rv$top_table)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                actionBttn(inputId = "next_statistics_microarray_raw",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          #======================================================================#
          # ORA
          #======================================================================#
          # Go to data upload tab
          observeEvent(input$next_statistics_microarray_raw,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_ORA_microarray_raw")
          })
          
          # Comparisons for which ORA should be performed
          observe({
            req(rv$top_table)
            output$UI_comparisons_view_ORA_microarray_raw <- renderUI({
              shinyWidgets::pickerInput(inputId = "comparisons_view_ORA_microarray_raw",
                                        label = NULL,
                                        choices = names(rv$top_table),
                                        selected = names(rv$top_table)[1],
                                        multiple = FALSE)
            })
          })
          
          # Which columns contain gene ids
          observe({
            
            # Require: selected compairson
            req(input$comparisons_view_ORA_microarray_raw)
            
            # Find which columns can contain gene IDs:
            # - The first column (i.e., probeset ids) 
            # - If present: the columns 7 to end (i.e., biomaRt annotation)
            col_choice <- 1
            if (length(colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])) > 6){
              col_choice <- c(1,7:ncol(rv$top_table[[input$comparisons_view_ORA_microarray_raw]]))
            }
            
            # Render UI
            output$UI_geneID_ORA_microarray_raw <- renderUI({
              tagList(
                
                # Select organism (needed for GO/WP annotations)
                selectInput(inputId = "organism_ORA_microarray_raw",
                            label = "Organism",
                            choices = c("Bos taurus",
                                        "Caenorhabditis elegans",
                                        "Homo sapiens",
                                        "Mus musculus", 
                                        "Rattus norvegicus"),
                            selected = rv$Organism),
                
                # Select column with gene IDs
                shinyWidgets::pickerInput(inputId = "geneID_ORA_microarray_raw",
                                          label = "Which column of the top table contains the gene IDs?",
                                          choices = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])[col_choice],
                                          selected = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])[1],
                                          multiple = FALSE),
                
                shinyWidgets::pickerInput(inputId = "selID_ORA_microarray_raw",
                                          label = "Which gene ID to use?",
                                          choices = c("Ensembl Gene ID" = "ENSEMBL", 
                                                      "Entrez Gene ID" = "ENTREZID", 
                                                      "Gene Symbol/Name" = "SYMBOL"),
                                          selected = selID(rv$ProbeAnnotation),
                                          multiple = FALSE),
              )
            })
          })
          
          # Show modal
          observeEvent(input$calculate_ORA_microarray_raw,{
            shinybusy::show_modal_spinner(text = "Overrepresentation analysis...",
                                          color="#0dc5c1")
            
            # Perform ORA
            if (input$topNorThres_microarray_raw == "Threshold"){
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                 geneset = input$geneset_ORA_microarray_raw,
                                 geneID_col = input$geneID_ORA_microarray_raw,
                                 geneID_type = input$selID_ORA_microarray_raw,
                                 organism = input$organism_ORA_microarray_raw,
                                 updown = input$updown_ORA_microarray_raw,
                                 topN = FALSE,
                                 N = NULL,
                                 rawadj = input$rawp_ORA_microarray_raw,
                                 p_thres = input$p_thres_ORA_microarray_raw,
                                 logFC_thres = input$logFC_thres_ORA_microarray_raw)
            } else{
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                 geneset = input$geneset_ORA_microarray_raw,
                                 geneID_col = input$geneID_ORA_microarray_raw,
                                 geneID_type = input$selID_ORA_microarray_raw,
                                 organism = input$organism_ORA_microarray_raw,
                                 updown = input$updown_ORA_microarray_raw,
                                 topN = TRUE,
                                 N = input$topN_microarray_raw,
                                 rawadj = NULL,
                                 p_thres = NULL,
                                 logFC_thres = NULL)
            }
            
            observe({
              
              # Remove modal
              remove_modal_spinner()
              
              # Show message
              if (is.null(rv$ORA_data)){
                sendSweetAlert(
                  session = session,
                  title = "Error!",
                  text = "No significant genes!",
                  type = "error")
                
              }else{
                
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = "Gene overrepresentation analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
                  type = "info")
                
                # print GO table
                output$ORA_table_microarray_raw <- renderDataTable({
                  req(input$geneset_ORA_microarray_raw)
                  output <- rv$ORA_data@result
                  
                  if (input$geneset_ORA_microarray_raw == "WikiPathways"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.wikipathways.org/pathways/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_microarray_raw == "KEGG"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.genome.jp/pathway/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_microarray_raw %in% c("GO-BP", "GO-MF", "GO-CC")){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.ebi.ac.uk/QuickGO/term/",
                        output$ID
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  
                  return(output)
                },options = list(pageLength = 6),
                selection = list(mode = "single", selected = 1), escape = FALSE)
                
                # Download button
                output$download_ORA_table_microarray_raw <- downloadHandler(
                  filename = paste0("ORATable_",input$comparisons_view_ORA_microarray_raw,"_",input$geneset_ORA_microarray_raw,".csv"),
                  content = function(file){
                    write.csv(rv$ORA_data@result, file, quote = FALSE, row.names = FALSE)
                  }
                )
                
                # Print statistics of genes in selected Term
                output$ORAgene_table_microarray_raw <- renderDataTable({
                  req(input$ORA_table_microarray_raw_rows_selected)
                  
                  # Make ORA gene table
                  output <- make_ORAgene_table(ORA_data = rv$ORA_data,
                                               top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                               geneID_col = input$geneID_ORA_microarray_raw,
                                               sel_row_ORA = input$ORA_table_microarray_raw_rows_selected)
                  
                  # Add links to gene ids
                  if (!is.null(rv$ProbeAnnotation)){
                    if(rv$ProbeAnnotation == "ENTREZG"){
                      output$GeneID <- paste0(
                        '<a ',
                        'href=',
                        paste(
                          "https://www.ncbi.nlm.nih.gov/gene/",
                          output$GeneID,
                          sep = ''
                        ),
                        ' target="_blank"',
                        '>',
                        output$GeneID,
                        '</a>'
                      )
                    }
                    
                    if(rv$ProbeAnnotation == "ENSG"){
                      output$GeneID <- paste0(
                        '<a ',
                        'href=',
                        paste(
                          "http://www.ensembl.org/id/",
                          output$GeneID,
                          sep = ''
                        ),
                        ' target="_blank"',
                        '>',
                        output$GeneID,
                        '</a>'
                      )
                    }
                  }
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  output$meanExpr <- round(output$meanExpr,3)
                  output$log2FC <- round(output$log2FC,3)
                  output$`log2FC SE` <- round(output$`log2FC SE`,3)
                  
                  return(output)
                }, options = list(pageLength = 6), escape = FALSE)
                
                # Text for gene table
                output$text_ORAgene_table_microarray_raw <- renderText({
                  text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_microarray_raw_rows_selected,"ID"],
                                 "</b></h3>")
                  return(text)
                })
                
                # ORA barchart
                observeEvent(input$plot_ORAplot_microarray_raw, {
                  req(input$nSets_ORAplot_microarray_raw)
                  p <- makeORAplot(rv$ORA_data,
                                   nSets = input$nSets_ORAplot_microarray_raw)
                  
                  output$ORAplot_microarray_raw <- renderPlotly(p)
                  
                }, ignoreNULL = FALSE)
                
                # ORA network diagram
                observeEvent(input$plot_ORAnetwork_microarray_raw, {
                  req(input$layout_network_microarray_raw)
                  req(input$nSets_network_microarray_raw)
                  p <- makeORAnetwork(ORA_data = rv$ORA_data,
                                      layout = input$layout_network_microarray_raw,
                                      nSets = input$nSets_network_microarray_raw)
                  
                  output$ORAnetwork_microarray_raw <- renderPlot(p,
                                                                 height = 500, 
                                                                 width = 800)
                  
                }, ignoreNULL = FALSE)
                
              }
            }) # Observe
            
            
            
            # Render output table
            observe({
              output$UI_output_ORA_microarray_raw <- renderUI({
                tagList(
                  tabsetPanel(
                    
                    #********************************************************************#
                    # top table tab
                    #********************************************************************#
                    
                    tabPanel("ORA table",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             h3(strong("Statistics table")),
                             h5("The ORA statistics table encompasses the output of the gene 
                              overrepresentation analysis."),
                             hr(),
                             dataTableOutput(outputId = "ORA_table_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1"),
                             downloadButton("download_ORA_table_microarray_raw", 
                                            "Download"),
                             br(),
                             htmlOutput("text_ORAgene_table_microarray_raw"),
                             h5(paste0("The gene table encompasses the statistics of all genes 
                              from the selected geneset.")),
                             hr(),
                             dataTableOutput(outputId = "ORAgene_table_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")
                    ),
                    
                    tabPanel("Bar chart",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             h3(strong("ORA bar chart")),
                             h5("The ORA bar chart visualizes the results from the overrepresentation analysis."),
                             hr(),
                             fluidRow(
                               column(3,
                                      # Number of genesets
                                      numericInput(
                                        inputId = "nSets_ORAplot_microarray_raw",
                                        label = "Number of genesets (5-20)",
                                        value = 10,
                                        min = 5,
                                        max = 20,
                                        step = 1)
                               )
                             ),
                             hr(),
                             actionBttn(inputId = "plot_ORAplot_microarray_raw", 
                                        label = "Plot",
                                        style = "jelly",
                                        color = "primary",
                                        icon = icon("sync")),
                             br(),
                             br(),
                             plotlyOutput("ORAplot_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")
                    ),
                    
                    tabPanel("Network diagram",
                             icon = icon("fas fa-mouse-pointer"),
                             br(),
                             h3(strong("ORA network diagram")),
                             h5("The ORA network diagram visualize the similarity between the most significant genesets."),
                             hr(),
                             fluidRow(
                               column(3,
                                      # Number of genesets
                                      numericInput(
                                        inputId = "nSets_network_microarray_raw",
                                        label = "Number of genesets (5-20)",
                                        value = 10,
                                        min = 5,
                                        max = 20,
                                        step = 1)
                               ),
                               column(3,
                                      # Network layout
                                      pickerInput(inputId = "layout_network_microarray_raw",
                                                  label = "Network layout",
                                                  choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                              'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                  selected = 'graphopt',
                                                  multiple = FALSE)
                               )
                             ),
                             hr(),
                             actionBttn(inputId = "plot_ORAnetwork_microarray_raw", 
                                        label = "Plot",
                                        style = "jelly",
                                        color = "primary",
                                        icon = icon("sync")),
                             br(),
                             br(),
                             plotOutput("ORAnetwork_microarray_raw") %>% 
                               withSpinner(color="#0dc5c1")
                    ),
                    
                  ) # tabsetpanel
                ) # taglist
              }) # renderUI
            }) # Observe
            
            
            
          }) #observe event
          
        } # raw or norm
        
        
        #************************************************************************#
        # Preprocessed microarray data
        #************************************************************************#
        if (raw_or_norm() == "Processed data"){
          
          
          #======================================================================#
          # Data Upload
          #======================================================================#
          
          #----------------------------------------------------------------------#
          # Go to data upload tab
          #----------------------------------------------------------------------#
          
          observeEvent(input$startAnalysis,{
            
            # Show microarray (norm) upload tab
            showTab("navbar", target = "panel_upload_microarray_norm")
            
            # Remove the other RNA-seq (raw) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_raw")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_raw")
            hideTab("navbar", target = "panel_statistics_rnaseq_raw" )
            hideTab("navbar", target = "panel_ORA_rnaseq_raw")
            
            # Remove the other RNA-seq (norm) tabs
            hideTab("navbar", target = "panel_upload_rnaseq_norm")
            hideTab("navbar", target = "panel_preprocessing_rnaseq_norm")
            hideTab("navbar", target = "panel_statistics_rnaseq_norm")
            hideTab("navbar", target = "panel_ORA_rnaseq_norm")
            
            # Remove the other microarray (raw) tabs
            hideTab("navbar", target = "panel_upload_microarray_raw")
            hideTab("navbar", target = "panel_preprocessing_microarray_raw")
            hideTab("navbar", target = "panel_statistics_microarray_raw" )
            hideTab("navbar", target = "panel_ORA_microarray_raw")
            
            # Remove the other microarray (norm) tabs
            hideTab("navbar", target = "panel_preprocessing_microarray_norm")
            hideTab("navbar", target = "panel_statistics_microarray_norm")
            hideTab("navbar", target = "panel_ORA_microarray_norm")
            
            # Go to microarray (norm) tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_upload_microarray_norm")
            
            # Example meta data file
            output$downloadmeta_example_microarray_norm <- downloadHandler(
              filename = "MetaData_example.csv",
              content = function(file){
                write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
              }
            )
          })
          
          #----------------------------------------------------------------------#
          # Upload files
          #----------------------------------------------------------------------#
          
          observeEvent(input$upload_microarray_norm,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Read expression data
            if (!is.null(input$uploadExprData_microarray_norm_smf)){
              rv$gxData <- getGEO(filename=input$uploadExprData_microarray_norm_smf$datapath)
              # Guess the organism
              rv$Organism <- getOrganism(gxData = rv$gxData)
            } else{
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "You forgot to upload an expression and/or metadata file",
                type = "error")
              shinybusy::remove_modal_spinner()
            }
            
            # Get metadata
            if (input$MetaFileType_microarray_norm==".tsv/.csv file"){
              if (!is.null(input$uploadMeta_microarray_norm_tsv)){
                rv$metaData <- getMetaData(path = input$uploadMeta_microarray_norm_tsv$datapath,
                                           celfiles = colnames(exprs(rv$gxData)),
                                           filetype = input$MetaFileType_microarray_norm)
                
              }else{
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Error!!",
                  text = "You forgot to upload an expression and/or metadata file",
                  type = "error")
                shinybusy::remove_modal_spinner()
              }
            }
            if (input$MetaFileType_microarray_norm=="Series Matrix File"){
              if (!is.null(input$uploadMeta_microarray_norm_smf)){
                rv$metaData <- getMetaData(path = input$uploadMeta_microarray_norm_smf$datapath,
                                           celfiles =  colnames(exprs(rv$gxData)),
                                           filetype = input$MetaFileType_microarray_norm)
                
              }else{
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Error!!",
                  text = "You forgot to upload an expression and/or metadata file",
                  type = "error")
                shinybusy::remove_modal_spinner()
              }
            }
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(Biobase::exprs(rv$gxData))){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- Biobase::ExpressionSet(assayData = Biobase::exprs(rv$gxData)[,rownames(rv$metaData)])
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_microarray_norm <- DT::renderDataTable(NULL)
              output$metaTable_microarray_norm <- DT::renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_microarray_norm <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(Biobase::exprs(rv$gxData),6)
                colnames(output) <- stringr::str_remove(colnames(output),"\\.CEL.*")
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              # output$metaTable_microarray_norm <- DT::renderDataTable({
              #   req(rv$metaData)
              #   return(rv$metaData)
              # }, options = list(pageLength = 6))
              
              output$metaTable_microarray_norm <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              
              observeEvent(input$metaTable_microarray_norm_cell_edit, {
                row  <- input$metaTable_microarray_norm_cell_edit$row
                clmn <- input$metaTable_microarray_norm_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_microarray_norm_cell_edit$value
              })
              
              
              # Render UI for main tab
              output$UI_upload_microarray_norm <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_microarray_norm") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_microarray_norm") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_microarray_norm <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_microarray_norm")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(exprs(rv$gxData))){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_microarray_norm",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # Observe event
          
          #----------------------------------------------------------------------#
          # Run Example
          #----------------------------------------------------------------------#
          
          observeEvent(input$example_microarray_norm,{
            
            # Show modal
            shinybusy::show_modal_spinner(text = "Reading data...",
                                          color="#0dc5c1")
            
            # Read expression data
            rv$gxData <- getGEO(filename="Data/Microarray/GSE6955_series_matrix.txt.gz")
            
            # Guess the organism
            rv$Organism <- getOrganism(gxData = rv$gxData)
            
            # Get metadata
            rv$metaData <- getMetaData(path = "Data/Microarray/metaData_GSE6955.csv",
                                       celfiles = colnames(Biobase::exprs(rv$gxData)),
                                       filetype = ".tsv/.csv file")
            
            
            # Read raw expression data
            req(rv$metaData)
            if(nrow(rv$metaData)>0){
              
              # check if some samples are removed
              if (nrow(rv$metaData) != ncol(exprs(rv$gxData))){
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Warning!",
                  text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                  type = "warning")
              }
              
              # Filter expression data for samples with metadata
              rv$gxData <- Biobase::ExpressionSet(assayData = Biobase::exprs(rv$gxData)[,rownames(rv$metaData)])
              
              
              #------------------------------------------------------------------#
              # Outputs
              #------------------------------------------------------------------#
              
              # Set tables to zero
              output$exprTable_upload_microarray_norm <- DT::renderDataTable(NULL)
              output$metaTable_microarray_norm <- DT::renderDataTable(NULL)
              
              # Print expression table
              output$exprTable_upload_microarray_norm <- DT::renderDataTable({
                req(rv$gxData)
                output <- head(Biobase::exprs(rv$gxData),6)
                colnames(output) <- stringr::str_remove(colnames(output),"\\.CEL.*")
                return(output)
                
              }, options = list(pageLength = 6))
              
              # Print meta table
              output$metaTable_microarray_norm <- DT::renderDT({
                DT::datatable(rv$metaData, editable = TRUE)
              })
              
              
              observeEvent(input$metaTable_microarray_norm_cell_edit, {
                row  <- input$metaTable_microarray_norm_cell_edit$row
                clmn <- input$metaTable_microarray_norm_cell_edit$col
                rv$metaData[row, clmn] <- input$metaTable_microarray_norm_cell_edit$value
              })
              
              
              # Render UI for main tab
              output$UI_upload_microarray_norm <- renderUI({
                tagList(
                  tabsetPanel(
                    tabPanel("Expression matrix",
                             h3(strong("Expression matrix")),
                             h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "exprTable_upload_microarray_norm") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                    tabPanel("Meta data",                  # Meta table
                             h3(strong("Meta data")),
                             h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                             hr(),
                             DT::dataTableOutput(outputId = "metaTable_microarray_norm") %>% 
                               shinycssloaders::withSpinner(color="#0dc5c1")),
                  )
                )
              })
              
              # Allow user to go to next tab
              output$next_upload_microarray_norm <- renderUI({
                req(rv$metaData)
                req(rv$gxData)
                
                # Remove modal
                shinybusy::remove_modal_spinner()
                
                # Show RNA-seq upload tab
                showTab("navbar", target = "panel_preprocessing_microarray_norm")
                
                # Show message
                if (nrow(rv$metaData) >= ncol(Biobase::exprs(rv$gxData))){
                  shinyWidgets::sendSweetAlert(
                    session = session,
                    title = "Info",
                    text = "The data has been uploaded. Please check the tables on this 
                page to see whether the data has been correctly uploaded.",
                    type = "info")
                }
                
                # Show "next" button
                tagList(
                  hr(),
                  h2(strong("Continue your analysis")),
                  shinyWidgets::actionBttn(inputId = "next_upload_microarray_norm",
                                           label = "Next",
                                           style = "jelly",
                                           color = "danger",
                                           icon = icon("arrow-right"))
                )
              })
              
            } else {
              # No common samples between meta data and expression data
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Error!!",
                text = "No common samples",
                type = "error")
              
              shinybusy::remove_modal_spinner()
            }
            
          }) # Observe event
          
          
          #======================================================================#
          # Data Pre-processing
          #======================================================================#
          
          # Go to data upload tab
          observeEvent(input$next_upload_microarray_norm,{
            
            # Go to RNA-seq tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_preprocessing_microarray_norm")
          })
          
          
          # 1. Select outliers
          output$UI_outlier_microarray_norm <- renderUI({
            if(!input$outlier_norm){
              samples <- rownames(rv$metaData)
              
              pickerInput(inputId = "select_outliers_microarray_norm",
                          label = tags$span(
                            "Select samples to be removed", 
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select one or more samples to exclude from the analysis.", 
                                         position = "right",
                                         size = "large")
                          ),
                          choices = samples,
                          multiple = TRUE)
            } else{
              NULL
            }
            
          })
          
          # 2. Select experimental groups
          output$UI_groupselect_microarray_norm <- renderUI({
            pickerInput(inputId = "groupselect_microarray_norm",
                        label = NULL,
                        choices = colnames(rv$metaData),
                        selected = autoGroup(rv$metaData),
                        multiple = TRUE)
          })
          
          # print the experimental levels
          output$experimentallevels_microarray_norm <- renderText({
            req(input$groupselect_microarray_norm)
            
            if(length(input$groupselect_microarray_norm) > 1){
              experimentFactor <- make.names(apply(rv$metaData[,input$groupselect_microarray_norm], 1, paste, collapse = "_" ))
              } else{
              experimentFactor <- make.names(rv$metaData[,input$groupselect_microarray_norm])
            }
            levels <- unique(experimentFactor)
            if (length(levels) <= 6){
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(levels, collapse = ", "),
                             "</p>")
            } else {
              text <- paste0("<p><b>Levels (",length(levels),"):</b> ", paste(head(levels), collapse = ", "),
                             "...</p>")
            }
            return(text)
            
          })
          
          # Automatically select transformation
          output$UI_transformation_microarray_norm <- renderUI({
            
            if (checkTransformation(gxMatrix = exprs(rv$gxData)) == "Count matrix"){
              tagList(
                prettyRadioButtons(inputId = "transformation_microarray_norm", 
                                   label = NULL, 
                                   choices = c("Log2-transformation",
                                               "None"),
                                   inline = TRUE, 
                                   status = "danger",
                                   fill = TRUE,
                                   selected = "Log2-transformation"),
                h5(strong("NOTE: "),"Data seems not to be log-transformed. So, 
                 log2-transformation is needed.")
                
              )
              
            }else{
              tagList(
                prettyRadioButtons(inputId = "transformation_microarray_norm", 
                                   label = NULL, 
                                   choices = c("Log2-transformation",
                                               "None"),
                                   inline = TRUE, 
                                   status = "danger",
                                   fill = TRUE,
                                   selected = "None"),
                h5(strong("NOTE: "),"Data seems to be already log-transformed. So, 
                   no additional transformation is needed.")
              )
            }
          })
          
          
          # 3. pre-process the raw data
          observeEvent(input$start_preprocessing_microarray_norm,{
            
            # Show modal
            show_modal_spinner(text = "Pre-processing data...",
                               color="#0dc5c1")
            
            # Select outlier
            if (!isTRUE(input$outier)){
              rv$outlier <- input$select_outliers_microarray_norm
            } else{
              rv$outlier <- NULL
            }
            
            # Remove outlier
            if (!is.null(rv$outlier)){
              gxMatrix_temp <- exprs(rv$gxData)
              rv$gxData_fil <- ExpressionSet(assayData = gxMatrix_temp[,setdiff(colnames(gxMatrix_temp),rv$outlier)])
              rm(gxMatrix_temp)
            } else
              rv$gxData_fil <- rv$gxData
            
            # Filter metadata and expression data (samples in correct order)
            rv$metaData_fil <- rv$metaData[stringr::str_remove(colnames(exprs(rv$gxData_fil)),"\\.CEL.*"),]
            
            # Experiment factor
            if(length(input$groupselect_microarray_norm) > 1){
              rv$experimentFactor <- factor(make.names(apply(rv$metaData_fil[,input$groupselect_microarray_norm], 1, paste, collapse = "_" )))
              rv$experimentName <- input$groupselect_microarray_norm
            } else{
              rv$experimentFactor <- factor(make.names(rv$metaData_fil[,input$groupselect_microarray_norm]))
              rv$experimentName <- input$groupselect_microarray_norm
              
            }
            
            # Normalization
            rv$normData <- microarrayNormalization_processed(gxData = rv$gxData_fil,
                                                             experimentFactor = rv$experimentFactor,
                                                             transMeth = input$transformation_microarray_norm,
                                                             normMeth = input$normMeth_microarray_norm,
                                                             perGroup_name = input$perGroup_microarray_norm)
            
            rv$normMatrix <- exprs(rv$normData)
            
            # Collect chosen pre-processing settings into a dataframe
            rv$processingSettings <- data.frame(
              Option = c("Removed samples",
                         "Experimental group",
                         "Experimental levels",
                         "Normalization method",
                         "Transformation"),
              Selected = c(paste(rv$outlier, collapse = "; "),
                           paste(input$groupselect_microarray_norm, collapse = "; "),
                           paste(unique(rv$experimentFactor), collapse = "; "),
                           paste(input$normMeth_microarray_norm,"; ",input$perGroup_microarray_norm),
                           input$transformation_microarray_norm
              )
            )
            
            #======================================================================#
            # QC
            #======================================================================#
            
            #********************************************************************#
            # Output 1: Expression values
            #********************************************************************#
            
            # Print expression table
            output$exprTable_microarray_norm <- DT::renderDataTable({
              
              # Remove modal
              shinybusy::remove_modal_spinner()
              
              # Show message
              shinyWidgets::sendSweetAlert(
                session = session,
                title = "Info",
                text = "The data has been pre-processed. Please check the different 
              QC plots on this page to assess the pre-processing quality.",
                type = "info")
              
              # Show microarray statistics tab
              showTab("navbar", target = "panel_statistics_microarray_norm")
              
              output <- rv$normMatrix
              return(output)
              
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$downloadNormalizedData_microarray_norm <- downloadHandler(
              filename = "normalizedData.csv",
              content = function(file){
                write.csv(rv$normMatrix, file, quote = FALSE, row.names = TRUE)
              }
            )
            
            # Change color by click on button
            rv$colorOrder <- 1:length(levels(rv$experimentFactor))
            observeEvent(input$geneboxplot_changeOrder_microarray_norm,{
              all_orders <- permute(1:length(levels(rv$experimentFactor))) 
              sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
              
              if (sel < length(all_orders)){
                rv$colorOrder <- all_orders[[sel+1]]
              } else{
                rv$colorOrder <- all_orders[[1]]
              }
            })
            
            
            # Boxplot of single gene (based on selected row in the expression matrix)
            output$ExprBoxplot_microarray_norm <- renderPlot({
              req(input$exprTable_microarray_norm_rows_selected)
              
              if (length(levels(rv$experimentFactor)) > 4){
                legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
              } else{
                legendColors <- c(input$geneboxplot_col1_microarray_norm,
                                  input$geneboxplot_col2_microarray_norm,
                                  input$geneboxplot_col3_microarray_norm,
                                  input$geneboxplot_col4_microarray_norm,
                                  input$geneboxplot_col5_microarray_norm,
                                  input$geneboxplot_col6_microarray_norm)
              }
              # Make boxplot
              rv$temp1 <- geneBoxplot(experimentFactor = rv$experimentFactor, 
                                     normMatrix = rv$normMatrix, 
                                     sel_row = input$exprTable_microarray_norm_rows_selected,
                                     legendColors = legendColors[rv$colorOrder],
                                     groupOrder = input$geneboxplot_order_microarray_norm,
                                     rnaseq = FALSE,
                                     jitter = input$jitter_geneboxplot_microarray_norm,
                                     seed = sample(1:1000,1))
              return(rv$temp1)
            })
            
            
            # Get number of experimental groups
            output$length_geneboxplot_microarray_norm <- reactive({
              length(levels(rv$experimentFactor))
            })
            outputOptions(output, "length_geneboxplot_microarray_norm", suspendWhenHidden = FALSE) 
            
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_geneboxplot_microarray_norm <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp1, 
                                filename = file,
                                width = input$width_geneboxplot_microarray_norm,
                                height = input$height_geneboxplot_microarray_norm,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_geneboxplot_microarray_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_geneboxplot_microarray_norm", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_geneboxplot_microarray_norm", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_geneboxplot_microarray_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Output 2: Boxplots
            #********************************************************************#
            
            # Boxplots of all genes together
            output$boxplots_microarray_norm <- renderImage({
              req(session$clientData$output_boxplots_microarray_norm_width)
              req(session$clientData$output_boxplots_microarray_norm_height)
              getBoxplots(experimentFactor = rv$experimentFactor,
                          normData = rv$normData,
                          RNASeq = FALSE,
                          width = session$clientData$output_boxplots_microarray_norm_width,
                          height = session$clientData$output_boxplots_microarray_norm_height)
            }, deleteFile = TRUE)
            
            #***************************#
            # Modal to download figure
            #***************************#
            
            # Download plot
            output$realdownload_boxplots_microarray_norm <- downloadHandler(
              filename = function(){"QC_Boxplots.png"},
              content = function(file){
                png(file,
                    width=input$width_boxplots_microarray_norm,
                    height=input$height_boxplots_microarray_norm,
                    pointsize=24)
                getBoxplots_download(experimentFactor = rv$experimentFactor,
                                     normData = rv$normData,
                                     RNASeq = FALSE)
                dev.off()
              }
            )
            
            
            # Make modal
            observeEvent(input$download_boxplots_microarray_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                size = "m",
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_boxplots_microarray_norm", 
                                       "Height",
                                       min = 1200, max = 1600,
                                       value = 1500, step = 1,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_boxplots_microarray_norm", 
                                       "Width",
                                       min = 800, max = 1200,
                                       value = 1000, step = 1,
                                       width = "100%"),
                    )
                  ),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_boxplots_microarray_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            
            #********************************************************************#
            # Output 3: Densityplots
            #********************************************************************#
            
            # Densityplots of all genes together
            output$densityplots_microarray_norm <- renderPlotly({
              getDensityplots(experimentFactor = rv$experimentFactor,
                              normMatrix = rv$normMatrix)
              
            })
            
            #********************************************************************#
            # Output 4: Heatmap
            #********************************************************************#
            
            # Heatmap of sample-sample correlations
            output$heatmap_microarray_norm  <- renderPlotly({
              
              # Make colors
              if(length(input$heatmapFactor_microarray_norm) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$heatmapFactor_microarray_norm], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$heatmapFactor_microarray_norm])
              }
              
              getHeatmap(experimentFactor = colorFactor, 
                         normMatrix = rv$normMatrix,
                         clusterOption1 = input$clusteroption1_microarray_norm,
                         clusterOption2 = input$clusteroption2_microarray_norm,
                         theme = input$heatmaptheme_microarray_norm)
            })
            
            
            #********************************************************************#
            # Output 5: PCA
            #********************************************************************#
            
            #Perform PCA
            rv$PCA_data <- prcomp(t(rv$normMatrix[apply(rv$normMatrix, 1, var) != 0,]),
                                  retx = TRUE, 
                                  center = TRUE,
                                  scale.= TRUE)
            
            
            # Make PCA plot
            output$PCA_microarray_norm <- renderPlotly({
              
              if(length(input$colorFactor_microarray_norm) > 1){
                colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_microarray_norm], 1, paste, collapse = "_" ))
              } else{
                colorFactor <- factor(rv$metaData_fil[,input$colorFactor_microarray_norm])
              }
              
              plot_PCA(PC_data = rv$PCA_data, 
                       colorFactor = colorFactor, 
                       xpc = as.numeric(str_remove(input$xpca_microarray_norm,"PC")), 
                       ypc = as.numeric(str_remove(input$ypca_microarray_norm,"PC")), 
                       zpc = ifelse(input$xyz_microarray_norm,as.numeric(str_remove(input$zpca_microarray_norm,"PC")),3), 
                       xyz = input$xyz_microarray_norm)
              
            })
            
            #********************************************************************#
            # Output 6: Overview of pre-processing settings
            #********************************************************************#
            
            # Print table with settings
            output$processingSettings_microarray_norm <- DT::renderDataTable({
              return(rv$processingSettings)
            },options = list(pageLength = 10),
            selection = "none")
            
            # Download button
            output$downloadProcessingSettings_microarray_norm <- downloadHandler(
              filename = "preprocessingSettings.csv",
              content = function(file){
                write.csv(rv$processingSettings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
            #********************************************************************#
            # UI
            #********************************************************************#
            
            # Render UI for main tab
            output$UI_QC_microarray_norm <- renderUI({
              tagList(
                tabsetPanel(
                  
                  # TAB1: expression values
                  tabPanel("Expression values",
                           icon = icon("fas fa-mouse-pointer"),
                           h3(strong("Normalized expression values")),
                           h5("Here you can view the normalized log intensity (expression) values."),
                           hr(),
                           DT::dataTableOutput(outputId = "exprTable_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadNormalizedData_microarray_norm", 
                                          "Download table"),
                           br(),
                           br(),
                           
                           # Dropdown Button to adjust the plot settings
                           shinyWidgets::dropdownButton(
                             tags$div(
                               style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                               
                               
                               # Change order of the boxplots:
                               tags$h4("Drag to change boxplot order"),
                               shinyjqui::orderInput(inputId = 'geneboxplot_order_microarray_norm', 
                                                     label = NULL, 
                                                     items = levels(rv$experimentFactor),
                                                     item_class = 'default'),
                               br(),
                               
                               # Change colour of the boxplots by button. 
                               # This is used when there are more than 6 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_microarray_norm > 6",
                                 tags$h4("Click to change boxplot colours"),
                                 shinyWidgets::actionBttn("geneboxplot_changeOrder_microarray_norm",
                                                          label = "Change color",
                                                          style = "simple",
                                                          color = "primary",
                                                          icon = icon("sync"))
                               ),
                               
                               # Change colour of the boxplots by colour picker
                               # This is used when there are less than 7 experimental groups
                               conditionalPanel(
                                 condition = "output.length_geneboxplot_microarray_norm < 7",
                                 tags$h4("Click to select boxplot colours"),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 0",
                                   colourpicker::colourInput("geneboxplot_col1_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[1])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 1",
                                   colourpicker::colourInput("geneboxplot_col2_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[2])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 2",
                                   colourpicker::colourInput("geneboxplot_col3_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[3])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 3",
                                   colourpicker::colourInput("geneboxplot_col4_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[4])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 4",
                                   colourpicker::colourInput("geneboxplot_col5_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[5])
                                 ),
                                 conditionalPanel(
                                   condition = "output.length_geneboxplot_microarray_norm > 5",
                                   colourpicker::colourInput("geneboxplot_col6_microarray_norm", 
                                                             NULL, 
                                                             colorsByFactor(rv$experimentFactor)$legendColors[6])
                                 )
                               ),
                               br(),
                               tags$h4("Drag to change jitter"),
                               sliderInput("jitter_geneboxplot_microarray_norm", 
                                           NULL,
                                           min = 0, max = 0.3,
                                           value = 0.1, step = 0.01),
                               br()
                             ),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"),
                             
                             tooltip = shinyWidgets::tooltipOptions(
                               title = "Click to personalize plot!")
                             
                           ), # EO dropdownButton
                           
                           # Boxplot of the selected gene's expression values
                           plotOutput("ExprBoxplot_microarray_norm")%>% 
                             shinycssloaders::withSpinner(color="#0dc5c1"),
                           
                           actionButton("download_geneboxplot_microarray_norm", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           br(),
                           br()
                  ),
                  
                  # TAB2: boxplots
                  tabPanel("Boxplots",
                           icon = icon("fas fa-file"),
                           br(),
                           actionButton("download_boxplots_microarray_norm", 
                                        "Download figure",
                                        icon = shiny::icon("download")),
                           plotOutput(outputId = "boxplots_microarray_norm",
                                      width = "65vw", height = "80vw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # TAB3: Density plots
                  tabPanel("Density plots",
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           h2(strong("Density plot of normalized intensities"), align = "center"),
                           h4("Distributions should be comparable between arrays", align = "center"),
                           plotly::plotlyOutput(outputId = "densityplots_microarray_norm",
                                                width = "65vw", height = "40vw") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1")
                  ),
                  
                  # TAB4: sample-sample correlations
                  tabPanel("Correlation plot", 
                           icon = icon("fas fa-mouse-pointer"),
                           
                           br(),
                           
                           fluidRow(
                             
                             #Distance method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption1_microarray_norm",
                                      label = "Distance calculation 
                                              method",
                                      choices = c("Pearson","Spearman",
                                                  "Euclidean"),
                                      options = list(
                                        style = "btn-primary"))
                             ),
                             
                             #Ckustering method
                             column(4,
                                    shinyWidgets::pickerInput(
                                      inputId = "clusteroption2_microarray_norm",
                                      label = "Clustering method",
                                      choices = c("ward.D2","single",
                                                  "complete","average",
                                                  "mcquitty","median",
                                                  "centroid"),
                                      options = list(
                                        style = "btn-info"))
                             )),
                           
                           hr(),
                           
                           #Theme
                           dropdownButton(
                             tags$h3("Theme"),
                             
                             selectInput(inputId = "heatmapFactor_microarray_norm",
                                         label = "Side colors",
                                         choices = colnames(rv$metaData_fil),
                                         selected = rv$experimentName,
                                         multiple = TRUE),
                             
                             selectInput(inputId = 'heatmaptheme_microarray_norm',
                                         label = "Heatmap theme",
                                         choices = c("Default", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             circle = TRUE, status = "info",
                             icon = icon("fas fa-cog"), width = "300px",
                             tooltip = tooltipOptions(
                               title = "Click to change colors!")
                             
                           ),
                           
                           plotly::plotlyOutput("heatmap_microarray_norm", 
                                                width = "1000px", 
                                                height="600px") %>% 
                             shinycssloaders::withSpinner(color="#0dc5c1", 
                                                          proxy.height = "400px")
                           
                  ),
                  
                  # TAB5: PCA
                  tabPanel("PCA",
                           icon = icon("fas fa-mouse-pointer"),
                           # Title + text
                           fluidRow(
                             h3(strong("Principal Component Analysis (PCA)")),
                             h5("Here you can view the PCA score plot."),
                             hr(),
                           ),
                           
                           # Set color + 3D/2D
                           fluidRow(
                             column(3,
                                    # Color by which factor?
                                    shinyWidgets::pickerInput(inputId = "colorFactor_microarray_norm",
                                                              label = "Color by:",
                                                              choices = colnames(rv$metaData_fil),
                                                              selected = rv$experimentName,
                                                              multiple = TRUE)
                             ),
                             column(3,
                                    br(),
                                    # 3D plot?
                                    shinyWidgets::materialSwitch(
                                      inputId = "xyz_microarray_norm",
                                      label = "3D",
                                      value = FALSE, 
                                      status = "danger")
                                    
                             )
                           ),
                           
                           # Set axes
                           fluidRow(
                             column(3,
                                    #X-axis
                                    selectInput(inputId = "xpca_microarray_norm", 
                                                label = "x-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC1")
                             ),
                             column(3,
                                    #Y-axis
                                    selectInput(inputId = "ypca_microarray_norm", 
                                                label = "y-axis",
                                                choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                            "PC6", "PC7", "PC8"),
                                                selected = "PC2")
                             ),
                             column(3,
                                    #Z-axis
                                    conditionalPanel(
                                      condition = "input.xyz_microarray_norm==true",
                                      selectInput(inputId = "zpca_microarray_norm", 
                                                  label = "z-axis",
                                                  choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                              "PC6", "PC7", "PC8"),
                                                  selected = "PC3")
                                    )
                             )
                           ),
                           
                           # Print plot
                           fluidRow(
                             hr(),
                             plotly::plotlyOutput("PCA_microarray_norm")%>% 
                               withSpinner(color="#0dc5c1")
                           ) 
                           
                  ), # PCA tabpanel
                  
                  # # TAB6: settings table
                  tabPanel("Settings overview",
                           icon = icon("fas fa-file"),
                           h3(strong("Pre-processing settings")),
                           h5("Here you can see an overview of the chosen pre-processing settings."),
                           hr(),
                           DT::dataTableOutput(outputId = "processingSettings_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("downloadProcessingSettings_microarray_norm", 
                                          "Download")
                           
                  ) # EO Settings tabPanel
                ) # EO tabsetPanel
              ) # EO tagList
            }) # EO renderUI
            
            
            # Allow user to go to next tab
            output$UI_next_preprocessing_microarray_norm <- renderUI({
              req(rv$normData)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                actionBttn(inputId = "next_preprocessing_microarray_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          
          #======================================================================#
          # Statistical analysis
          #======================================================================#
          # Go to data upload tab
          observeEvent(input$next_preprocessing_microarray_norm,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_statistics_microarray_norm")
          })
          
          # SIDEPANEL:
          
          
          # Select continuous covariates
          output$UI_covGroups_num_microarray_norm <- renderUI({
            tagList(
              pickerInput(inputId = "covGroups_num_microarray_norm",
                          label = "Continuous covariates (e.g., age):",
                          choices = setdiff(colnames(rv$metaData_fil),
                                            rv$experimentName),
                          selected = NULL,
                          multiple = TRUE)
            )
          })
          
          # Select discrete covariates
          output$UI_covGroups_char_microarray_norm <- renderUI({
            tagList(
              pickerInput(inputId = "covGroups_char_microarray_norm",
                          label = "Discrete covariates (e.g., sex):",
                          choices = setdiff(colnames(rv$metaData_fil),
                                            rv$experimentName),
                          selected = NULL,
                          multiple = TRUE)
            )
          })
          
          # Select comparisons
          output$UI_comparisons_microarray_norm <- renderUI({
            
            tagList(
              multiInput(
                inputId = "comparisons_microarray_norm",
                label = "Comparisons:", 
                choices = makeComparisons(make.names(unique(rv$experimentFactor))),
                selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1]
              )
            )
          })
          
          
          # Select comparisons
          output$UI_biomart_dataset_microarray_norm <- renderUI({
            req(input$addAnnotation_microarray_norm)
            pickerInput(inputId = "biomart_dataset_microarray_norm",
                        label = "Organism:",
                        choices = c("Homo sapiens" = "hsapiens_gene_ensembl" ,
                                    "Bos taurus" = "btaurus_gene_ensembl",
                                    "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                    "Mus musculus" = "mmusculus_gene_ensembl",
                                    "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                        selected = switch(rv$Organism,
                                          "Homo sapiens" = "hsapiens_gene_ensembl" ,
                                          "Bos taurus" = "btaurus_gene_ensembl",
                                          "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                          "Mus musculus" = "mmusculus_gene_ensembl",
                                          "Rattus norvegicus" = "rnorvegicus_gene_ensembl"),
                        multiple = FALSE)
          })
          
          output$UI_addAnnotations_microarray_norm <- renderUI({
            req(input$addAnnotation_microarray_norm)
            req(input$biomart_dataset_microarray_norm)
            
            tagList(
              
              pickerInput(inputId = "biomart_filter_microarray_norm",
                          label = "Probeset ID",
                          choices = filterList[[input$biomart_dataset_microarray_norm]],
                          selected = selFilter(rv$ProbeAnnotation),
                          multiple = FALSE),
              
              pickerInput(inputId = "biomart_attributes_microarray_norm",
                          label = "Output",
                          choices = c("Ensembl Gene ID",
                                      "Entrez Gene ID",
                                      "Gene Symbol/Name"),
                          selected = "Gene Symbol/Name",
                          multiple = TRUE)
            )
            
          })
          
          
          
          #======================================================================#
          # Output of statistical analysis
          #======================================================================#
          observeEvent(input$calculate_statistics_microarray_norm,{
            show_modal_spinner(text = "Statistical analysis...",
                               color="#0dc5c1")
            
            # Calculate statistics
            if (isTRUE(input$addAnnotation_microarray_norm)){
              rv$top_table_list <- getStatistics(normMatrix = rv$normMatrix, 
                                            metaData = rv$metaData_fil, 
                                            expFactor = rv$experimentName, 
                                            covGroups_num = input$covGroups_num_microarray_norm,
                                            covGroups_char = input$covGroups_char_microarray_norm,
                                            comparisons = input$comparisons_microarray_norm,
                                            addAnnotation = input$addAnnotation_microarray_norm,
                                            biomart_dataset = input$biomart_dataset_microarray_norm,
                                            biomart_attributes = unique(c(input$biomart_filter_microarray_norm,
                                                                          input$biomart_attributes_microarray_norm)),
                                            biomart_filters = input$biomart_filter_microarray_norm)
              
            } else{
              rv$top_table_list <- getStatistics(normMatrix = rv$normMatrix, 
                                            metaData = rv$metaData_fil, 
                                            expFactor = rv$experimentName, 
                                            covGroups_num = input$covGroups_num_microarray_norm,
                                            covGroups_char = input$covGroups_char_microarray_norm,
                                            comparisons = input$comparisons_microarray_norm,
                                            addAnnotation = input$addAnnotation_microarray_norm,
                                            biomart_dataset = NULL,
                                            biomart_attributes = NULL,
                                            biomart_filters = NULL)
              
            }
            
            # Select comparison for output
            observe({
              rv$top_table <- rv$top_table_list[[1]]
              
              if (!is.null(rv$top_table_list)){
                # Remove modal
                remove_modal_spinner()
                
                # Show comparisons
                output$UI_comparisons_view_microarray_norm <- renderUI({
                  pickerInput(inputId = "comparisons_view_microarray_norm",
                              label = "Select comparison:",
                              choices = names(rv$top_table),
                              selected = names(rv$top_table)[1],
                              multiple = FALSE)
                })
                
                # Show message
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = rv$top_table_list[[2]],
                  type = "info")
                
                # Show microarray statistics tab
                showTab("navbar", target = "panel_ORA_microarray_norm")
              } else{
                remove_modal_spinner()
                # Show message
                sendSweetAlert(
                  session = session,
                  title = "Error",
                  text = "Oops...something went wrong! Please try again!",
                  type = "error")
              }
              
            })
            
            
            # Make outputs:
            
            #********************************************************************#
            # top table
            #********************************************************************#
            
            # print top table
            output$top_table_microarray_norm <- renderDataTable({
              req(input$comparisons_view_microarray_norm)
              
              output <- rv$top_table[[input$comparisons_view_microarray_norm]]
              
              return(output)
            },options = list(pageLength = 6),
            selection = list(mode = "single", selected = 1), escape = FALSE)
            
            # Download button
            output$download_top_table_microarray_norm <- downloadHandler(
              filename = paste0("topTable_",input$comparisons_view_microarray_norm,".csv"),
              content = function(file){
                write.csv(rv$top_table[[input$comparisons_view_microarray_norm]], file, quote = FALSE, row.names = FALSE)
              }
            )
            
            # Change plotting data depending on whether all experimental groups will be plotted  
            observe({
              req(input$comparisons_view_microarray_norm)
              req(rv$top_table)
              
              if (!is.null(input$boxplotAll_microarray_norm)){
                if (input$boxplotAll_microarray_norm){
                  rv$newFactor <- rv$experimentFactor
                  rv$newData <- rv$normMatrix
                }
                if (!input$boxplotAll_microarray_norm){
                  if(length(rv$experimentName) > 1){
                    t <- make.names(apply(rv$metaData_fil[,rv$experimentName], 1, paste, collapse = "_" ))
                  } else{
                    t <- make.names(rv$metaData_fil[,rv$experimentName])
                  }
                  rv$newData <- rv$normMatrix[,t %in% (unlist(stringr::str_split(input$comparisons_view_microarray_norm, " - ")))]
                  rv$newFactor <- factor(as.character(rv$experimentFactor)[t %in% (unlist(stringr::str_split(input$comparisons_view_microarray_norm, " - ")))])
                }
              }
              
              # Change color by click on button
              rv$colorOrder <- 1:length(levels(rv$newFactor))
              observeEvent(input$statboxplot_changeOrder_microarray_norm,{
                all_orders <- permute(1:length(levels(rv$newFactor))) 
                sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder))))
                
                if (sel < length(all_orders)){
                  rv$colorOrder <- all_orders[[sel+1]]
                } else{
                  rv$colorOrder <- all_orders[[1]]
                }
              })
              
              
              # Boxplot of single gene (based on selected row in the top table)
              output$ExprBoxplot_statistics_microarray_norm <- renderPlot({
                req(input$top_table_microarray_norm_rows_selected)
                req(rv$top_table)
                
                if (length(levels(rv$newFactor)) > 6){
                  legendColors <- colorsByFactor(rv$newFactor)$legendColors
                } else{
                  legendColors <- c(input$statboxplot_col1_microarray_norm,
                                    input$statboxplot_col2_microarray_norm,
                                    input$statboxplot_col3_microarray_norm,
                                    input$statboxplot_col4_microarray_norm,
                                    input$statboxplot_col5_microarray_norm,
                                    input$statboxplot_col6_microarray_norm)
                }
                
                gene <- rv$top_table[[input$comparisons_view_microarray_norm]]$GeneID[input$top_table_microarray_norm_rows_selected]
                sel_row <- which(as.character(rownames(rv$normMatrix)) %in% as.character(gene))
                
                # Make boxplot
                rv$temp <- geneBoxplot(experimentFactor = rv$newFactor, 
                                       normMatrix = rv$newData, 
                                       sel_row = sel_row,
                                       legendColors = legendColors[rv$colorOrder],
                                       groupOrder = input$statboxplot_order_microarray_norm,
                                       jitter = input$jitter_statboxplot_microarray_norm,
                                       rnaseq=FALSE,
                                       seed = sample(1:1000,1))
                
                return(rv$temp)
                
              })
              
              # Get number of experimental groups
              output$length_statboxplot_microarray_norm <- reactive({
                length(levels(rv$experimentFactor))
              })
              outputOptions(output, "length_statboxplot_microarray_norm", suspendWhenHidden = FALSE) 
              
            })
            
            #***************************#
            # Modal to download boxplot
            #***************************#
            
            # Download plot
            output$realdownload_statboxplot_microarray_norm <- downloadHandler(
              filename = "GeneBoxplot.png",
              content = function(file){
                ggplot2::ggsave(plot = rv$temp, 
                                filename = file,
                                width = input$width_statboxplot_microarray_norm,
                                height = input$height_statboxplot_microarray_norm,
                                units = "px")
              }
            )
            
            
            # Make modal
            observeEvent(input$download_statboxplot_microarray_norm, {
              showModal(modalDialog(
                title = NULL,
                easyClose = TRUE,
                footer = tagList(
                  fluidRow(
                    column(6,
                           sliderInput("height_statboxplot_microarray_norm", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    ),
                    column(6,
                           sliderInput("width_statboxplot_microarray_norm", 
                                       "Width",
                                       min = 800, max = 3000,
                                       value = 2100, step = 10,
                                       width = "100%"),
                    )
                  ),
                  hr(),
                  fluidRow(
                    column(12, align = "left",
                           downloadButton('realdownload_statboxplot_microarray_norm', 
                                          'Download')
                    )
                  )
                  
                )
                
              ))
            })
            
            #********************************************************************#
            # Histograms
            #********************************************************************#
            observe({
              req(input$comparisons_view_microarray_norm)
              
              if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
                output$Phistogram_microarray_norm <- renderPlotly({
                  req(rv$top_table)
                  p <- makePHistogram(rv$top_table[[input$comparisons_view_microarray_norm]][,"p-value"])
                  return(p)
                })
                
                output$logFChistogram_microarray_norm <- renderPlotly({
                  req(rv$top_table)
                  p <- makelogFCHistogram(rv$top_table[[input$comparisons_view_microarray_norm]][,"log2FC"])
                  return(p)
                })
              }
            })
            
            #********************************************************************#
            # Volcano plot
            #********************************************************************#
            observeEvent(input$plot_volcano_microarray_norm, {
              req(rv$top_table)
              req(input$rawp_volcano_microarray_norm)
              req(input$p_thres_volcano_microarray_norm)
              req(input$logFC_thres_volcano_microarray_norm)
              
              req(input$comparisons_view_microarray_norm)
              
              if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
                p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                                 p = input$rawp_volcano_microarray_norm, 
                                 p_threshold = input$p_thres_volcano_microarray_norm, 
                                 logFC_threshold = input$logFC_thres_volcano_microarray_norm)
                
                output$volcano_microarray_norm <- renderPlotly(p)
              }
            }, ignoreNULL = FALSE)
            
            
            #=========================================#
            #  # UI: Output in different tabs
            #=========================================#
            
            observe({
              req(rv$top_table)
              
              output$UI_boxplotAll_microarray_norm <- renderUI({
                tagList(
                  # Change order of the boxplots:
                  tags$h4("Drag to change boxplot order"),
                  shinyjqui::orderInput(inputId = 'statboxplot_order_microarray_norm', 
                                        label = NULL, 
                                        items = levels(rv$newFactor),
                                        item_class = 'default'),
                  br(),
                  
                  # Change colour of the boxplots by button. 
                  # This is used when there are more than 6 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_microarray_norm > 6",
                    tags$h4("Click to change boxplot colours"),
                    shinyWidgets::actionBttn("statboxplot_changeOrder_microarray_norm",
                                             label = "Change color",
                                             style = "simple",
                                             color = "primary",
                                             icon = icon("sync"))
                  ),
                  
                  # Change colour of the boxplots by colour picker
                  # This is used when there are less than 7 experimental groups
                  conditionalPanel(
                    condition = "output.length_statboxplot_microarray_norm < 7",
                    tags$h4("Click to select boxplot colours"),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 0",
                      colourpicker::colourInput("statboxplot_col1_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[1])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 1",
                      colourpicker::colourInput("statboxplot_col2_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[2])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 2",
                      colourpicker::colourInput("statboxplot_col3_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[3])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 3",
                      colourpicker::colourInput("statboxplot_col4_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[4])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 4",
                      colourpicker::colourInput("statboxplot_col5_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[5])
                    ),
                    conditionalPanel(
                      condition = "output.length_statboxplot_microarray_norm > 5",
                      colourpicker::colourInput("statboxplot_col6_microarray_norm", 
                                                NULL, 
                                                colorsByFactor(rv$newFactor)$legendColors[6])
                    )
                  ),
                  br(),
                  tags$h4("Drag to change jitter"),
                  sliderInput("jitter_statboxplot_microarray_norm", 
                              NULL,
                              min = 0, max = 0.3,
                              value = 0.1, step = 0.01),
                  br()
                )
              })
            })
            
            observe({
              if (is.null(rv$top_table)){
                output$UI_output_statistics_microarray_norm <- renderUI(NULL)
              } else{
                output$UI_output_statistics_microarray_norm <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #********************************************************************#
                      # top table tab
                      #********************************************************************#
                      
                      tabPanel("Top table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("Top Table")),
                               h5("The Top Table includes the output of the selected statistical analysis."),
                               hr(),
                               dataTableOutput(outputId = "top_table_microarray_norm") %>% 
                                 withSpinner(color="#0dc5c1"),
                               downloadButton("download_top_table_microarray_norm", 
                                              "Download table"),
                               br(),
                               br(),
                               
                               # Dropdown Button to adjust the plot settings
                               shinyWidgets::dropdownButton(
                                 tags$div(
                                   style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                                   
                                   # Plot all experimental groups?
                                   conditionalPanel(
                                     condition = "output.length_statboxplot_microarray_norm > 2",
                                     tags$h4("Plot all experimental groups?"),
                                     shinyWidgets::materialSwitch(inputId = "boxplotAll_microarray_norm",
                                                                  label = NULL, 
                                                                  value = TRUE,
                                                                  status = "primary"),
                                     
                                     br()
                                   ),
                                   uiOutput("UI_boxplotAll_microarray_norm")
                                   
                                   
                                   
                                 ),
                                 circle = TRUE, status = "info",
                                 icon = icon("fas fa-cog"),
                                 
                                 tooltip = shinyWidgets::tooltipOptions(
                                   title = "Click to personalize plot!")
                                 
                               ), # EO dropdownButton
                               plotOutput("ExprBoxplot_statistics_microarray_norm")%>% 
                                 shinycssloaders::withSpinner(color="#0dc5c1"),
                               actionButton("download_statboxplot_microarray_norm", 
                                            "Download figure",
                                            icon = shiny::icon("download")),
                               br(),
                               br()
                      ),
                      #********************************************************************#
                      # histogram tab
                      #********************************************************************#
                      
                      tabPanel("Histograms",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               plotlyOutput("Phistogram_microarray_norm")%>% 
                                 withSpinner(color="#0dc5c1"),
                               hr(),
                               plotlyOutput("logFChistogram_microarray_norm")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      ),
                      
                      #********************************************************************#
                      # volcano tab
                      #********************************************************************#
                      
                      tabPanel("Volcano plot",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               
                               fluidRow(
                                 column(3,
                                        # Raw or adjusted P value?
                                        prettyRadioButtons(
                                          inputId = "rawp_volcano_microarray_norm",
                                          label = "P value", 
                                          choices = 
                                            c("Raw P value" = "raw", 
                                              "Adjusted P value" = "adj"))
                                 ),
                                 column(3,
                                        #P value threshold
                                        numericInput(
                                          inputId = "p_thres_volcano_microarray_norm",
                                          label = "P threshold",
                                          value = 0.05)
                                 ),
                                 column(3,
                                        #logFC threshold
                                        numericInput(
                                          inputId = "logFC_thres_volcano_microarray_norm",
                                          label = "logFC threshold",
                                          value = 1)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_volcano_microarray_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               # Volcano plot
                               plotlyOutput("volcano_microarray_norm")%>% 
                                 withSpinner(color="#0dc5c1")
                               
                      )
                      
                    ) # tabsetpanel
                  ) # taglist
                }) # renderUI
              }
            }) # Observe
            
            # Allow user to go to next tab
            output$UI_next_statistics_microarray_norm <- renderUI({
              req(rv$top_table)
              tagList(
                hr(),
                h2(strong("Continue your analysis")),
                actionBttn(inputId = "next_statistics_microarray_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          }) # observeEvent
          
          
          #======================================================================#
          # ORA
          #======================================================================#
          # Go to data upload tab
          observeEvent(input$next_statistics_microarray_norm,{
            
            # Go to microarray statistics tab
            updateNavbarPage(session, "navbar",
                             selected = "panel_ORA_microarray_norm")
          })
          
          # Comparisons for which ORA should be performed
          observe({
            req(rv$top_table)
            output$UI_comparisons_view_ORA_microarray_norm <- renderUI({
              pickerInput(inputId = "comparisons_view_ORA_microarray_norm",
                          label = NULL,
                          choices = names(rv$top_table),
                          selected = names(rv$top_table)[1],
                          multiple = FALSE)
            })
          })
          
          # Which columns contain gene ids
          observe({
            req(input$comparisons_view_ORA_microarray_norm)
            col_choice <- 1
            if (length(colnames(rv$top_table[[input$comparisons_view_ORA_microarray_norm]])) > 6){
              col_choice <- c(1,7:ncol(rv$top_table[[input$comparisons_view_ORA_microarray_norm]]))
            }
            output$UI_geneID_ORA_microarray_norm <- renderUI({
              tagList(
                selectInput(inputId = "organism_ORA_microarray_norm",
                            label = "Organism",
                            choices = c("Bos taurus",
                                        "Caenorhabditis elegans",
                                        "Homo sapiens",
                                        "Mus musculus", 
                                        "Rattus norvegicus"),
                            selected = rv$Organism),
                
                pickerInput(inputId = "geneID_ORA_microarray_norm",
                            label = "Which column of the top table contains the gene IDs?",
                            choices = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_norm]])[col_choice],
                            selected = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_norm]])[1],
                            multiple = FALSE),
                
                pickerInput(inputId = "selID_ORA_microarray_norm",
                            label = "Which gene ID to use?",
                            choices = c("Ensembl Gene ID" = "ENSEMBL", 
                                        "Entrez Gene ID" = "ENTREZID", 
                                        "Gene Symbol/Name" = "SYMBOL"),
                            selected = selID(rv$ProbeAnnotation),
                            multiple = FALSE),
              )
            })
          })
          
          # Show modal
          observeEvent(input$calculate_ORA_microarray_norm,{
            show_modal_spinner(text = "Overrepresentation analysis...",
                               color="#0dc5c1")
            
            # Perform ORA
            if (input$topNorThres_microarray_norm == "Threshold"){
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                                 geneset = input$geneset_ORA_microarray_norm,
                                 geneID_col = input$geneID_ORA_microarray_norm,
                                 geneID_type = input$selID_ORA_microarray_norm,
                                 organism = input$organism_ORA_microarray_norm,
                                 updown = input$updown_ORA_microarray_norm,
                                 topN = FALSE,
                                 N = NULL,
                                 rawadj = input$rawp_ORA_microarray_norm,
                                 p_thres = input$p_thres_ORA_microarray_norm,
                                 logFC_thres = input$logFC_thres_ORA_microarray_norm)
            } else{
              rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                                 geneset = input$geneset_ORA_microarray_norm,
                                 geneID_col = input$geneID_ORA_microarray_norm,
                                 geneID_type = input$selID_ORA_microarray_norm,
                                 organism = input$organism_ORA_microarray_norm,
                                 updown = input$updown_ORA_microarray_norm,
                                 topN = TRUE,
                                 N = input$topN_microarray_norm,
                                 rawadj = NULL,
                                 p_thres = NULL,
                                 logFC_thres = NULL)
            }
            
            observe({
              
              # Remove modal
              remove_modal_spinner()
              
              # Show message
              if (is.null(rv$ORA_data)){
                sendSweetAlert(
                  session = session,
                  title = "Error!",
                  text = "No significant genes!",
                  type = "error")
                
              }else{
                
                sendSweetAlert(
                  session = session,
                  title = "Info",
                  text = "Gene overrepresentation analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
                  type = "info")
                
                # print GO table
                output$ORA_table_microarray_norm <- renderDataTable({
                  req(input$geneset_ORA_microarray_norm)
                  output <- rv$ORA_data@result
                  
                  if (input$geneset_ORA_microarray_norm == "WikiPathways"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.wikipathways.org/pathways/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  } 
                  
                  if (input$geneset_ORA_microarray_norm == "KEGG"){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.genome.jp/pathway/",
                        output$ID, ".html"
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  
                  if (input$geneset_ORA_microarray_norm %in% c("GO-BP", "GO-MF", "GO-CC")){
                    output$ID <- paste0(
                      '<a ',
                      'href=',
                      paste0(
                        "https://www.ebi.ac.uk/QuickGO/term/",
                        output$ID
                      ),
                      ' target="_blank"',
                      '>',
                      output$ID,
                      '</a>'
                    )
                  }
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  
                  return(output)
                },options = list(pageLength = 6),
                selection = list(mode = "single", selected = 1), escape = FALSE)
                
                # Download button
                output$download_ORA_table_microarray_norm <- downloadHandler(
                  filename = paste0("ORATable_",input$comparisons_view_ORA_microarray_norm,"_",input$geneset_ORA_microarray_norm,".csv"),
                  content = function(file){
                    write.csv(rv$ORA_data@result, file, quote = FALSE, row.names = FALSE)
                  }
                )
                
                # Print statistics of genes in selected Term
                output$ORAgene_table_microarray_norm <- renderDataTable({
                  req(input$ORA_table_microarray_norm_rows_selected)
                  
                  # Make ORA gene table
                  output <- make_ORAgene_table(ORA_data = rv$ORA_data,
                                               top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                                               geneID_col = input$geneID_ORA_microarray_norm,
                                               sel_row_ORA = input$ORA_table_microarray_norm_rows_selected)
                  
                  output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
                  output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
                  output$meanExpr <- round(output$meanExpr,3)
                  output$log2FC <- round(output$log2FC,3)
                  output$`log2FC SE` <- round(output$`log2FC SE`,3)
                  
                  return(output)
                }, options = list(pageLength = 6), escape = FALSE)
                
                # Text for gene table
                output$text_ORAgene_table_microarray_norm <- renderText({
                  text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_microarray_norm_rows_selected,"ID"],
                                 "</b></h3>")
                  return(text)
                })
                
                # ORA barchart
                observeEvent(input$plot_ORAplot_microarray_norm, {
                  req(input$nSets_ORAplot_microarray_norm)
                  p <- makeORAplot(rv$ORA_data,
                                   nSets = input$nSets_ORAplot_microarray_norm)
                  
                  output$ORAplot_microarray_norm <- renderPlotly(p)
                  
                }, ignoreNULL = FALSE)
                
                # ORA network diagram
                observeEvent(input$plot_ORAnetwork_microarray_norm, {
                  req(input$layout_network_microarray_norm)
                  req(input$nSets_network_microarray_norm)
                  p <- makeORAnetwork(ORA_data = rv$ORA_data,
                                      layout = input$layout_network_microarray_norm,
                                      nSets = input$nSets_network_microarray_norm)
                  
                  output$ORAnetwork_microarray_norm <- renderPlot(p,
                                                                  height = 500, 
                                                                  width = 800)
                  
                }, ignoreNULL = FALSE)
                
              }
            }) # Observe
            
            
            
            # Render output table
            observe({
              if (!is.null(rv$ORA_data)){
                output$UI_output_ORA_microarray_norm <- renderUI({
                  tagList(
                    tabsetPanel(
                      
                      #********************************************************************#
                      # top table tab
                      #********************************************************************#
                      
                      tabPanel("ORA table",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("Statistics table")),
                               h5("The ORA statistics table encompasses the output of the gene 
                              overrepresentation analysis."),
                               hr(),
                               dataTableOutput(outputId = "ORA_table_microarray_norm") %>% 
                                 withSpinner(color="#0dc5c1"),
                               downloadButton("download_ORA_table_microarray_norm", 
                                              "Download"),
                               br(),
                               htmlOutput("text_ORAgene_table_microarray_norm"),
                               h5(paste0("The gene table encompasses the statistics of all genes 
                              from the selected geneset.")),
                               hr(),
                               dataTableOutput(outputId = "ORAgene_table_microarray_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                      tabPanel("Bar chart",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("ORA bar chart")),
                               h5("The ORA bar chart visualizes the results from the overrepresentation analysis."),
                               hr(),
                               fluidRow(
                                 column(3,
                                        # Number of genesets
                                        numericInput(
                                          inputId = "nSets_ORAplot_microarray_norm",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_ORAplot_microarray_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               plotlyOutput("ORAplot_microarray_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                      tabPanel("Network diagram",
                               icon = icon("fas fa-mouse-pointer"),
                               br(),
                               h3(strong("ORA network diagram")),
                               h5("The ORA network diagram visualize the similarity between the most significant genesets."),
                               hr(),
                               fluidRow(
                                 column(3,
                                        # Number of genesets
                                        numericInput(
                                          inputId = "nSets_network_microarray_norm",
                                          label = "Number of genesets (5-20)",
                                          value = 10,
                                          min = 5,
                                          max = 20,
                                          step = 1)
                                 ),
                                 column(3,
                                        # Network layout
                                        pickerInput(inputId = "layout_network_microarray_norm",
                                                    label = "Network layout",
                                                    choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                                'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                    selected = 'graphopt',
                                                    multiple = FALSE)
                                 )
                               ),
                               hr(),
                               actionBttn(inputId = "plot_ORAnetwork_microarray_norm", 
                                          label = "Plot",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("sync")),
                               br(),
                               br(),
                               plotOutput("ORAnetwork_microarray_norm") %>% 
                                 withSpinner(color="#0dc5c1")
                      ),
                      
                    ) # tabsetpanel
                  ) # taglist
                })
              }# renderUI
            }) # Observe
            
            
            
          }) #observe event
          
          
        } # raw or norm
        
        
      } # microarray or rnaseq
    })
  }) # Observe
  
} # server