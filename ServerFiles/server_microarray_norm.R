observe({
  # Make list for reactive values
  rv <- reactiveValues()
  
  # Download session info
  output$downloadSessionInfo_microarray_norm <- downloadHandler(
    filename = "sessionInfo.txt",
    content = function(file){
      writeLines(capture.output(sessionInfo()), file)
    }
  )
  
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
    
    # Series matrix file
    if (input$ExprDataFileType_microarray_norm == "Series Matrix File"){
      if (!is.null(input$uploadExprData_microarray_norm_smf)){
        
        # Read file
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
    }
    
    # .tsv/.csv file
    if (input$ExprDataFileType_microarray_norm == ".tsv/.csv file"){
      if (!is.null(input$uploadExprData_microarray_norm_tsv)){
        
        # Read file
        m <- readExprMatrix(input$uploadExprData_microarray_norm_tsv$datapath)
        rv$gxData <- Biobase::ExpressionSet(assayData = m)
        
        # Guess the organism
        rv$Organism <- "Homo sapiens"
        
      } else{
        shinyWidgets::sendSweetAlert(
          session = session,
          title = "Error!!",
          text = "You forgot to upload an expression and/or metadata file",
          type = "error")
        shinybusy::remove_modal_spinner()
      }
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
    shinyWidgets::pickerInput(inputId = "groupselect_microarray_norm",
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
                 log-transformation is needed.")
        
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
                   paste0(input$normMeth_microarray_norm,"; ",input$perGroup_microarray_norm),
                   input$transformation_microarray_norm
      )
    )
    
    #------------------------------------------------------------------#
    # Generate output: QC tables/plots
    #------------------------------------------------------------------#
    
    #******************************************************************#
    # Output 1: Table with normalized expression values
    #******************************************************************#
    
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
      
      if (length(levels(rv$experimentFactor)) > 5){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      } else{
        legendColors <- c(input$boxplots_col1_microarray_norm,
                          input$boxplots_col2_microarray_norm,
                          input$boxplots_col3_microarray_norm,
                          input$boxplots_col4_microarray_norm,
                          input$boxplots_col5_microarray_norm)
      }
      if (length(legendColors) != length(levels(rv$experimentFactor))){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      }
      names(legendColors) <- levels(rv$experimentFactor)
      
      getBoxplots(experimentFactor = rv$experimentFactor,
                  legendColors = legendColors,
                  normData = rv$normData,
                  RNASeq = FALSE,
                  width = session$clientData$output_boxplots_microarray_norm_width,
                  height = session$clientData$output_boxplots_microarray_norm_height)
    }, deleteFile = TRUE)
    
    # Allow users to set colors
    observe({
      test <- length(levels(rv$experimentFactor))
      
      if (test == 2){
        output$UI_color_boxplots_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("boxplots_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("boxplots_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2])
          )
        })
      }
      if (test == 3){
        output$UI_color_boxplots_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("boxplots_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("boxplots_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("boxplots_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3])
          )
        })
      }
      if (test == 4){
        output$UI_color_boxplots_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("boxplots_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("boxplots_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("boxplots_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3]),
            colourpicker::colourInput("boxplots_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[4])
          )
        })
        
      }
      if (test == 5){
        output$UI_color_boxplots_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("boxplots_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("boxplots_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("boxplots_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3]),
            colourpicker::colourInput("boxplots_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[4]),
            colourpicker::colourInput("boxplots_col5_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[5])
          )
        })
        
      }
      if (test > 5){
        output$UI_color_boxplots_microarray_norm <- renderUI({
          h5("There are too many experimental groups to select colour manually")
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_boxplots_microarray_norm <- downloadHandler(
      filename = function(){"QC_Boxplots.png"},
      content = function(file){
        png(file,
            width=input$width_boxplots_microarray_norm*3,
            height=input$height_boxplots_microarray_norm*3,
            pointsize=24,
            res = 300)
        
        if (length(levels(rv$experimentFactor)) > 5){
          legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
        } else{
          legendColors <- c(input$boxplots_col1_microarray_norm,
                            input$boxplots_col2_microarray_norm,
                            input$boxplots_col3_microarray_norm,
                            input$boxplots_col4_microarray_norm,
                            input$boxplots_col5_microarray_norm)
        }
        if (length(legendColors) != length(levels(rv$experimentFactor))){
          legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
        }
        names(legendColors) <- levels(rv$experimentFactor)
        
        getBoxplots_download(experimentFactor = rv$experimentFactor,
                             legendColors = legendColors,
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
    
    # Density plots of all genes together
    output$densityplots_microarray_norm <- plotly::renderPlotly({
      
      if (length(levels(rv$experimentFactor)) > 5){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      } else{
        legendColors <- c(input$density_col1_microarray_norm,
                          input$density_col2_microarray_norm,
                          input$density_col3_microarray_norm,
                          input$density_col4_microarray_norm,
                          input$density_col5_microarray_norm)
      }
      if (length(legendColors) != length(levels(rv$experimentFactor))){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      }
      names(legendColors) <- levels(rv$experimentFactor)
      
      rv$density <- getDensityplots(experimentFactor = rv$experimentFactor,
                                    legendColors = legendColors,
                                    normMatrix = rv$normMatrix,
                                    RNASeq = FALSE)
      return(rv$density)
      
    })
    
    # Allow users to set colors
    observe({
      test <- length(levels(rv$experimentFactor))
      
      if (test == 2){
        output$UI_color_density_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("density_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("density_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2])
          )
        })
      }
      if (test == 3){
        output$UI_color_density_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("density_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("density_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("density_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3])
          )
        })
      }
      if (test == 4){
        output$UI_color_density_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("density_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("density_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("density_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3]),
            colourpicker::colourInput("density_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[4])
          )
        })
        
      }
      if (test == 5){
        output$UI_color_density_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("density_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[1]),
            colourpicker::colourInput("density_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[2]),
            colourpicker::colourInput("density_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[3]),
            colourpicker::colourInput("density_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[4]),
            colourpicker::colourInput("density_col5_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(rv$experimentFactor)$legendColors[5])
          )
        })
        
      }
      if (test > 5){
        output$UI_color_density_microarray_norm <- renderUI({
          h5("There are too many experimental groups to select colour manually")
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_density_microarray_norm <- downloadHandler(
      filename = function(){ifelse(input$static_density_microarray_norm, "QC_Density.png", "QC_Density.html")},
      content = function(file){
        
        if (input$static_density_microarray_norm){
          
          if (length(levels(rv$experimentFactor)) > 5){
            legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
          } else{
            legendColors <- c(input$density_col1_microarray_norm,
                              input$density_col2_microarray_norm,
                              input$density_col3_microarray_norm,
                              input$density_col4_microarray_norm,
                              input$density_col5_microarray_norm)
          }
          if (length(legendColors) != length(levels(rv$experimentFactor))){
            legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
          }
          names(legendColors) <- levels(rv$experimentFactor)
          
          p <- getDensityplots_static(experimentFactor = rv$experimentFactor,
                                      legendColors = legendColors,
                                      normMatrix = rv$normMatrix,
                                      RNASeq = FALSE)
          ggplot2::ggsave(plot = p, 
                          filename = file,
                          width = input$width_density_microarray_norm,
                          height = input$height_density_microarray_norm,
                          units = "px")
        } else{
          htmlwidgets::saveWidget(rv$density, 
                                  file)
        }
      }
    )
    
    
    # Make modal
    observeEvent(input$download_density_microarray_norm, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   shinyWidgets::materialSwitch(
                     inputId = "static_density_microarray_norm",
                     label = "Click to make static plot",
                     value = FALSE, 
                     status = "primary"))
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.static_density_microarray_norm==true",
                     sliderInput("height_density_microarray_norm", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.static_density_microarray_norm==true",
                     sliderInput("width_density_microarray_norm", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_density_microarray_norm', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    
    #********************************************************************#
    # Output 4: Heatmap
    #********************************************************************#
    
    output$heatmap_microarray_norm  <- plotly::renderPlotly({
      
      # Make color factor
      if(length(input$colorFactor_heatmap_microarray_norm) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm])
      }
      
      # Set colors
      if (length(levels(colorFactor)) > 5){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      } else{
        legendColors <- c(input$heatmap_col1_microarray_norm,
                          input$heatmap_col2_microarray_norm,
                          input$heatmap_col3_microarray_norm,
                          input$heatmap_col4_microarray_norm,
                          input$heatmap_col5_microarray_norm)
      }
      if (length(legendColors) != length(levels(colorFactor))){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      }
      names(legendColors) <- levels(colorFactor)
      
      # Make heatmap
      rv$heatmap <- getHeatmap(experimentFactor = colorFactor,
                               legendColors = legendColors,
                               normMatrix = rv$normMatrix,
                               clusterOption1 = input$clusteroption1_microarray_norm,
                               clusterOption2 = input$clusteroption2_microarray_norm,
                               theme = input$heatmaptheme_microarray_norm)
    })
    
    # Allow users to set colors
    observe({
      req(input$colorFactor_heatmap_microarray_norm)
      if(length(input$colorFactor_heatmap_microarray_norm) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm])
      }
      test <- length(levels(colorFactor))
      
      if (test == 2){
        output$UI_color_heatmap_microarray_norm <- renderUI({
          tagList(
            h4("Side color"),
            colourpicker::colourInput("heatmap_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("heatmap_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2])
          )
        })
      }
      
      if (test == 3){
        output$UI_color_heatmap_microarray_norm <- renderUI({
          tagList(
            h4("Side color"),
            colourpicker::colourInput("heatmap_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("heatmap_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("heatmap_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3])
          )
        })
      }
      if (test == 4){
        output$UI_color_heatmap_microarray_norm <- renderUI({
          tagList(
            h4("Side color"),
            colourpicker::colourInput("heatmap_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("heatmap_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("heatmap_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3]),
            colourpicker::colourInput("heatmap_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[4])
          )
        })
        
      }
      if (test == 5){
        output$UI_color_heatmap_microarray_norm <- renderUI({
          tagList(
            h4("Side color"),
            colourpicker::colourInput("heatmap_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("heatmap_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("heatmap_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3]),
            colourpicker::colourInput("heatmap_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[4]),
            colourpicker::colourInput("heatmap_col5_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[5])
          )
        })
        
      }
      if (test > 5){
        output$UI_color_heatmap_microarray_norm <- renderUI({
          NULL
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_heatmap_microarray_norm <- downloadHandler(
      filename = function(){ifelse(input$static_heatmap_microarray_norm, "QC_Heatmap.png", "QC_Heatmap.html")},
      content = function(file){
        
        if (input$static_heatmap_microarray_norm){
          
          # Make color factor
          if(length(input$colorFactor_heatmap_microarray_norm) > 1){
            colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm], 1, paste, collapse = "_" ))
          } else{
            colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_norm])
          }
          
          # Set colors
          if (length(levels(colorFactor)) > 5){
            legendColors <- colorsByFactor(colorFactor)$legendColors
          } else{
            legendColors <- c(input$heatmap_col1_microarray_norm,
                              input$heatmap_col2_microarray_norm,
                              input$heatmap_col3_microarray_norm,
                              input$heatmap_col4_microarray_norm,
                              input$heatmap_col5_microarray_norm)
          }
          if (length(legendColors) != length(levels(colorFactor))){
            legendColors <- colorsByFactor(colorFactor)$legendColors
          }
          names(legendColors) <- levels(colorFactor)
          
          # Make heatmap
          getHeatmap_static(experimentFactor = colorFactor,
                            legendColors = legendColors,
                            normMatrix = rv$normMatrix,
                            clusterOption1 = input$clusteroption1_microarray_norm,
                            clusterOption2 = input$clusteroption2_microarray_norm,
                            theme = input$heatmaptheme_microarray_norm,
                            width = input$width_heatmap_microarray_norm,
                            height = input$height_heatmap_microarray_norm,
                            file)
        } else{
          htmlwidgets::saveWidget(rv$heatmap, 
                                  file)
        }
      }
    )
    
    
    # Make modal
    observeEvent(input$download_heatmap_microarray_norm, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   shinyWidgets::materialSwitch(
                     inputId = "static_heatmap_microarray_norm",
                     label = "Click to make static plot",
                     value = FALSE, 
                     status = "primary"))
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.static_heatmap_microarray_norm==true",
                     sliderInput("height_heatmap_microarray_norm", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.static_heatmap_microarray_norm==true",
                     sliderInput("width_heatmap_microarray_norm", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_heatmap_microarray_norm', 
                                  'Download')
            )
          )
          
        )
        
      ))
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
    output$PCA_microarray_norm <- plotly::renderPlotly({
      
      # Get factor to color by
      if(length(input$colorFactor_PCA_microarray_norm) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm])
      }
      
      # Set colors
      if (length(levels(colorFactor)) > 5){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      } else{
        legendColors <- c(input$PCA_col1_microarray_norm,
                          input$PCA_col2_microarray_norm,
                          input$PCA_col3_microarray_norm,
                          input$PCA_col4_microarray_norm,
                          input$PCA_col5_microarray_norm)
      }
      if (length(legendColors) != length(levels(colorFactor))){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      }
      
      # Make PCA score plot
      rv$PCAplot <- plot_PCA(PC_data = rv$PCA_data, 
                             colorFactor = colorFactor,
                             legendColors = legendColors, 
                             xpc = as.numeric(stringr::str_remove(input$xpca_microarray_norm,"PC")), 
                             ypc = as.numeric(stringr::str_remove(input$ypca_microarray_norm,"PC")), 
                             zpc = ifelse(input$xyz_microarray_norm,as.numeric(stringr::str_remove(input$zpca_microarray_norm,"PC")),3), 
                             xyz = input$xyz_microarray_norm)
      return(rv$PCAplot)
    })
    
    
    # Allow users to set colors
    observe({
      req(input$colorFactor_PCA_microarray_norm)
      if(length(input$colorFactor_PCA_microarray_norm) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm])
      }
      test <- length(levels(colorFactor))
      
      if (test == 2){
        output$UI_color_PCA_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("PCA_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("PCA_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2])
          )
        })
      }
      
      if (test == 3){
        output$UI_color_PCA_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("PCA_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("PCA_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("PCA_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3])
          )
        })
      }
      if (test == 4){
        output$UI_color_PCA_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("PCA_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("PCA_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("PCA_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3]),
            colourpicker::colourInput("PCA_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[4])
          )
        })
        
      }
      if (test == 5){
        output$UI_color_PCA_microarray_norm <- renderUI({
          tagList(
            h4("Colors"),
            colourpicker::colourInput("PCA_col1_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[1]),
            colourpicker::colourInput("PCA_col2_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[2]),
            colourpicker::colourInput("PCA_col3_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[3]),
            colourpicker::colourInput("PCA_col4_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[4]),
            colourpicker::colourInput("PCA_col5_microarray_norm", 
                                      NULL, 
                                      colorsByFactor(colorFactor)$legendColors[5])
          )
        })
        
      }
      if (test > 5){
        output$UI_color_PCA_microarray_norm <- renderUI({
          NULL
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_pca_microarray_norm <- downloadHandler(
      filename = function(){ifelse(input$static_pca_microarray_norm, "QC_PCA.png", "QC_PCA.html")},
      content = function(file){
        
        if (input$static_pca_microarray_norm){
          
          # Get factor to color by
          if(length(input$colorFactor_PCA_microarray_norm) > 1){
            colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm], 1, paste, collapse = "_" ))
          } else{
            colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_norm])
          }
          
          # Set colors
          if (length(levels(colorFactor)) > 5){
            legendColors <- colorsByFactor(colorFactor)$legendColors
          } else{
            legendColors <- c(input$PCA_col1_microarray_norm,
                              input$PCA_col2_microarray_norm,
                              input$PCA_col3_microarray_norm,
                              input$PCA_col4_microarray_norm,
                              input$PCA_col5_microarray_norm)
          }
          if (length(legendColors) != length(levels(colorFactor))){
            legendColors <- colorsByFactor(colorFactor)$legendColors
          }
          
          # Make PCA score plot
          p <- plot_PCA_static(PC_data = rv$PCA_data, 
                               colorFactor = colorFactor,
                               legendColors = legendColors, 
                               xpc = as.numeric(stringr::str_remove(input$xpca_microarray_norm,"PC")), 
                               ypc = as.numeric(stringr::str_remove(input$ypca_microarray_norm,"PC")))
          
          ggplot2::ggsave(plot = p, 
                          filename = file,
                          width = input$width_pca_microarray_norm,
                          height = input$height_pca_microarray_norm,
                          units = "px")
        } else{
          htmlwidgets::saveWidget(rv$PCAplot, 
                                  file)
        }
      }
    )
    
    
    # Make modal
    observeEvent(input$download_pca_microarray_norm, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   shinyWidgets::materialSwitch(
                     inputId = "static_pca_microarray_norm",
                     label = "Click to make static plot",
                     value = FALSE, 
                     status = "primary"))
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.static_pca_microarray_norm==true",
                     sliderInput("height_pca_microarray_norm", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.static_pca_microarray_norm==true",
                     sliderInput("width_pca_microarray_norm", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_pca_microarray_norm', 
                                  'Download')
            )
          )
          
        )
        
      ))
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
    # QC report
    #********************************************************************#
    output$QCreport_microarray_norm <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "QCreport.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "QCreport_microarray_norm.Rmd")
        file.copy("Reports/QCreport_microarray_norm.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(processingSettings = rv$processingSettings,
                       experimentFactor = rv$experimentFactor,
                       legendColors = colorsByFactor(rv$experimentFactor)$legendColors,
                       normData = rv$normData,
                       normMatrix = rv$normMatrix,
                       PCAData = rv$PCA_data)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
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
                   h5("Here you can view the normalized and log-transformed intensities. 
                              Click on the table to explore the data!"),
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
                   br(),
                   br(),
                   # customize heatmap
                   shinyWidgets::dropdownButton(
                     uiOutput("UI_color_boxplots_microarray_norm"),
                     circle = TRUE, status = "info",
                     icon = icon("fas fa-cog"), width = "300px",
                     tooltip = shinyWidgets::tooltipOptions(
                       title = "Click to change colors!")
                   ),
                   h2(strong("Boxplot of normalized intensities"), align = "center"),
                   h4("Distributions should be comparable between samples", align = "center"),
                   plotOutput(outputId = "boxplots_microarray_norm",
                              width = "65vw", height = "80vw")%>% 
                     shinycssloaders::withSpinner(color="#0dc5c1")
          ),
          
          # TAB3: Density plots
          tabPanel("Density plots",
                   icon = icon("fas fa-mouse-pointer"),
                   br(),
                   actionButton('download_density_microarray_norm', 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   # customize heatmap
                   shinyWidgets::dropdownButton(
                     uiOutput("UI_color_density_microarray_norm"),
                     circle = TRUE, status = "info",
                     icon = icon("fas fa-cog"), width = "300px",
                     tooltip = shinyWidgets::tooltipOptions(
                       title = "Click to change colors!")
                   ),
                   h2(strong("Density plot of normalized intensities"), align = "center"),
                   h4("Distributions should be comparable between samples", align = "center"),
                   plotly::plotlyOutput(outputId = "densityplots_microarray_norm",
                                        width = "65vw", height = "40vw") %>% 
                     shinycssloaders::withSpinner(color="#0dc5c1")
          ),
          
          # TAB4: Heatmap of sample-sample correlations
          tabPanel("Correlation plot", 
                   icon = icon("fas fa-mouse-pointer"),
                   h3(strong("Heatmap of sample-sample correlations")),
                   hr(),
                   fluidRow(
                     column(3,
                            shinyWidgets::pickerInput(
                              inputId = "colorFactor_heatmap_microarray_norm",
                              label = "Side color",
                              choices = colnames(rv$metaData_fil),
                              selected = rv$experimentName,
                              multiple = TRUE,
                              options = list(
                                style = "btn-info"))
                     ),
                     # Select distance method
                     column(3,
                            shinyWidgets::pickerInput(
                              inputId = "clusteroption1_microarray_norm",
                              label = "Distance",
                              choices = c("Pearson","Spearman",
                                          "Euclidean"))
                     ),
                     
                     #Clustering method
                     column(3,
                            shinyWidgets::pickerInput(
                              inputId = "clusteroption2_microarray_norm",
                              label = "Clustering",
                              choices = c("ward.D2","single",
                                          "complete","average",
                                          "mcquitty","median",
                                          "centroid"))
                     )
                   ),
                   hr(),
                   actionButton('download_heatmap_microarray_norm', 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   
                   # customize heatmap
                   shinyWidgets::dropdownButton(
                     h4("Heatmap theme"),
                     selectInput(inputId = 'heatmaptheme_microarray_norm',
                                 label = NULL,
                                 choices = c("Default", 
                                             "Yellow-red", 
                                             "Blues", 
                                             "Reds")),
                     br(),
                     uiOutput("UI_color_heatmap_microarray_norm"),
                     
                     circle = TRUE, status = "info",
                     icon = icon("fas fa-cog"), width = "300px",
                     tooltip = shinyWidgets::tooltipOptions(
                       title = "Click to change colors!")
                   ),
                   
                   # Make heatmap
                   plotly::plotlyOutput("heatmap_microarray_norm", 
                                        width = "65vw", height = "35vw") %>% 
                     shinycssloaders::withSpinner(color="#0dc5c1", 
                                                  proxy.height = "400px")
                   
          ),
          
          # TAB5: PCA score plot
          tabPanel("PCA",
                   icon = icon("fas fa-mouse-pointer"),
                   # Title + text
                   h3(strong("Principal Component Analysis (PCA)")),
                   shinyWidgets::materialSwitch(
                     inputId = "xyz_microarray_norm",
                     label = "3D",
                     value = FALSE, 
                     status = "info"),
                   
                   hr(),
                   
                   # Set color + 3D/2D
                   fluidRow(
                     column(3,
                            # Color by which factor?
                            shinyWidgets::pickerInput(inputId = "colorFactor_PCA_microarray_norm",
                                                      label = "Color by:",
                                                      choices = colnames(rv$metaData_fil),
                                                      selected = rv$experimentName,
                                                      multiple = TRUE,
                                                      options = list(
                                                        style = "btn-info"))
                     ),
                     column(3,
                            #X-axis
                            shinyWidgets::pickerInput(inputId = "xpca_microarray_norm", 
                                                      label = "x-axis",
                                                      choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                                  "PC6", "PC7", "PC8"),
                                                      selected = "PC1")
                     ),
                     column(3,
                            #Y-axis
                            shinyWidgets::pickerInput(inputId = "ypca_microarray_norm", 
                                                      label = "y-axis",
                                                      choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                                  "PC6", "PC7", "PC8"),
                                                      selected = "PC2")
                     ),
                     column(3,
                            #Z-axis
                            conditionalPanel(
                              condition = "input.xyz_microarray_norm==true",
                              shinyWidgets::pickerInput(inputId = "zpca_microarray_norm", 
                                                        label = "z-axis",
                                                        choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                                    "PC6", "PC7", "PC8"),
                                                        selected = "PC3")
                            )
                     )
                   ),
                   
                   hr(),
                   actionButton("download_pca_microarray_norm", 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   
                   fluidRow(
                     
                     # Customize PCA plot
                     shinyWidgets::dropdownButton(
                       
                       uiOutput("UI_color_PCA_microarray_norm"),
                       circle = TRUE, status = "info",
                       icon = icon("fas fa-cog"), width = "300px",
                       tooltip = shinyWidgets::tooltipOptions(
                         title = "Click to change colors!")
                     ),
                     
                     # Make PCA plot
                     plotly::plotlyOutput("PCA_microarray_norm",
                                          width = "55vw", height = "30vw")%>% 
                       shinycssloaders::withSpinner(color="#0dc5c1"),
                   ),
                   
                   
          ), # EO PCA tabpanel
          
          # # TAB6: settings table
          tabPanel("Settings overview",
                   icon = icon("fas fa-file"),
                   h3(strong("Pre-processing settings")),
                   h5("Here you can see an overview of the chosen pre-processing settings."),
                   hr(),
                   DT::dataTableOutput(outputId = "processingSettings_microarray_norm") %>% 
                     withSpinner(color="#0dc5c1"),
                   br(),
                   downloadButton("downloadProcessingSettings_microarray_norm", 
                                  "Download table"),
                   downloadButton("downloadSessionInfo_microarray_norm", 
                                  "Session info")
                   
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
        shinyWidgets::downloadBttn(outputId = "QCreport_microarray_norm",
                                   label = "Get report",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("download")),
        shinyWidgets::actionBttn(inputId = "next_preprocessing_microarray_norm",
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
      
      # Make top table
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
      
      # Make settings table
      rv$statSettings <- list()
      for (c in 1:length(input$comparisons_microarray_norm)){
        rv$statSettings[[c]] <- data.frame(
          Option = c("Selected comparison",
                     "Continuous covariate(s)",
                     "Discrete covariate(s)",
                     "Gene annotation dataset",
                     "Gene annotation attribute(s)",
                     "Gene annotation filter"),
          Selected = c(input$comparisons_microarray_norm[c],
                       ifelse(is.null(input$covGroups_num_microarray_norm), " ",
                              paste(input$covGroups_num_microarray_norm, collapse = "; ")),
                       ifelse(is.null(input$covGroups_char_microarray_norm), " ",
                              paste(input$covGroups_char_microarray_norm, collapse = "; ")),
                       rv$top_table_list[[3]],
                       paste(input$biomart_attributes_microarray_norm, collapse = "; "),
                       input$biomart_filter_microarray_norm
          )
        )
      }
      names(rv$statSettings) <- input$comparisons_microarray_norm
      
    } else{
      
      # Make top table
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
      
      # Make settings table
      for (c in 1:length(input$comparisons_microarray_norm)){
        rv$statSettings[[c]] <- data.frame(
          Option = c("Comparison",
                     "Continuous covariate(s)",
                     "Discrete covariate(s)"),
          Selected = c(input$comparisons_microarray_norm[c],
                       ifelse(is.null(input$covGroups_num_microarray_norm), " ",
                              paste(input$covGroups_num_microarray_norm, collapse = "; ")),
                       ifelse(is.null(input$covGroups_char_microarray_norm), " ",
                              paste(input$covGroups_char_microarray_norm, collapse = "; ")))
        )
      }
      names(rv$statSettings) <- input$comparisons_microarray_norm
      
    }
    
    # Select comparison for output
    observe({
      rv$top_table <- rv$top_table_list[[1]]
      
      if (!is.null(rv$top_table_list)){
        
        # Remove modal
        shinybusy::remove_modal_spinner()
        
        # Show comparisons
        output$UI_comparisons_view_microarray_norm <- renderUI({
          shinyWidgets::pickerInput(inputId = "comparisons_view_microarray_norm",
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
        showTab("navbar", target = "panel_ORA_microarray_norm")
        
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
    
    
    # Make outputs:
    
    #********************************************************************#
    # top table
    #********************************************************************#
    
    # print top table
    output$top_table_microarray_norm <- DT::renderDataTable({
      
      # Required: selected comparison
      req(input$comparisons_view_microarray_norm)
      
      # Get statistics of selected comparison
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
    # TAB2: P value and logFC histograms
    #********************************************************************#
    
    observe({
      req(input$comparisons_view_microarray_norm)
      req(rv$top_table)
      
      if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
        
        # Make p-value histogram
        rv$Phistogram <- makePHistogram(P = rv$top_table[[input$comparisons_view_microarray_norm]][,"p-value"],
                                        color = input$histogram_color_microarray_norm,
                                        bins = input$histogram_bins_microarray_norm)
        
        # Plot p-value histogram
        output$Phistogram_microarray_norm <- plotly::renderPlotly({
          return(rv$Phistogram)
        })
        
        #***************************#
        # Modal to P value histogram
        #***************************#
        
        # Download plot
        output$realdownload_Phistogram_microarray_norm <- downloadHandler(
          filename = function(){ifelse(input$static_Phistogram_microarray_norm, "Phistogram.png", "Phistogram.html")},
          content = function(file){
            
            if (input$static_Phistogram_microarray_norm){
              
              # Make volcano plot
              p <- makePHistogram(P = rv$top_table[[input$comparisons_view_microarray_norm]][,"p-value"],
                                  color = input$histogram_color_microarray_norm,
                                  bins = input$histogram_bins_microarray_norm,
                                  static = TRUE)
              
              ggplot2::ggsave(plot = p, 
                              filename = file,
                              width = input$width_Phistogram_microarray_norm,
                              height = input$height_Phistogram_microarray_norm,
                              units = "px")
            } else{
              htmlwidgets::saveWidget(rv$Phistogram, 
                                      file)
            }
          }
        )
        
        
        # Make modal
        observeEvent(input$download_Phistogram_microarray_norm, {
          showModal(modalDialog(
            title = NULL,
            easyClose = TRUE,
            size = "m",
            footer = tagList(
              fluidRow(
                column(12, align = "left",
                       shinyWidgets::materialSwitch(
                         inputId = "static_Phistogram_microarray_norm",
                         label = "Click to make static plot",
                         value = FALSE, 
                         status = "primary"))
              ),
              fluidRow(
                column(6,
                       conditionalPanel(
                         condition = "input.static_Phistogram_microarray_norm==true",
                         sliderInput("height_Phistogram_microarray_norm", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                       )
                ),
                column(6,
                       conditionalPanel(
                         condition = "input.static_Phistogram_microarray_norm==true",
                         sliderInput("width_Phistogram_microarray_norm", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                       )
                )
              ),
              
              fluidRow(
                column(12, align = "left",
                       downloadButton('realdownload_Phistogram_microarray_norm', 
                                      'Download')
                )
              )
              
            )
            
          ))
        })
        
        
        # Make logFC histogram
        rv$logFChistogram <- makelogFCHistogram(logFC = rv$top_table[[input$comparisons_view_microarray_norm]][,"log2FC"],
                                                color = input$histogram_color_microarray_norm,
                                                bins = input$histogram_bins_microarray_norm)
        
        # Plot logFC histogram
        output$logFChistogram_microarray_norm <- plotly::renderPlotly({
          return(rv$logFChistogram)
        })
        
        #********************************#
        # Modal to logFC value histogram
        #********************************#
        
        # Download plot
        output$realdownload_logFChistogram_microarray_norm <- downloadHandler(
          filename = function(){ifelse(input$static_logFChistogram_microarray_norm, "logFChistogram.png", "logFChistogram.html")},
          content = function(file){
            
            if (input$static_logFChistogram_microarray_norm){
              
              # Make volcano plot
              p <- makelogFCHistogram(logFC = rv$top_table[[input$comparisons_view_microarray_norm]][,"log2FC"],
                                      color = input$histogram_color_microarray_norm,
                                      bins = input$histogram_bins_microarray_norm,
                                      static = TRUE)
              
              ggplot2::ggsave(plot = p, 
                              filename = file,
                              width = input$width_logFChistogram_microarray_norm,
                              height = input$height_logFChistogram_microarray_norm,
                              units = "px")
            } else{
              htmlwidgets::saveWidget(rv$logFChistogram, 
                                      file)
            }
          }
        )
        
        
        # Make modal
        observeEvent(input$download_logFChistogram_microarray_norm, {
          showModal(modalDialog(
            title = NULL,
            easyClose = TRUE,
            size = "m",
            footer = tagList(
              fluidRow(
                column(12, align = "left",
                       shinyWidgets::materialSwitch(
                         inputId = "static_logFChistogram_microarray_norm",
                         label = "Click to make static plot",
                         value = FALSE, 
                         status = "primary"))
              ),
              fluidRow(
                column(6,
                       conditionalPanel(
                         condition = "input.static_logFChistogram_microarray_norm==true",
                         sliderInput("height_logFChistogram_microarray_norm", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                       )
                ),
                column(6,
                       conditionalPanel(
                         condition = "input.static_logFChistogram_microarray_norm==true",
                         sliderInput("width_logFChistogram_microarray_norm", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                       )
                )
              ),
              
              fluidRow(
                column(12, align = "left",
                       downloadButton('realdownload_logFChistogram_microarray_norm', 
                                      'Download')
                )
              )
              
            )
            
          ))
        })
        
        
      }
    })
    
    #********************************************************************#
    # TAB3: Volcano plot
    #********************************************************************#
    
    observeEvent(input$plot_volcano_microarray_norm, {
      req(rv$top_table)
      req(input$rawp_volcano_microarray_norm)
      req(input$p_thres_volcano_microarray_norm)
      req(input$logFC_thres_volcano_microarray_norm)
      req(input$comparisons_view_microarray_norm)
      
      if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
        rv$volcano <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                                  p = input$rawp_volcano_microarray_norm, 
                                  p_threshold = input$p_thres_volcano_microarray_norm, 
                                  logFC_threshold = input$logFC_thres_volcano_microarray_norm,
                                  unchanged_color = input$volcano_unchanged_color_microarray_norm,
                                  down_color = input$volcano_down_color_microarray_norm,
                                  up_color = input$volcano_up_color_microarray_norm
                                  )
        
        output$volcano_microarray_norm <- plotly::renderPlotly(rv$volcano)
      }
    }, ignoreNULL = FALSE) 
    # ignoreNULL: generate plot even if action button is not pressed
    
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_volcano_microarray_norm <- downloadHandler(
      filename = function(){ifelse(input$static_volcano_microarray_norm, "Volcano.png", "Volcano.html")},
      content = function(file){
        
        if (input$static_volcano_microarray_norm){
          
          
          # Make PCA score plot
          p <- makeVolcano_static(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                                  p = input$rawp_volcano_microarray_norm, 
                                  p_threshold = input$p_thres_volcano_microarray_norm, 
                                  logFC_threshold = input$logFC_thres_volcano_microarray_norm,
                                  unchanged_color = input$volcano_unchanged_color_microarray_norm,
                                  down_color = input$volcano_down_color_microarray_norm,
                                  up_color = input$volcano_up_color_microarray_norm)
          
          ggplot2::ggsave(plot = p, 
                          filename = file,
                          width = input$width_volcano_microarray_norm,
                          height = input$height_volcano_microarray_norm,
                          units = "px")
        } else{
          htmlwidgets::saveWidget(rv$volcano, 
                                  file)
        }
      }
    )
    
    
    # Make modal
    observeEvent(input$download_volcano_microarray_norm, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   shinyWidgets::materialSwitch(
                     inputId = "static_volcano_microarray_norm",
                     label = "Click to make static plot",
                     value = FALSE, 
                     status = "primary"))
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.static_volcano_microarray_norm==true",
                     sliderInput("height_volcano_microarray_norm", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.static_volcano_microarray_norm==true",
                     sliderInput("width_volcano_microarray_norm", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_volcano_microarray_norm', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    #********************************************************************#
    # TAB4: MA plot
    #********************************************************************#
    
    observeEvent(input$plot_MA_microarray_norm, {
      req(rv$top_table)
      req(input$rawp_MA_microarray_norm)
      req(input$p_thres_MA_microarray_norm)
      req(input$logFC_thres_MA_microarray_norm)
      req(input$comparisons_view_microarray_norm)
      
      if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
        rv$MA <- makeMAplot(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                            p = input$rawp_MA_microarray_norm, 
                            p_threshold = input$p_thres_MA_microarray_norm, 
                            logFC_threshold = input$logFC_thres_MA_microarray_norm,
                            unchanged_color = input$MA_unchanged_color_microarray_norm,
                            down_color = input$MA_down_color_microarray_norm,
                            up_color = input$MA_up_color_microarray_norm,
                            RNAseq = FALSE)
        
        output$MA_microarray_norm <- plotly::renderPlotly(rv$MA)
      }
    }, ignoreNULL = FALSE) 
    # ignoreNULL: generate plot even if action button is not pressed
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    output$realdownload_MA_microarray_norm <- downloadHandler(
      filename = function(){ifelse(input$static_MA_microarray_norm, "MA.png", "MA.html")},
      content = function(file){
        
        if (input$static_MA_microarray_norm){
          
          
          # Make PCA score plot
          p <- makeMAplot_static(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                                 p = input$rawp_MA_microarray_norm, 
                                 p_threshold = input$p_thres_MA_microarray_norm, 
                                 logFC_threshold = input$logFC_thres_MA_microarray_norm,
                                 unchanged_color = input$MA_unchanged_color_microarray_norm,
                                 down_color = input$MA_down_color_microarray_norm,
                                 up_color = input$MA_up_color_microarray_norm,
                                 RNAseq = FALSE)
          
          ggplot2::ggsave(plot = p, 
                          filename = file,
                          width = input$width_MA_microarray_norm,
                          height = input$height_MA_microarray_norm,
                          units = "px")
        } else{
          htmlwidgets::saveWidget(rv$MA, 
                                  file)
        }
      }
    )
    
    
    # Make modal
    observeEvent(input$download_MA_microarray_norm, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(12, align = "left",
                   shinyWidgets::materialSwitch(
                     inputId = "static_MA_microarray_norm",
                     label = "Click to make static plot",
                     value = FALSE, 
                     status = "primary"))
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.static_MA_microarray_norm==true",
                     sliderInput("height_MA_microarray_norm", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.static_MA_microarray_norm==true",
                     sliderInput("width_MA_microarray_norm", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_MA_microarray_norm', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    #********************************************************************#
    # TAB5: Summary
    #********************************************************************#
    observeEvent(input$get_summary_microarray_norm, {
      req(rv$top_table)
      req(input$rawp_summary_microarray_norm)
      req(input$p_thres_summary_microarray_norm)
      req(input$logFC_thres_summary_microarray_norm)
      req(input$comparisons_view_microarray_norm)
      
      if (input$comparisons_view_microarray_norm %in% names(rv$top_table)){
        # DEG table
        temp <- matrix(NA,2,3)
        
        if (input$rawp_summary_microarray_norm == "raw"){
          rownames(temp) <- c(paste0("p-value < ", input$p_thres_summary_microarray_norm),
                              paste0("p-value > ", input$p_thres_summary_microarray_norm))
          colnames(temp) <- c(paste0("log2FC < -", input$logFC_thres_summary_microarray_norm),
                              paste0("-",input$logFC_thres_summary_microarray_norm, " < log2FC < ", input$logFC_thres_summary_microarray_norm),
                              paste0("log2FC > ", input$logFC_thres_summary_microarray_norm))
          
          temp[1,1] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` < input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC < -1*input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[1,2] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` < input$p_thres_summary_microarray_norm) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC) < input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[1,3] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` < input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC > input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,1] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` > input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC < -1*input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,2] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` > input$p_thres_summary_microarray_norm) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC) < input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,3] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` > input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC > input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
        } else{
          rownames(temp) <- c(paste0("adj. p-value < ", input$p_thres_summary_microarray_norm),
                              paste0("adj. p-value > ", input$p_thres_summary_microarray_norm))
          colnames(temp) <- c(paste0("log2FC < -", input$logFC_thres_summary_microarray_norm),
                              paste0("-",input$logFC_thres_summary_microarray_norm, " < log2FC < ", input$logFC_thres_summary_microarray_norm),
                              paste0("log2FC > ", input$logFC_thres_summary_microarray_norm))
          
          temp[1,1] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`adj. p-value` < input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC < -1*input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[1,2] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`adj. p-value` < input$p_thres_summary_microarray_norm) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC) < input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[1,3] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`adj. p-value` < input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC > input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,1] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`adj. p-value` > input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC < -1*input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,2] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`p-value` > input$p_thres_summary_microarray_norm) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC) < input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
          temp[2,3] <- sum((rv$top_table[[input$comparisons_view_microarray_norm]]$`adj. p-value` > input$p_thres_summary_microarray_norm) &
                             (rv$top_table[[input$comparisons_view_microarray_norm]]$log2FC > input$logFC_thres_summary_microarray_norm), na.rm = TRUE)
        }
        rv$summaryTable <- temp
        rm(temp)
        
        output$summaryTable_microarray_norm <- DT::renderDataTable({
          return(rv$summaryTable)
        },options = list(pageLength = 10,
                         dom = 't'),
        selection = "none")
        
      }
    }, ignoreNULL = FALSE)
    
    
    # Download button
    output$downloadSummaryTable_microarray_norm <- downloadHandler(
      filename = "SummaryTable.csv",
      content = function(file){
        write.csv(rv$summaryTable, file, quote = FALSE, row.names = TRUE)
      }
    )
    
    #********************************************************************#
    # TAB6: Settings
    #********************************************************************#
    observe({
      req(input$comparisons_view_microarray_norm)
      
      # Print table with settings
      output$statSettings_microarray_norm <- DT::renderDataTable({
        return(rv$statSettings[[input$comparisons_view_microarray_norm]])
      },options = list(pageLength = 10,
                       dom = 't',
                       rownames = FALSE),
      selection = "none")
      
      # Download button
      output$downloadStatSettings_microarray_norm <- downloadHandler(
        filename = "StatisticalAnalysis_Settings.csv",
        content = function(file){
          write.csv(rv$statSettings[[input$comparisons_view_microarray_norm]], file, quote = FALSE, row.names = FALSE)
        }
      )
      
    })
    
    #********************************************************************#
    # Statistical analysis report
    #********************************************************************#
    output$SAreport_microarray_norm <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "SAreport.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "SAreport_microarray_norm.Rmd")
        file.copy("Reports/SAreport_microarray_norm.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(statSettings = rv$statSettings[[input$comparisons_view_microarray_norm]],
                       topTable = rv$top_table[[input$comparisons_view_microarray_norm]],
                       volcanoTable = rv$summaryTable)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
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
                       h3(strong("Top table")),
                       h5("The top table includes the output of the statistical analysis. 
                                  Click on the table to explore the data!"),
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
                       # Dropdown Button to adjust the plot settings
                       shinyWidgets::dropdownButton(
                         tags$div(
                           style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                           tags$h4("Select color"),
                           colourpicker::colourInput("histogram_color_microarray_norm", 
                                                     NULL, 
                                                     "#d3d3d3"),
                           br(),
                           tags$h4("Number of bins"),
                           numericInput(inputId = "histogram_bins_microarray_norm",
                                        label = NULL,
                                        value = 100),
                           br(),br()
                         ),
                         circle = TRUE, status = "info",
                         icon = icon("fas fa-cog"),
                         
                         tooltip = shinyWidgets::tooltipOptions(
                           title = "Click to personalize histograms!")
                         
                       ), # EO dropdownButton
                       plotly::plotlyOutput("Phistogram_microarray_norm")%>% 
                         shinycssloaders::withSpinner(color="#0dc5c1"),
                       actionButton("download_Phistogram_microarray_norm", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       hr(),
                       plotly::plotlyOutput("logFChistogram_microarray_norm")%>% 
                         shinycssloaders::withSpinner(color="#0dc5c1"),
                       actionButton("download_logFChistogram_microarray_norm", 
                                    "Download figure",
                                    icon = shiny::icon("download"))
                       
              ),
              
              #********************************************************************#
              # volcano tab
              #********************************************************************#
              
              tabPanel("Volcano plot",
                       icon = icon("fas fa-mouse-pointer"),
                       br(),
                       
                       # Dropdown Button to adjust the plot settings
                       shinyWidgets::dropdownButton(
                         tags$div(
                           style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                           tags$h4("Color of unchanged genes"),
                           colourpicker::colourInput("volcano_unchanged_color_microarray_norm", 
                                                     NULL, 
                                                     "darkgrey"),
                           tags$h4("Color of downregulated genes"),
                           colourpicker::colourInput("volcano_down_color_microarray_norm", 
                                                     NULL, 
                                                     "blue"),
                           tags$h4("Color of upregulated genes"),
                           colourpicker::colourInput("volcano_up_color_microarray_norm", 
                                                     NULL, 
                                                     "red"),
                           br()
                         ),
                         circle = TRUE, status = "info",
                         icon = icon("fas fa-cog"),
                         
                         tooltip = shinyWidgets::tooltipOptions(
                           title = "Click to personalize the volcano plot!")
                         
                       ), # EO dropdownButton
                       
                       # Volcano plot output
                       plotly::plotlyOutput("volcano_microarray_norm")%>% 
                         withSpinner(color="#0dc5c1"),
                       
                       
                       br(),
                       actionButton("download_volcano_microarray_norm", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                prettyRadioButtons(
                                  inputId = "rawp_volcano_microarray_norm",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_volcano_microarray_norm",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "plot_volcano_microarray_norm", 
                                                         label = "Load",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_volcano_microarray_norm",
                                  label = HTML(paste0("log", tags$sub("2"),"FC threshold")),
                                  value = 1,
                                  width = "80%")
                         )
                       ),
                       
                       hr()
                       
              ),
              
              #********************************************************************#
              # MA plot tab
              #********************************************************************#
              
              tabPanel("MA plot",
                       icon = icon("fas fa-mouse-pointer"),
                       br(),
                       
                       # Dropdown Button to adjust the plot settings
                       shinyWidgets::dropdownButton(
                         tags$div(
                           style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                           tags$h4("Color of unchanged genes"),
                           colourpicker::colourInput("MA_unchanged_color_microarray_norm", 
                                                     NULL, 
                                                     "darkgrey"),
                           tags$h4("Color of downregulated genes"),
                           colourpicker::colourInput("MA_down_color_microarray_norm", 
                                                     NULL, 
                                                     "blue"),
                           tags$h4("Color of upregulated genes"),
                           colourpicker::colourInput("MA_up_color_microarray_norm", 
                                                     NULL, 
                                                     "red"),
                           br()
                         ),
                         circle = TRUE, status = "info",
                         icon = icon("fas fa-cog"),
                         
                         tooltip = shinyWidgets::tooltipOptions(
                           title = "Click to personalize the MA plot!")
                         
                       ), # EO dropdownButton
                       
                       # MA plot output
                       plotly::plotlyOutput("MA_microarray_norm")%>% 
                         withSpinner(color="#0dc5c1"),
                       br(),
                       actionButton("download_MA_microarray_norm", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                prettyRadioButtons(
                                  inputId = "rawp_MA_microarray_norm",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_MA_microarray_norm",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "plot_MA_microarray_norm", 
                                                         label = "Load",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_MA_microarray_norm",
                                  label = HTML(paste0("log", tags$sub("2"),"FC threshold")),
                                  value = 1,
                                  width = "80%")
                         )
                       ),
                       
                       hr()
                       
              ),
              
              #********************************************************************#
              # Summary tab
              #********************************************************************#
              tabPanel("Summary",
                       icon = icon("fas fa-file"),
                       h3(strong("Summary table")),
                       h5(HTML(paste0("The table provides a summary of the 
                                      number of genes per p-value and log",tags$sub("2"),"FC category"))),
                       hr(),
                       DT::dataTableOutput(outputId = "summaryTable_microarray_norm") %>%
                         withSpinner(color="#0dc5c1"),
                       br(),
                       downloadButton("downloadSummaryTable_microarray_norm",
                                      "Download table"),
                       
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                prettyRadioButtons(
                                  inputId = "rawp_summary_microarray_norm",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_summary_microarray_norm",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "get_summary_microarray_norm", 
                                                         label = "Calculate",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_summary_microarray_norm",
                                  label = HTML(paste0("log", tags$sub("2"),"FC threshold")),
                                  value = 1,
                                  width = "80%")
                         )
                       ),
                       
                       hr()
              ),
              
              #********************************************************************#
              # Settings tab
              #********************************************************************#
              tabPanel("Settings overview",
                       icon = icon("fas fa-file"),
                       h3(strong("Statistical analysis settings")),
                       h5("To enhance reproducibility, view and download the overview of the chosen statistical analysis settings."),
                       hr(),
                       DT::dataTableOutput(outputId = "statSettings_microarray_norm") %>% 
                         withSpinner(color="#0dc5c1"),
                       br(),
                       downloadButton("downloadStatSettings_microarray_norm", 
                                      "Download table"),
                       downloadButton("downloadSessionInfo_microarray_norm", 
                                      "Session info")
                       
              ) # EO Settings tabPanel
              
              
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
        shinyWidgets::downloadBttn(outputId = "SAreport_microarray_norm",
                                   label = "Get report",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("download")),
        actionBttn(inputId = "next_statistics_microarray_norm",
                   label = "Next",
                   style = "jelly",
                   color = "danger",
                   icon = icon("arrow-right"))
      )
    })
    
  }) # observeEvent
  
  
  #======================================================================#
  # TAB5: ORA/GSEA
  #======================================================================#
  # Go to data upload tab
  observeEvent(input$next_statistics_microarray_norm,{
    
    # Go to microarray statistics tab
    updateNavbarPage(session, "navbar",
                     selected = "panel_ORA_microarray_norm")
  })
  
  #********************************************************************#
  # Options for geneset analysis
  #********************************************************************#
  
  # Comparisons for which geneset analysis should be performed
  observe({
    req(rv$top_table)
    output$UI_comparisons_view_ORA_microarray_norm <- renderUI({
      shinyWidgets::pickerInput(inputId = "comparisons_view_ORA_microarray_norm",
                                label = NULL,
                                choices = names(rv$top_table),
                                selected = names(rv$top_table)[1],
                                multiple = FALSE)
    })
  })
  
  # Gene ID information
  observe({
    
    # Get all columns of the top table that can possibly contain the gene IDs
    #req(input$comparisons_view_ORA_microarray_norm)
    col_choice <- 1
    if (length(colnames(rv$top_table[[1]])) > 6){
      col_choice <- c(1,7:ncol(rv$top_table[[1]]))
    }
    
    # Several options need to be selected before ORA can be performed
    output$UI_geneID_ORA_microarray_norm <- renderUI({
      tagList(
        
        # Which columns of the top table contains the gene ids?
        shinyWidgets::pickerInput(inputId = "geneID_ORA_microarray_norm",
                                  label = "Which column of the top table contains the gene IDs?",
                                  choices = colnames(rv$top_table[[1]])[col_choice],
                                  selected = colnames(rv$top_table[[1]])[1],
                                  multiple = FALSE)
      )
    })
  })
  
  
  observeEvent(input$calculate_ORA_microarray_norm,{
    
    #==========================================================================#
    
    # ORA
    
    #==========================================================================#
    
    if (input$ORA_or_GSEA_microarray_norm == "ORA"){
      # Show modal
      shinybusy::show_modal_spinner(text = "Overrepresentation analysis...",
                                    color="#0dc5c1")
      
      # Perform ORA:
      
      # Perform ORA based on logFC/P value threshold(s)
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
        
        
        rv$ORA_settings <- data.frame(
          Option = c("Comparison",
                     "Geneset",
                     "Gene ID",
                     "Organism",
                     "Method",
                     "p-value threshold",
                     "logFC threshold"
          ),
          Selected = c(input$comparisons_view_ORA_microarray_norm,
                       input$geneset_ORA_microarray_norm,
                       paste0(input$selID_ORA_microarray_norm, " (top table column: ", input$geneID_ORA_microarray_norm, ")"),
                       input$organism_ORA_microarray_norm,
                       paste0("ORA (", input$topNorThres_microarray_norm, "; ", input$updown_ORA_microarray_norm, ")"),
                       paste0(input$p_thres_ORA_microarray_norm, " (", input$rawp_ORA_microarray_norm, ")"),
                       input$logFC_thres_ORA_microarray_norm
          )
        )
        
        # Perform ORA based on top N most significant genes
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
        
        
        rv$ORA_settings <- data.frame(
          Option = c("Comparison",
                     "Geneset",
                     "Gene ID",
                     "Organism",
                     "Method",
                     "Top N"
          ),
          Selected = c(input$comparisons_view_ORA_microarray_norm,
                       input$geneset_ORA_microarray_norm,
                       paste0(input$selID_ORA_microarray_norm, " (top table column: ", input$geneID_ORA_microarray_norm, ")"),
                       input$organism_ORA_microarray_norm,
                       paste0("ORA (top ",input$topN_microarray_norm, "; ", input$updown_ORA_microarray_norm, ")"),
                       input$topN_microarray_norm
          )
        )
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
            text = "Overrepresentation analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
            type = "info")
          
          #--------------------------------------------------------------#
          # ORA statistics table
          #--------------------------------------------------------------#
          
          output$ORA_table_microarray_norm <- DT::renderDataTable({
            req(input$geneset_ORA_microarray_norm)
            req(rv$ORA_data)
            output <- rv$ORA_data@result
            
            # Link to wikipathways website if geneset ID is from WikiPathways
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
            
            # Link to QuickGO website if geneset ID is a GO term
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
          output$ORAgene_table_microarray_norm <- DT::renderDataTable({
            req(input$ORA_table_microarray_norm_rows_selected)
            req(rv$ORA_data)
            
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
            req(rv$ORA_data)
            text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_microarray_norm_rows_selected,"ID"],
                           "</b></h3>")
            return(text)
          })
          
          #--------------------------------------------------------------#
          # ORA barchart
          #--------------------------------------------------------------#
          
          observe({
            req(input$nSets_ORAplot_microarray_norm)
            req(rv$ORA_data)
            rv$ORAplot <- makeORAplot(rv$ORA_data,
                                      nSets = input$nSets_ORAplot_microarray_norm,
                                      color = input$color_ORAplot_microarray_norm)
            
            output$ORAplot_microarray_norm <- plotly::renderPlotly(rv$ORAplot)
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          output$realdownload_ORAplot_microarray_norm <- downloadHandler(
            filename = function(){ifelse(input$static_ORAplot_microarray_norm, "ORA barchart.png", "ORA barchart.html")},
            content = function(file){
              
              if (input$static_ORAplot_microarray_norm){
                
                
                # Make MA plot
                p <- makeORAplot(rv$ORA_data,
                                 nSets = input$nSets_ORAplot_microarray_norm,
                                 color = input$color_ORAplot_microarray_norm,
                                 static = TRUE)
                
                ggplot2::ggsave(plot = p, 
                                filename = file,
                                width = input$width_ORAplot_microarray_norm,
                                height = input$height_ORAplot_microarray_norm,
                                units = "px")
              } else{
                htmlwidgets::saveWidget(rv$ORAplot, 
                                        file)
              }
            }
          )
          
          
          # Make modal
          observeEvent(input$download_ORAplot_microarray_norm, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(12, align = "left",
                         shinyWidgets::materialSwitch(
                           inputId = "static_ORAplot_microarray_norm",
                           label = "Click to make static plot",
                           value = FALSE, 
                           status = "primary"))
                ),
                fluidRow(
                  column(6,
                         conditionalPanel(
                           condition = "input.static_ORAplot_microarray_norm==true",
                           sliderInput("height_ORAplot_microarray_norm", 
                                       "Height",
                                       min = 800, max = 2000,
                                       value = 1200, step = 10,
                                       width = "100%")
                         )
                  ),
                  column(6,
                         conditionalPanel(
                           condition = "input.static_ORAplot_microarray_norm==true",
                           sliderInput("width_ORAplot_microarray_norm", 
                                       "Width",
                                       min = 800, max = 2000,
                                       value = 1500, step = 10,
                                       width = "100%")
                         )
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_ORAplot_microarray_norm', 
                                        'Download')
                  )
                )
                
              )
              
            ))
          })
          
          #--------------------------------------------------------------#
          # ORA network diagram
          #--------------------------------------------------------------#
          
          observe({
            req(input$layout_ORAnetwork_microarray_norm)
            req(input$nSets_ORAnetwork_microarray_norm)
            req(rv$ORA_data)
            rv$ORAnetwork <- makeORAnetwork(ORA_data = rv$ORA_data,
                                            layout = input$layout_ORAnetwork_microarray_norm,
                                            nSets = input$nSets_ORAnetwork_microarray_norm,
                                            color = input$color_ORAnetwork_microarray_norm)
            
            output$ORAnetwork_microarray_norm <- renderPlot(rv$ORAnetwork,
                                                            height = 500, 
                                                            width = 800)
            
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          output$realdownload_ORAnetwork_microarray_norm <- downloadHandler(
            filename = "ORA network.png",
            content = function(file){
              
              ggplot2::ggsave(plot = rv$ORAnetwork, 
                              filename = file,
                              width = input$width_ORAnetwork_microarray_norm*2,
                              height = input$height_ORAnetwork_microarray_norm*2,
                              units = "px")
            }
          )
          
          
          # Make modal
          observeEvent(input$download_ORAnetwork_microarray_norm, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,
                         sliderInput("height_ORAnetwork_microarray_norm", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                         
                  ),
                  column(6,
                         sliderInput("width_ORAnetwork_microarray_norm", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_ORAnetwork_microarray_norm', 
                                        'Download')
                  )
                )
                
              )
              
            ))
          })
          
          #-------------------------------------------------------------#
          # ORA settings
          #-------------------------------------------------------------#
          observe({
            
            # Print table with settings
            output$ORASettings_microarray_norm <- DT::renderDataTable({
              return(rv$ORA_settings)
            },options = list(pageLength = 10,
                             dom = 't',
                             rownames = FALSE),
            selection = "none")
            
            # Download button
            output$downloadORASettings_microarray_norm <- downloadHandler(
              filename = "ORA_Settings.csv",
              content = function(file){
                write.csv(rv$ORA_settings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
          })
          
          #--------------------------------------------------------------#
          # ORA report
          #--------------------------------------------------------------#
          
          output$ORAreport_microarray_norm <- downloadHandler(
            # For PDF output, change this to "report.pdf"
            filename = "ORAreport.html",
            content = function(file) {
              # Copy the report file to a temporary directory before processing it, in
              # case we don't have write permissions to the current working dir (which
              # can happen when deployed).
              tempReport <- file.path(tempdir(), "ORAreport_microarray_norm.Rmd")
              file.copy("Reports/ORAreport_microarray_norm.Rmd", tempReport, overwrite = TRUE)
              
              # Set up parameters to pass to Rmd document
              params <- list(ORASettings = rv$ORA_settings,
                             ORATable = rv$ORA_data)
              
              # Knit the document, passing in the `params` list, and eval it in a
              # child of the global environment (this isolates the code in the document
              # from the code in this app).
              rmarkdown::render(tempReport, output_file = file,
                                params = params,
                                envir = new.env(parent = globalenv())
              )
            }
          )
          
        } # EO ifelse
      }) # EO observe
      
      
      
      
      
      #********************************************************************#
      # UI: ORA output
      #********************************************************************#
      observe({
        if (!is.null(rv$ORA_data)){
          output$UI_output_ORA_microarray_norm <- renderUI({
            tagList(
              tabsetPanel(
                
                #--------------------------------------------------------------#
                # ORA statistics table
                #--------------------------------------------------------------#
                tabPanel("Statistics",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of statistics table
                         h3(strong("Statistics table")),
                         h5("The ORA statistics table encompasses the output of the gene 
                              overrepresentation analysis."),
                         hr(),
                         
                         # Statistics table
                         DT::dataTableOutput(outputId = "ORA_table_microarray_norm") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         
                         # Download button
                         downloadButton("download_ORA_table_microarray_norm", 
                                        "Download"),
                         br(),
                         
                         # Title + description of gene table
                         htmlOutput("text_ORAgene_table_microarray_norm"),
                         h5(paste0("The gene table encompasses the statistics of all genes 
                              from the selected geneset.")),
                         hr(),
                         
                         # Gene table
                         DT::dataTableOutput(outputId = "ORAgene_table_microarray_norm") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                
                #--------------------------------------------------------------#
                # ORA bar chart
                #--------------------------------------------------------------#
                tabPanel("Bar chart",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of bar chart
                         h3(strong("Bar chart")),
                         h5("The bar chart visualizes the results from the overrepresentation analysis."),
                         hr(),
                         actionButton("download_ORAplot_microarray_norm", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         shinyWidgets::dropdownButton(
                           tags$div(
                             style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                             
                             # Color gradient
                             tags$h4("Color gradient"),
                             selectInput(inputId = "color_ORAplot_microarray_norm",
                                         label = NULL,
                                         choices = c("Viridis", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             br(),
                             
                             # Number of genesets
                             tags$h4("# Genesets"),
                             sliderInput(
                               inputId = "nSets_ORAplot_microarray_norm",
                               label = NULL,
                               value = 10,
                               min = 5,
                               max = 20,
                               step = 1),
                             
                             br(),br()
                           ),
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"),
                           
                           tooltip = shinyWidgets::tooltipOptions(
                             title = "Click to personalize the barchart!")
                           
                         ), # EO dropdownButton
                         
                         # Interactive plot output
                         plotly::plotlyOutput("ORAplot_microarray_norm") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                #--------------------------------------------------------------#
                # ORA network diagram
                #--------------------------------------------------------------#
                tabPanel("Network diagram",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of the network diagram
                         h3(strong("Network diagram")),
                         h5("The network diagram visualize the similarity between the most significant genesets."),
                         hr(),
                         actionButton("download_ORAnetwork_microarray_norm", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         shinyWidgets::dropdownButton(
                           tags$div(
                             style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                             
                             # Color gradient
                             tags$h4("Color gradient"),
                             selectInput(inputId = "color_ORAnetwork_microarray_norm",
                                         label = NULL,
                                         choices = c("Viridis", 
                                                     "Yellow-red", 
                                                     "Blues", 
                                                     "Reds")),
                             br(),
                             
                             # Network layout
                             tags$h4("Network layout"),
                             selectInput(inputId = "layout_ORAnetwork_microarray_norm",
                                         label = NULL,
                                         choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                     'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                         selected = 'graphopt',
                                         multiple = FALSE),
                             br(),
                             
                             # Number of genesets
                             tags$h4("# Genesets"),
                             sliderInput(
                               inputId = "nSets_ORAnetwork_microarray_norm",
                               label = NULL,
                               value = 10,
                               min = 5,
                               max = 20,
                               step = 1),
                             br(),br()
                           ),
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"),
                           
                           tooltip = shinyWidgets::tooltipOptions(
                             title = "Click to personalize the network!")
                           
                         ), # EO dropdownButton
                         
                         
                         
                         # Make plot
                         plotOutput("ORAnetwork_microarray_norm") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         br(),
                         
                ),
                
                #--------------------------------------------------------------#
                # Settings tab
                #--------------------------------------------------------------#
                tabPanel("Settings overview",
                         icon = icon("fas fa-file"),
                         h3(strong("ORA settings")),
                         h5("To enhance reproducibility, view and download the overview of the chosen ORA settings."),
                         hr(),
                         DT::dataTableOutput(outputId = "ORASettings_microarray_norm") %>% 
                           withSpinner(color="#0dc5c1"),
                         br(),
                         downloadButton("downloadORASettings_microarray_norm", 
                                        "Download table"),
                         downloadButton("downloadSessionInfo_microarray_norm", 
                                        "Session info")
                         
                )
                
              ) # EO tabSetPanel
            ) # EO tagList
          }) # EO renderUI
        }# EO if !is.null(rv$ORA_data)
        
        # Allow user to download ORA report
        output$UI_ORAreport_microarray_norm <- renderUI({
          req(rv$ORA_data)
          tagList(
            shinyWidgets::downloadBttn(outputId = "ORAreport_microarray_norm",
                                       label = "Get ORA report",
                                       style = "jelly",
                                       color = "primary",
                                       icon = icon("download")),
            br(),br()
          )
        })
        
      }) # EO observe
      
      
    } # EO ORA
    
    
    #==========================================================================#
    
    # GSEA
    
    #==========================================================================#
    
    if (input$ORA_or_GSEA_microarray_norm == "GSEA"){
      
      # Show modal
      shinybusy::show_modal_spinner(text = "Gene set enrichment analysis...",
                                    color="#0dc5c1")
      
      # Perform GSEA:
      rv$GSEA_data <- performGSEA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                                  geneset = input$geneset_ORA_microarray_norm,
                                  geneID_col = input$geneID_ORA_microarray_norm,
                                  geneID_type = input$selID_ORA_microarray_norm,
                                  organism = input$organism_ORA_microarray_norm,
                                  rankingVar = input$ranking_GSEA_microarray_norm)
      
      if (input$ranking_GSEA_microarray_norm == "pvalue"){
        rankvar <- "-log p-value"
      }
      if (input$ranking_GSEA_microarray_norm == "signed_pvalue"){
        rankvar <- "-log p-value x sign logFC"
      }
      if (input$ranking_GSEA_microarray_norm == "logFC"){
        rankvar <- "logFC"
      }
      
      rv$GSEA_settings <- data.frame(
        Option = c("Comparison",
                   "Geneset",
                   "Gene ID",
                   "Organism",
                   "Method",
                   "Ranking variable"
        ),
        Selected = c(input$comparisons_view_ORA_microarray_norm,
                     input$geneset_ORA_microarray_norm,
                     paste0(input$selID_ORA_microarray_norm, " (top table column: ", input$geneID_ORA_microarray_norm, ")"),
                     input$organism_ORA_microarray_norm,
                     "GSEA",
                     rankvar
        )
      )
      
      
      
      #********************************************************************#
      # Generate GSEA output
      #********************************************************************#
      
      observe({
        
        # Remove modal
        shinybusy::remove_modal_spinner()
        
        # Show error message
        if (is.null(rv$GSEA_data)){
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
            text = "Geneset enrichment analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
            type = "info")
          
          #--------------------------------------------------------------#
          # GSEA statistics table
          #--------------------------------------------------------------#
          
          output$GSEA_table_microarray_norm <- DT::renderDataTable({
            req(input$geneset_ORA_microarray_norm)
            req(rv$GSEA_data)
            output <- rv$GSEA_data@result
            
            # Link to wikipathways website if geneset ID is from WikiPathways
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
            
            # Link to QuickGO website if geneset ID is a GO term
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
          output$download_GSEA_table_microarray_norm <- downloadHandler(
            filename = paste0("GSEATable_",input$comparisons_view_ORA_microarray_norm,"_",input$geneset_ORA_microarray_norm,".csv"),
            content = function(file){
              write.csv(rv$GSEA_data@result, file, quote = FALSE, row.names = FALSE)
            }
          )
          
          # Print statistics of genes in selected Term
          output$GSEAgene_table_microarray_norm <- DT::renderDataTable({
            req(input$GSEA_table_microarray_norm_rows_selected)
            req(rv$GSEA_data)
            
            # Make GSEA gene table
            output <- make_ORAgene_table(ORA_data = rv$GSEA_data,
                                         top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                                         geneID_col = input$geneID_ORA_microarray_norm,
                                         sel_row_ORA = input$GSEA_table_microarray_norm_rows_selected)
            
            output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
            output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
            output$meanExpr <- round(output$meanExpr,3)
            output$log2FC <- round(output$log2FC,3)
            output$`log2FC SE` <- round(output$`log2FC SE`,3)
            
            return(output)
          }, options = list(pageLength = 6), escape = FALSE)
          
          # Text for gene table
          output$text_GSEAgene_table_microarray_norm <- renderText({
            req(rv$GSEA_data)
            text <- paste0("<h3><b>Gene table: ",rv$GSEA_data@result[input$GSEA_table_microarray_norm_rows_selected,"ID"],
                           "</b></h3>")
            return(text)
          })
          
          #--------------------------------------------------------------#
          # GSEA barchart
          #--------------------------------------------------------------#
          
          observe({
            req(input$nSets_GSEAplot_microarray_norm)
            req(rv$GSEA_data)
            rv$GSEAplot <- makeGSEAplot(rv$GSEA_data,
                                        nSets = input$nSets_GSEAplot_microarray_norm,
                                        color = c(input$lowcol_GSEAplot_microarray_norm,
                                                  input$midcol_GSEAplot_microarray_norm,
                                                  input$highcol_GSEAplot_microarray_norm))
            
            output$GSEAplot_microarray_norm <- plotly::renderPlotly(rv$GSEAplot)
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          output$realdownload_GSEAplot_microarray_norm <- downloadHandler(
            filename = function(){ifelse(input$static_GSEAplot_microarray_norm, "GSEA barchart.png", "GSEA barchart.html")},
            content = function(file){
              
              if (input$static_GSEAplot_microarray_norm){
                
                
                # Make MA plot
                p <- makeGSEAplot(rv$GSEA_data,
                                  nSets = input$nSets_GSEAplot_microarray_norm,
                                  color = c(input$lowcol_GSEAplot_microarray_norm,
                                            input$midcol_GSEAplot_microarray_norm,
                                            input$highcol_GSEAplot_microarray_norm),
                                  static = TRUE)
                
                ggplot2::ggsave(plot = p, 
                                filename = file,
                                width = input$width_GSEAplot_microarray_norm,
                                height = input$height_GSEAplot_microarray_norm,
                                units = "px")
              } else{
                htmlwidgets::saveWidget(rv$GSEAplot, 
                                        file)
              }
            }
          )
          
          
          # Make modal
          observeEvent(input$download_GSEAplot_microarray_norm, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(12, align = "left",
                         shinyWidgets::materialSwitch(
                           inputId = "static_GSEAplot_microarray_norm",
                           label = "Click to make static plot",
                           value = FALSE, 
                           status = "primary"))
                ),
                fluidRow(
                  column(6,
                         conditionalPanel(
                           condition = "input.static_GSEAplot_microarray_norm==true",
                           sliderInput("height_GSEAplot_microarray_norm", 
                                       "Height",
                                       min = 800, max = 2000,
                                       value = 1200, step = 10,
                                       width = "100%")
                         )
                  ),
                  column(6,
                         conditionalPanel(
                           condition = "input.static_GSEAplot_microarray_norm==true",
                           sliderInput("width_GSEAplot_microarray_norm", 
                                       "Width",
                                       min = 800, max = 2000,
                                       value = 1500, step = 10,
                                       width = "100%")
                         )
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_GSEAplot_microarray_norm', 
                                        'Download')
                  )
                )
                
              )
              
            ))
          })
          
          #--------------------------------------------------------------#
          # GSEA network diagram
          #--------------------------------------------------------------#
          
          observe({
            req(input$layout_GSEAnetwork_microarray_norm)
            req(input$nSets_GSEAnetwork_microarray_norm)
            req(rv$GSEA_data)
            rv$GSEAnetwork <- makeGSEAnetwork(GSEA_data = rv$GSEA_data,
                                              layout = input$layout_GSEAnetwork_microarray_norm,
                                              nSets = input$nSets_GSEAnetwork_microarray_norm,
                                              color = c(input$lowcol_GSEAnetwork_microarray_norm,
                                                        input$midcol_GSEAnetwork_microarray_norm,
                                                        input$highcol_GSEAnetwork_microarray_norm))
            
            output$GSEAnetwork_microarray_norm <- renderPlot(rv$GSEAnetwork,
                                                            height = 500, 
                                                            width = 800)
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          output$realdownload_GSEAnetwork_microarray_norm <- downloadHandler(
            filename = "GSEA network.png",
            content = function(file){
              
              ggplot2::ggsave(plot = rv$GSEAnetwork, 
                              filename = file,
                              width = input$width_GSEAnetwork_microarray_norm*2,
                              height = input$height_GSEAnetwork_microarray_norm*2,
                              units = "px")
            }
          )
          
          
          # Make modal
          observeEvent(input$download_GSEAnetwork_microarray_norm, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,
                         sliderInput("height_GSEAnetwork_microarray_norm", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                         
                  ),
                  column(6,
                         sliderInput("width_GSEAnetwork_microarray_norm", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_GSEAnetwork_microarray_norm', 
                                        'Download')
                  )
                )
                
              )
              
            ))
          })
          
          #-------------------------------------------------------------#
          # GSEA settings
          #-------------------------------------------------------------#
          observe({
            
            # Print table with settings
            output$GSEASettings_microarray_norm <- DT::renderDataTable({
              return(rv$GSEA_settings)
            },options = list(pageLength = 10,
                             dom = 't',
                             rownames = FALSE),
            selection = "none")
            
            # Download button
            output$downloadGSEASettings_microarray_norm <- downloadHandler(
              filename = "GSEA_Settings.csv",
              content = function(file){
                write.csv(rv$GSEA_settings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
          })
          
          
          #--------------------------------------------------------------#
          # GSEA report
          #--------------------------------------------------------------#
          
          output$GSEAreport_microarray_norm <- downloadHandler(
            # For PDF output, change this to "report.pdf"
            filename = "GSEAreport.html",
            content = function(file) {
              # Copy the report file to a temporary directory before processing it, in
              # case we don't have write permissions to the current working dir (which
              # can happen when deployed).
              tempReport <- file.path(tempdir(), "GSEAreport_microarray_norm.Rmd")
              file.copy("Reports/GSEAreport_microarray_norm.Rmd", tempReport, overwrite = TRUE)
              
              # Set up parameters to pass to Rmd document
              params <- list(GSEASettings = rv$GSEA_settings,
                             GSEATable = rv$GSEA_data)
              
              # Knit the document, passing in the `params` list, and eval it in a
              # child of the global environment (this isolates the code in the document
              # from the code in this app).
              rmarkdown::render(tempReport, output_file = file,
                                params = params,
                                envir = new.env(parent = globalenv())
              )
            }
          )
          
          
        } # EO ifelse
      }) # EO observe
      
      
      
      #********************************************************************#
      # UI: GSEA output
      #********************************************************************#
      observe({
        if (!is.null(rv$GSEA_data)){
          output$UI_output_GSEA_microarray_norm <- renderUI({
            tagList(
              tabsetPanel(
                
                #--------------------------------------------------------------#
                # GSEA statistics table
                #--------------------------------------------------------------#
                tabPanel("Statistics",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of statistics table
                         h3(strong("Statistics table")),
                         h5("The statistics table encompasses the output of the gene set enrichment analysis."),
                         hr(),
                         
                         # Statistics table
                         DT::dataTableOutput(outputId = "GSEA_table_microarray_norm") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         
                         # Download button
                         downloadButton("download_GSEA_table_microarray_norm",
                                        "Download"),
                         br(),
                         
                         # Title + description of gene table
                         htmlOutput("text_GSEAgene_table_microarray_norm"),
                         h5(paste0("The gene table encompasses the statistics of all genes
                              from the selected geneset.")),
                         hr(),
                         
                         # Gene table
                         DT::dataTableOutput(outputId = "GSEAgene_table_microarray_norm") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                
                #--------------------------------------------------------------#
                # GSEA bar chart
                #--------------------------------------------------------------#
                tabPanel("Bar chart",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of bar chart
                         h3(strong("Bar chart")),
                         h5("The bar chart visualizes the results from the gene set enrichment analysis."),
                         hr(),
                         actionButton("download_GSEAplot_microarray_norm", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         shinyWidgets::dropdownButton(
                           tags$div(
                             style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                             
                             # Color gradient
                             tags$h4("Color gradient"),
                             colourpicker::colourInput("lowcol_GSEAplot_microarray_norm",
                                                       NULL,
                                                       "#000072"),
                             colourpicker::colourInput("midcol_GSEAplot_microarray_norm",
                                                       NULL,
                                                       "#FEE6CE"),
                             colourpicker::colourInput("highcol_GSEAplot_microarray_norm",
                                                       NULL,
                                                       "red"),
                             br(),
                             
                             # Number of genesets
                             tags$h4("# Genesets"),
                             sliderInput(
                               inputId = "nSets_GSEAplot_microarray_norm",
                               label = NULL,
                               value = 10,
                               min = 5,
                               max = 20,
                               step = 1),
                             br(),br()
                           ),
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"),
                           
                           tooltip = shinyWidgets::tooltipOptions(
                             title = "Click to personalize the barchart!")
                           
                         ), # EO dropdownButton
                         
                         # Interactive plot output
                         plotly::plotlyOutput("GSEAplot_microarray_norm") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                #--------------------------------------------------------------#
                # GSEA network diagram
                #--------------------------------------------------------------#
                tabPanel("Network diagram",
                         icon = icon("fas fa-mouse-pointer"),
                         br(),
                         
                         # Title + description of the network diagram
                         h3(strong("Network diagram")),
                         h5("The network diagram visualize the similarity between the most significant gene sets."),
                         hr(),
                         actionButton("download_GSEAnetwork_microarray_norm", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         shinyWidgets::dropdownButton(
                           tags$div(
                             style = "max-height: 300px; overflow-y: auto; padding: 5px;",
                             
                             # Color gradient
                             tags$h4("Color gradient"),
                             colourpicker::colourInput("lowcol_GSEAnetwork_microarray_norm",
                                                       NULL,
                                                       "#000072"),
                             colourpicker::colourInput("midcol_GSEAnetwork_microarray_norm",
                                                       NULL,
                                                       "#FEE6CE"),
                             colourpicker::colourInput("highcol_GSEAnetwork_microarray_norm",
                                                       NULL,
                                                       "red"),
                             br(),
                             
                             # Network layout
                             tags$h4("Network layout"),
                             selectInput(inputId = "layout_GSEAnetwork_microarray_norm",
                                         label = NULL,
                                         choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds',
                                                     'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                         selected = 'graphopt',
                                         multiple = FALSE),
                             br(),
                             
                             # Number of genesets
                             tags$h4("# Gene sets"),
                             sliderInput(
                               inputId = "nSets_GSEAnetwork_microarray_norm",
                               label = NULL,
                               value = 10,
                               min = 5,
                               max = 20,
                               step = 1),
                             br(),br()
                           ),
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"),
                           
                           tooltip = shinyWidgets::tooltipOptions(
                             title = "Click to personalize the network!")
                           
                         ), # EO dropdownButton
                         
                         # Make plot
                         plotOutput("GSEAnetwork_microarray_norm") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                #--------------------------------------------------------------#
                # Settings tab
                #--------------------------------------------------------------#
                tabPanel("Settings overview",
                         icon = icon("fas fa-file"),
                         h3(strong("GSEA settings")),
                         h5("To enhance reproducibility, view and download the overview of the chosen GSEA settings."),
                         hr(),
                         DT::dataTableOutput(outputId = "GSEASettings_microarray_norm") %>% 
                           withSpinner(color="#0dc5c1"),
                         br(),
                         downloadButton("downloadGSEASettings_microarray_norm", 
                                        "Download table"),
                         downloadButton("downloadSessionInfo_microarray_norm", 
                                        "Session info")
                         
                )
                
              ) # EO tabSetPanel
            ) # EO tagList
          }) # EO renderUI
        }# EO if !is.null(rv$GSEA_data)
        
        
        # Allow user to download GSEA report
        output$UI_GSEAreport_microarray_norm <- renderUI({
          req(rv$GSEA_data)
          tagList(
            shinyWidgets::downloadBttn(outputId = "GSEAreport_microarray_norm",
                                       label = "Get GSEA report",
                                       style = "jelly",
                                       color = "primary",
                                       icon = icon("download")),
          )
        })
        
      }) # EO observe
      
      
    } # EO GSEA
    
  }) # EO observeEvent
  
})