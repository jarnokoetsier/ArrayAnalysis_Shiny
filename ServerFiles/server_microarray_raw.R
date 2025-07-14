observe({
  # Make list for reactive values
  rv <- reactiveValues()
  
  # Download session info
  output$downloadSessionInfo_microarray_raw <- downloadHandler(
    filename = "sessionInfo.txt",
    content = function(file){
      writeLines(capture.output(sessionInfo()), file)
    }
  )
  
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
    
    # Example metadata file
    output$downloadmeta_example_microarray_raw <- downloadHandler(
      filename = "MetaData_example.csv",
      content = function(file){
        file.copy("Data/Microarray/metaData_GSE6955.csv", file)
      }
    )
    
    # Example expression matrix file
    output$downloadexpr_example_microarray_raw <- downloadHandler(
      filename = "ExprData_example.zip",
      content = function(file){
        file.copy("Data/Microarray/GSE6955_RAW.zip", file)
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
      
      # Check whether all expression samples have metadata available
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
      
      # Allow the user to adjust metadata
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
                     h5("The table shows the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                     hr(),
                     DT::dataTableOutput(outputId = "exprTable_upload_microarray_raw") %>% 
                       withSpinner(color="#0dc5c1")),
            tabPanel("Metadata",                  # Meta table
                     h3(strong("Metadata")),
                     h5("Please check if the metadata has been correctly imported."),
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
      # No common samples between metadata and expression data
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
      
      # Check whether all expression samples have metadata available
      if (nrow(rv$metaData) < length(rv$celfiles)){
        shinyWidgets::sendSweetAlert(
          session = session,
          title = "Warning!",
          text = "One or more sample(s) in the expression data file do(es) not have 
          metadata available. These samples are excluded from the analysis.",
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
      
      # Allow the user to adjust the metadata
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
                     h5("The table shows the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                     hr(),
                     DT::dataTableOutput(outputId = "exprTable_upload_microarray_raw") %>% 
                       withSpinner(color="#0dc5c1")),
            
            # Tab with the metadata
            tabPanel("Metadata",                  
                     h3(strong("Metadata")),
                     h5("Please check if the metadata has been correctly imported."),
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
        
        # Show success message (only if the metadata is available for 
        # each expression file)
        if (nrow(rv$metaData) >= length(rv$celfiles)){
          
          # Show information tab
          shinyWidgets::sendSweetAlert(
            session = session,
            title = "Info",
            text = "The data has been succesfully uploaded. Please check 
                  the tables on this page to see whether the data has been correctly uploaded.",
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
      
      # No common samples between metadata and expression data
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
      selectInput(inputId = "select_outliers_microarray_raw",
                  label = tags$span(
                    "Select samples to be removed", 
                    tags$span(
                      icon(
                        name = "question-circle",
                      ) 
                    ) |>
                      prompter::add_prompt(message = "Select one or more sample(s) to exclude from the analysis.", 
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
    selectInput(inputId = "groupselect_microarray_raw",
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
    # Generate output: QC tables/plots
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
      
      if (length(levels(rv$experimentFactor)) > 6){
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
    observe({
      req(input$geneboxplot_file_microarray_raw)
      output$realdownload_geneboxplot_microarray_raw <- downloadHandler(
        filename = ifelse(input$geneboxplot_file_microarray_raw == "PNG", "GeneBoxplot.png",
                          ifelse(input$geneboxplot_file_microarray_raw == "PDF", "GeneBoxplot.pdf",
                                 "GeneBoxplot.tif")),
        content = function(file){
          ggplot2::ggsave(plot = rv$temp1, 
                          filename = file,
                          width = input$width_geneboxplot_microarray_raw,
                          height = input$height_geneboxplot_microarray_raw,
                          units = "px")
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_geneboxplot_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "geneboxplot_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF"),
                     selected = "PNG",
                     inline = TRUE,
                     fill = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   sliderInput("height_geneboxplot_microarray_raw", 
                               "Height",
                               min = 800, max = 5000,
                               value = 2100, step = 10,
                               width = "100%"),
            ),
            column(6,
                   sliderInput("width_geneboxplot_microarray_raw", 
                               "Width",
                               min = 800, max = 5000,
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
      
      if (length(levels(rv$experimentFactor)) > 5){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      } else{
        legendColors <- c(input$boxplots_col1_microarray_raw,
                          input$boxplots_col2_microarray_raw,
                          input$boxplots_col3_microarray_raw,
                          input$boxplots_col4_microarray_raw,
                          input$boxplots_col5_microarray_raw)
      }
      if (length(legendColors) != length(levels(rv$experimentFactor))){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      }
      names(legendColors) <- levels(rv$experimentFactor)
      
      getBoxplots(experimentFactor = rv$experimentFactor,
                  legendColors = legendColors,
                  normData = rv$normData,
                  RNASeq = FALSE,
                  width = session$clientData$output_boxplots_microarray_raw_width,
                  height = session$clientData$output_boxplots_microarray_raw_height)
    }, deleteFile = TRUE)
    
    # Allow users to set colors
    observe({
      test <- length(levels(rv$experimentFactor))
      
      if (test == 2){
        output$UI_color_boxplots_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("boxplots_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("boxplots_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2])
            )
          )
        })
      }
      if (test == 3){
        output$UI_color_boxplots_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("boxplots_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("boxplots_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("boxplots_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3])
            )
          )
        })
      }
      if (test == 4){
        output$UI_color_boxplots_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("boxplots_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("boxplots_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("boxplots_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3]),
                     colourpicker::colourInput("boxplots_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[4])
            )
          )
        })
        
      }
      if (test == 5){
        output$UI_color_boxplots_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("boxplots_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("boxplots_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("boxplots_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3]),
                     colourpicker::colourInput("boxplots_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[4]),
                     colourpicker::colourInput("boxplots_col5_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[5])
            )
          )
        })
        
      }
      if (test > 5){
        output$UI_color_boxplots_microarray_raw <- renderUI({
          tags$div(class = "dropdown-content",
                   h5("There are too many experimental groups to select colour manually")
          )
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$boxplots_file_microarray_raw)
      output$realdownload_boxplots_microarray_raw <- downloadHandler(
        filename = ifelse(input$boxplots_file_microarray_raw == "PNG", "QC_Boxplots.png",
                          ifelse(input$boxplots_file_microarray_raw == "PDF", "QC_Boxplots.pdf",
                                 "QC_Boxplots.tif")),
        content = function(file){
          
          if (input$boxplots_file_microarray_raw == "PNG"){
            png(file,
                width=input$width_boxplots_microarray_raw*3,
                height=input$height_boxplots_microarray_raw*3,
                pointsize=24,
                res = 300)
          }
          
          if (input$boxplots_file_microarray_raw == "PDF"){
            pdf(file,
                width=input$width_boxplots_microarray_raw/32,
                height=input$height_boxplots_microarray_raw/32,
                pointsize=24*3)
          }
          if (input$boxplots_file_microarray_raw == "TIF"){
            tiff(file,
                 width=input$width_boxplots_microarray_raw*3,
                 height=input$height_boxplots_microarray_raw*3,
                 pointsize=24,
                 res = 300)
          }
          
          
          if (length(levels(rv$experimentFactor)) > 5){
            legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
          } else{
            legendColors <- c(input$boxplots_col1_microarray_raw,
                              input$boxplots_col2_microarray_raw,
                              input$boxplots_col3_microarray_raw,
                              input$boxplots_col4_microarray_raw,
                              input$boxplots_col5_microarray_raw)
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
    })
    
    
    # Make modal
    observeEvent(input$download_boxplots_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "boxplots_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
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
    
    # Density plots of all genes together
    output$densityplots_microarray_raw <- plotly::renderPlotly({
      
      if (length(levels(rv$experimentFactor)) > 5){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      } else{
        legendColors <- c(input$density_col1_microarray_raw,
                          input$density_col2_microarray_raw,
                          input$density_col3_microarray_raw,
                          input$density_col4_microarray_raw,
                          input$density_col5_microarray_raw)
      }
      if (length(legendColors) != length(levels(rv$experimentFactor))){
        legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
      }
      names(legendColors) <- levels(rv$experimentFactor)
      
      rv$density <- getDensityplots(experimentFactor = rv$experimentFactor,
                                    legendColors = legendColors,
                                    normMatrix = rv$normMatrix,
                                    RNASeq = FALSE)
      
    })
    
    # Allow users to set colors
    observe({
      test <- length(levels(rv$experimentFactor))
      
      if (test == 2){
        output$UI_color_density_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("density_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("density_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2])
            )
          )
        })
      }
      if (test == 3){
        output$UI_color_density_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("density_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("density_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("density_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3])
            )
          )
        })
      }
      if (test == 4){
        output$UI_color_density_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("density_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("density_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("density_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3]),
                     colourpicker::colourInput("density_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[4])
            )
          )
        })
        
      }
      if (test == 5){
        output$UI_color_density_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("density_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[1]),
                     colourpicker::colourInput("density_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[2]),
                     colourpicker::colourInput("density_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[3]),
                     colourpicker::colourInput("density_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[4]),
                     colourpicker::colourInput("density_col5_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(rv$experimentFactor)$legendColors[5])
            )
          )
        })
        
      }
      if (test > 5){
        output$UI_color_density_microarray_raw <- renderUI({
          tags$div(class = "dropdown-content",
                   h5("There are too many experimental groups to select colour manually")
          )
        })
        
      }
      
    })
    
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$density_file_microarray_raw)
      output$realdownload_density_microarray_raw <- downloadHandler(
        filename = ifelse(input$density_file_microarray_raw == "HTML", "QC_Density.html",
                          ifelse(input$density_file_microarray_raw == "PNG", "QC_Density.png",
                                 ifelse(input$density_file_microarray_raw == "PDF", "QC_Density.pdf",
                                        "QC_Density.tif"))),
        content = function(file){
          
          if (input$density_file_microarray_raw != "HTML"){
            
            if (length(levels(rv$experimentFactor)) > 5){
              legendColors <- colorsByFactor(rv$experimentFactor)$legendColors
            } else{
              legendColors <- c(input$density_col1_microarray_raw,
                                input$density_col2_microarray_raw,
                                input$density_col3_microarray_raw,
                                input$density_col4_microarray_raw,
                                input$density_col5_microarray_raw)
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
                            width = input$width_density_microarray_raw,
                            height = input$height_density_microarray_raw,
                            units = "px")
          } else{
            htmlwidgets::saveWidget(rv$density, 
                                    file)
          }
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_density_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "density_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF", "HTML"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.density_file_microarray_raw!=`HTML`",
                     sliderInput("height_density_microarray_raw", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.density_file_microarray_raw!=`HTML`",
                     sliderInput("width_density_microarray_raw", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_density_microarray_raw', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    
    #********************************************************************#
    # Output 4: Sample-sample correlation heatmap
    #********************************************************************#
    
    # Heatmap of sample-sample correlations
    output$heatmap_microarray_raw  <- plotly::renderPlotly({
      
      # Make color factor
      if(length(input$colorFactor_heatmap_microarray_raw) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw])
      }
      
      # Set colors
      if (length(levels(colorFactor)) > 5){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      } else{
        legendColors <- c(input$heatmap_col1_microarray_raw,
                          input$heatmap_col2_microarray_raw,
                          input$heatmap_col3_microarray_raw,
                          input$heatmap_col4_microarray_raw,
                          input$heatmap_col5_microarray_raw)
      }
      if (length(legendColors) != length(levels(colorFactor))){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      }
      names(legendColors) <- levels(colorFactor)
      
      # Make heatmap
      rv$heatmap <- getHeatmap(experimentFactor = colorFactor,
                               legendColors = legendColors,
                               normMatrix = rv$normMatrix,
                               clusterOption1 = input$clusteroption1_microarray_raw,
                               clusterOption2 = input$clusteroption2_microarray_raw,
                               theme = input$heatmaptheme_microarray_raw)
      return(rv$heatmap)
    })
    
    # Allow users to set colors
    observe({
      req(input$colorFactor_heatmap_microarray_raw)
      if(length(input$colorFactor_heatmap_microarray_raw) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw])
      }
      test <- length(levels(colorFactor))
      
      if (test == 2){
        output$UI_color_heatmap_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Heatmap theme"),
                     selectInput(inputId = 'heatmaptheme_microarray_raw',
                                 label = NULL,
                                 choices = c("Default", 
                                             "Yellow-red", 
                                             "Blues", 
                                             "Reds")),
                     br(),
                     h4("Side color"),
                     colourpicker::colourInput("heatmap_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("heatmap_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2])
            )
          )
        })
      }
      
      if (test == 3){
        output$UI_color_heatmap_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Heatmap theme"),
                     selectInput(inputId = 'heatmaptheme_microarray_raw',
                                 label = NULL,
                                 choices = c("Default", 
                                             "Yellow-red", 
                                             "Blues", 
                                             "Reds")),
                     br(),
                     h4("Side color"),
                     colourpicker::colourInput("heatmap_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("heatmap_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("heatmap_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3])
            )
          )
        })
      }
      if (test == 4){
        output$UI_color_heatmap_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Heatmap theme"),
                     selectInput(inputId = 'heatmaptheme_microarray_raw',
                                 label = NULL,
                                 choices = c("Default", 
                                             "Yellow-red", 
                                             "Blues", 
                                             "Reds")),
                     br(),
                     h4("Side color"),
                     colourpicker::colourInput("heatmap_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("heatmap_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("heatmap_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3]),
                     colourpicker::colourInput("heatmap_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[4])
            )
          )
        })
        
      }
      if (test == 5){
        output$UI_color_heatmap_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Heatmap theme"),
                     selectInput(inputId = 'heatmaptheme_microarray_raw',
                                 label = NULL,
                                 choices = c("Default", 
                                             "Yellow-red", 
                                             "Blues", 
                                             "Reds")),
                     br(),
                     h4("Side color"),
                     colourpicker::colourInput("heatmap_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("heatmap_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("heatmap_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3]),
                     colourpicker::colourInput("heatmap_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[4]),
                     colourpicker::colourInput("heatmap_col5_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[5])
            )
          )
        })
        
      }
      if (test > 5){
        output$UI_color_heatmap_microarray_raw <- renderUI({
          tags$div(class = "dropdown-content",
                   h4("Heatmap theme"),
                   selectInput(inputId = 'heatmaptheme_microarray_raw',
                               label = NULL,
                               choices = c("Default", 
                                           "Yellow-red", 
                                           "Blues", 
                                           "Reds")),
                   br()
          )
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$heatmap_file_microarray_raw)
      output$realdownload_heatmap_microarray_raw <- downloadHandler(
        filename = ifelse(input$heatmap_file_microarray_raw == "HTML", "QC_Heatmap.html",
                          ifelse(input$heatmap_file_microarray_raw == "PNG", "QC_Heatmap.png",
                                 ifelse(input$heatmap_file_microarray_raw == "PDF", "QC_Heatmap.pdf",
                                        "QC_Heatmap.tif"))),
        content = function(file){
          
          if (input$heatmap_file_microarray_raw != "HTML"){
            
            # Make color factor
            if(length(input$colorFactor_heatmap_microarray_raw) > 1){
              colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw], 1, paste, collapse = "_" ))
            } else{
              colorFactor <- factor(rv$metaData_fil[,input$colorFactor_heatmap_microarray_raw])
            }
            
            # Set colors
            if (length(levels(colorFactor)) > 5){
              legendColors <- colorsByFactor(colorFactor)$legendColors
            } else{
              legendColors <- c(input$heatmap_col1_microarray_raw,
                                input$heatmap_col2_microarray_raw,
                                input$heatmap_col3_microarray_raw,
                                input$heatmap_col4_microarray_raw,
                                input$heatmap_col5_microarray_raw)
            }
            if (length(legendColors) != length(levels(colorFactor))){
              legendColors <- colorsByFactor(colorFactor)$legendColors
            }
            names(legendColors) <- levels(colorFactor)
            
            # Make heatmap
            
            getHeatmap_static(experimentFactor = colorFactor,
                              legendColors = legendColors,
                              normMatrix = rv$normMatrix,
                              clusterOption1 = input$clusteroption1_microarray_raw,
                              clusterOption2 = input$clusteroption2_microarray_raw,
                              theme = input$heatmaptheme_microarray_raw,
                              width = input$width_heatmap_microarray_raw,
                              height = input$height_heatmap_microarray_raw,
                              filetype = input$heatmap_file_microarray_raw,
                              file)
          } else{
            htmlwidgets::saveWidget(rv$heatmap, 
                                    file)
          }
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_heatmap_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "heatmap_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF", "HTML"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.heatmap_file_microarray_raw!=`HTML`",
                     sliderInput("height_heatmap_microarray_raw", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.heatmap_file_microarray_raw!=`HTML`",
                     sliderInput("width_heatmap_microarray_raw", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_heatmap_microarray_raw', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    
    #********************************************************************#
    # Output 5: PCA score plot
    #********************************************************************#
    
    #Perform PCA
    rv$PCA_data <- prcomp(t(rv$normMatrix[apply(rv$normMatrix, 1, var) != 0,]),
                          retx = TRUE, 
                          center = TRUE,
                          scale.= TRUE)
    
    
    # Make PCA plot
    output$PCA_microarray_raw <- plotly::renderPlotly({
      
      # Get factor to color by
      if(length(input$colorFactor_PCA_microarray_raw) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw])
      }
      
      # Set colors
      if (length(levels(colorFactor)) > 5){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      } else{
        legendColors <- c(input$PCA_col1_microarray_raw,
                          input$PCA_col2_microarray_raw,
                          input$PCA_col3_microarray_raw,
                          input$PCA_col4_microarray_raw,
                          input$PCA_col5_microarray_raw)
      }
      if (length(legendColors) != length(levels(colorFactor))){
        legendColors <- colorsByFactor(colorFactor)$legendColors
      }
      
      # Make PCA score plot
      rv$PCAplot <- plot_PCA(PC_data = rv$PCA_data, 
                             colorFactor = colorFactor,
                             legendColors = legendColors, 
                             xpc = as.numeric(stringr::str_remove(input$xpca_microarray_raw,"PC")), 
                             ypc = as.numeric(stringr::str_remove(input$ypca_microarray_raw,"PC")), 
                             zpc = ifelse(input$xyz_microarray_raw,as.numeric(stringr::str_remove(input$zpca_microarray_raw,"PC")),3), 
                             xyz = input$xyz_microarray_raw)
      return(rv$PCAplot)
    })
    
    
    # Allow users to set colors
    observe({
      req(input$colorFactor_PCA_microarray_raw)
      if(length(input$colorFactor_PCA_microarray_raw) > 1){
        colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw], 1, paste, collapse = "_" ))
      } else{
        colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw])
      }
      test <- length(levels(colorFactor))
      
      if (test == 2){
        output$UI_color_PCA_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("PCA_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("PCA_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2])
            )
          )
        })
      }
      
      if (test == 3){
        output$UI_color_PCA_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("PCA_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("PCA_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("PCA_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3])
            )
          )
        })
      }
      if (test == 4){
        output$UI_color_PCA_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("PCA_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("PCA_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("PCA_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3]),
                     colourpicker::colourInput("PCA_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[4])
            )
          )
        })
        
      }
      if (test == 5){
        output$UI_color_PCA_microarray_raw <- renderUI({
          tagList(
            tags$div(class = "dropdown-content",
                     h4("Colors"),
                     colourpicker::colourInput("PCA_col1_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[1]),
                     colourpicker::colourInput("PCA_col2_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[2]),
                     colourpicker::colourInput("PCA_col3_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[3]),
                     colourpicker::colourInput("PCA_col4_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[4]),
                     colourpicker::colourInput("PCA_col5_microarray_raw", 
                                               NULL, 
                                               colorsByFactor(colorFactor)$legendColors[5])
            )
          )
        })
        
      }
      if (test > 5){
        output$UI_color_PCA_microarray_raw <- renderUI({
          tags$div(class = "dropdown-content",
                   NULL
          )
        })
        
      }
      
    })
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$pca_file_microarray_raw)
      output$realdownload_pca_microarray_raw <- downloadHandler(
        filename = ifelse(input$pca_file_microarray_raw == "HTML", "QC_PCA.html",
                          ifelse(input$pca_file_microarray_raw == "PNG", "QC_PCA.png",
                                 ifelse(input$pca_file_microarray_raw == "PDF", "QC_PCA.pdf",
                                        "QC_PCA.tif"))),
        content = function(file){
          
          if (input$pca_file_microarray_raw != "HTML"){
            
            # Get factor to color by
            if(length(input$colorFactor_PCA_microarray_raw) > 1){
              colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw], 1, paste, collapse = "_" ))
            } else{
              colorFactor <- factor(rv$metaData_fil[,input$colorFactor_PCA_microarray_raw])
            }
            
            # Set colors
            if (length(levels(colorFactor)) > 5){
              legendColors <- colorsByFactor(colorFactor)$legendColors
            } else{
              legendColors <- c(input$PCA_col1_microarray_raw,
                                input$PCA_col2_microarray_raw,
                                input$PCA_col3_microarray_raw,
                                input$PCA_col4_microarray_raw,
                                input$PCA_col5_microarray_raw)
            }
            if (length(legendColors) != length(levels(colorFactor))){
              legendColors <- colorsByFactor(colorFactor)$legendColors
            }
            
            # Make PCA score plot
            p <- plot_PCA_static(PC_data = rv$PCA_data, 
                                 colorFactor = colorFactor,
                                 legendColors = legendColors, 
                                 xpc = as.numeric(stringr::str_remove(input$xpca_microarray_raw,"PC")), 
                                 ypc = as.numeric(stringr::str_remove(input$ypca_microarray_raw,"PC")))
            
            ggplot2::ggsave(plot = p, 
                            filename = file,
                            width = input$width_pca_microarray_raw,
                            height = input$height_pca_microarray_raw,
                            units = "px")
          } else{
            htmlwidgets::saveWidget(rv$PCAplot, 
                                    file)
          }
        }
      )
    })
    
    
    
    # Make modal
    observeEvent(input$download_pca_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "pca_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF", "HTML"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.pca_file_microarray_raw!=`HTML`",
                     sliderInput("height_pca_microarray_raw", 
                                 "Height",
                                 min = 800, max = 2000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.pca_file_microarray_raw!=`HTML`",
                     sliderInput("width_pca_microarray_raw", 
                                 "Width",
                                 min = 800, max = 2000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_pca_microarray_raw', 
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
    
    #********************************************************************#
    # QC report
    #********************************************************************#
    output$QCreport_microarray_raw <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "QCreport.html",
      content = function(file) {
        shinybusy::show_modal_spinner(text = "Making QC report...",
                                      color="#0dc5c1")
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "QCreport_microarray_raw.Rmd")
        file.copy("Reports/QCreport_microarray_raw.Rmd", tempReport, overwrite = TRUE)
        
        tempLogo <- file.path(tempdir(), "logo_main.PNG")
        file.copy("www/logo_main.PNG", tempLogo, overwrite = TRUE)
        
        tempHeader <- file.path(tempdir(), "header.html")
        file.copy("www/header.html", tempHeader, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(processingSettings = rv$processingSettings,
                       experimentFactor = rv$experimentFactor,
                       legendColors = colorsByFactor(rv$experimentFactor)$legendColors,
                       normData = rv$normData,
                       normMatrix = rv$normMatrix,
                       PCAData = rv$PCA_data,
                       dir = tempdir(),
                       ArrayAnalysis_version = ArrayAnalysis_version)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
        shinybusy::remove_modal_spinner()
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
                   h5("Here you can view the normalized and log-transformed intensities. 
                              Click on the table to explore the data!"),
                   hr(),
                   DT::dataTableOutput(outputId = "exprTable_microarray_raw") %>% 
                     withSpinner(color="#0dc5c1"),
                   downloadButton("downloadNormalizedData_microarray_raw", 
                                  "Download"),
                   br(),
                   br(),
                   
                   # Dropdown Button to adjust the plot settings
                   tags$div(class = "dropdown",
                            tags$button(class = "dropbtn", icon("cog")),
                            tags$div(class = "dropdown-content",
                                     
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
                                       tags$h4("Click to change boxplot colors"),
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
                            )), # EO dropdownButton
                   
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
                   br(),
                   br(),
                   tags$div(class = "dropdown",
                            tags$button(class = "dropbtn", icon("cog")),
                            uiOutput("UI_color_boxplots_microarray_raw")
                            
                   ),
                   h2(strong("Boxplot of normalized intensities"), align = "center"),
                   h4("Distributions should be comparable between samples", align = "center"),
                   plotOutput(outputId = "boxplots_microarray_raw",
                              width = "65vw", height = "80vw")%>% 
                     shinycssloaders::withSpinner(color="#0dc5c1")
          ),
          
          # Density plot
          tabPanel("Density plots",
                   icon = icon("fas fa-mouse-pointer"),
                   br(),
                   actionButton('download_density_microarray_raw', 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   # customize heatmap
                   tags$div(class = "dropdown",
                            tags$button(class = "dropbtn", icon("cog")),
                            uiOutput("UI_color_density_microarray_raw")
                   ),
                   h2(strong("Density plot of normalized intensities"), align = "center"),
                   h4("Distributions should be comparable between samples", align = "center"),
                   plotly::plotlyOutput(outputId = "densityplots_microarray_raw",
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
                            selectInput(
                              inputId = "colorFactor_heatmap_microarray_raw",
                              label = "Side color",
                              choices = colnames(rv$metaData_fil),
                              selected = rv$experimentName,
                              multiple = TRUE)
                     ),
                     # Select distance method
                     column(3,
                            selectInput(
                              inputId = "clusteroption1_microarray_raw",
                              label = "Distance",
                              choices = c("Pearson","Spearman",
                                          "Euclidean"))
                     ),
                     
                     #Clustering method
                     column(3,
                            selectInput(
                              inputId = "clusteroption2_microarray_raw",
                              label = "Clustering",
                              choices = c("ward.D2","single",
                                          "complete","average",
                                          "mcquitty","median",
                                          "centroid"))
                     )
                   ),
                   hr(),
                   actionButton('download_heatmap_microarray_raw', 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   
                   # customize heatmap
                   tags$div(class = "dropdown",
                            tags$button(class = "dropbtn", icon("cog")),
                            uiOutput("UI_color_heatmap_microarray_raw")
                   ),
                   
                   # Make heatmap
                   plotly::plotlyOutput("heatmap_microarray_raw", 
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
                     inputId = "xyz_microarray_raw",
                     label = "3D",
                     value = FALSE, 
                     status = "info"),
                   
                   hr(),
                   
                   # Set color + 3D/2D
                   fluidRow(
                     column(3,
                            # Color by which factor?
                            selectInput(inputId = "colorFactor_PCA_microarray_raw",
                                        label = "Color by:",
                                        choices = colnames(rv$metaData_fil),
                                        selected = rv$experimentName,
                                        multiple = TRUE)
                     ),
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
                   
                   hr(),
                   
                   actionButton("download_pca_microarray_raw", 
                                "Download figure",
                                icon = shiny::icon("download")),
                   br(),
                   br(),
                   
                   fluidRow(
                     
                     # Customize PCA plot
                     tags$div(class = "dropdown",
                              tags$button(class = "dropbtn", icon("cog")),
                              uiOutput("UI_color_PCA_microarray_raw"),
                     ),
                     
                     # Make PCA plot
                     plotly::plotlyOutput("PCA_microarray_raw",
                                          width = "55vw", height = "30vw")%>% 
                       shinycssloaders::withSpinner(color="#0dc5c1"),
                   ),
                   
                   
          ), # EO PCA tabpanel
          
          # Settings table
          tabPanel("Settings overview",
                   icon = icon("fas fa-file"),
                   h3(strong("Pre-processing settings")),
                   h5("To enhance reproducibility, download the overview of chosen pre-processing settings."),
                   hr(),
                   DT::dataTableOutput(outputId = "processingSettings_microarray_raw") %>% 
                     withSpinner(color="#0dc5c1"),
                   br(),
                   downloadButton("downloadProcessingSettings_microarray_raw", 
                                  "Download table"),
                   downloadButton("downloadSessionInfo_microarray_raw", 
                                  "Session info")
                   
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
        shinyWidgets::downloadBttn(outputId = "QCreport_microarray_raw",
                                   label = "Get report",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("download")),
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
      selectInput(inputId = "covGroups_num_microarray_raw",
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
      selectInput(inputId = "covGroups_char_microarray_raw",
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
      # shinyWidgets::multiInput(
      #   inputId = "comparisons_microarray_raw",
      #   label = "Comparisons:", 
      #   choices = makeComparisons(make.names(unique(rv$experimentFactor))),
      #   selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1]
      # )
      selectInput(
        inputId = "comparisons_microarray_raw",
        label = "Comparisons:", 
        choices = makeComparisons(make.names(unique(rv$experimentFactor))),
        selected = makeComparisons(make.names(unique(rv$experimentFactor)))[1],
        multiple = TRUE
      )
    )
  })
  
  
  
  # Select comparisons
  output$UI_biomart_dataset_microarray_raw <- renderUI({
    req(input$addAnnotation_microarray_raw)
    selectInput(inputId = "biomart_dataset_microarray_raw",
                label = tags$span(
                  "Organism", 
                  tags$span(
                    icon(
                      name = "question-circle",
                    ) 
                  ) |>
                    prompter::add_prompt(message = "Select the organism. 
                                               This information is needed to match the probset IDs to the gene IDs.", 
                                         position = "right",
                                         size = "large")
                ),
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
      
      selectInput(inputId = "biomart_filter_microarray_raw",
                  label = tags$span(
                    "Probeset ID", 
                    tags$span(
                      icon(
                        name = "question-circle",
                      ) 
                    ) |>
                      prompter::add_prompt(message = "Select which probeset ID is 
                                               used in the expression matrix.", 
                                           position = "right",
                                           size = "large")
                  ),
                  choices = filterList[[input$biomart_dataset_microarray_raw]],
                  selected = selFilter(rv$ProbeAnnotation),
                  multiple = FALSE),
      
      selectInput(inputId = "biomart_attributes_microarray_raw",
                  label = tags$span(
                    "Output", 
                    tags$span(
                      icon(
                        name = "question-circle",
                      ) 
                    ) |>
                      prompter::add_prompt(message = "Select which gene IDs should be added 
                                               to the output.", 
                                           position = "right",
                                           size = "large")
                  ),
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
      
      # Make top table
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
      
      # Make settings table
      rv$statSettings <- list()
      for (c in 1:length(input$comparisons_microarray_raw)){
        rv$statSettings[[c]] <- data.frame(
          Option = c("Selected comparison",
                     "Continuous covariate(s)",
                     "Discrete covariate(s)",
                     "Gene annotation dataset",
                     "Gene annotation attribute(s)",
                     "Gene annotation filter"),
          Selected = c(input$comparisons_microarray_raw[c],
                       ifelse(is.null(input$covGroups_num_microarray_raw), " ",
                              paste(input$covGroups_num_microarray_raw, collapse = "; ")),
                       ifelse(is.null(input$covGroups_char_microarray_raw), " ",
                              paste(input$covGroups_char_microarray_raw, collapse = "; ")),
                       ifelse(is.null(rv$top_table_list[[3]]), "N/A", rv$top_table_list[[3]]),
                       paste(input$biomart_attributes_microarray_raw, collapse = "; "),
                       input$biomart_filter_microarray_raw
          )
        )
      }
      names(rv$statSettings) <- input$comparisons_microarray_raw
      
    } else{
      
      # Make top table
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
      
      # Make settings table
      for (c in 1:length(input$comparisons_microarray_raw)){
        rv$statSettings[[c]] <- data.frame(
          Option = c("Comparison",
                     "Continuous covariate(s)",
                     "Discrete covariate(s)"),
          Selected = c(input$comparisons_microarray_raw[c],
                       ifelse(is.null(input$covGroups_num_microarray_raw), " ",
                              paste(input$covGroups_num_microarray_raw, collapse = "; ")),
                       ifelse(is.null(input$covGroups_char_microarray_raw), " ",
                              paste(input$covGroups_char_microarray_raw, collapse = "; ")))
        )
      }
      names(rv$statSettings) <- input$comparisons_microarray_raw
      
    }
    
    rv$newFactor <- rv$experimentFactor
    rv$newData <- rv$normData
    
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
          selectInput(inputId = "comparisons_view_microarray_raw",
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
      rv$colorOrder_stat <- 1:length(levels(rv$newFactor))
      observeEvent(input$statboxplot_changeOrder_microarray_raw,{
        all_orders <- permute(1:length(levels(rv$newFactor))) 
        sel <- which(unlist(lapply(all_orders, function(x) all(x == rv$colorOrder_stat))))
        
        if (sel < length(all_orders)){
          rv$colorOrder_stat <- all_orders[[sel+1]]
        } else{
          rv$colorOrder_stat <- all_orders[[1]]
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
                               legendColors = legendColors[rv$colorOrder_stat],
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
    observe({
      req(input$statboxplot_file_microarray_raw)
      output$realdownload_statboxplot_microarray_raw <- downloadHandler(
        filename = ifelse(input$statboxplot_file_microarray_raw == "PNG", "GeneBoxplot.png",
                          ifelse(input$statboxplot_file_microarray_raw == "PDF", "GeneBoxplot.pdf",
                                 "GeneBoxplot.tif")),
        content = function(file){
          ggplot2::ggsave(plot = rv$temp, 
                          filename = file,
                          width = input$width_statboxplot_microarray_raw,
                          height = input$height_statboxplot_microarray_raw,
                          units = "px")
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_statboxplot_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets:prettyRadioButtons(
                     inputId = "statboxplot_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   sliderInput("height_statboxplot_microarray_raw", 
                               "Height",
                               min = 800, max = 5000,
                               value = 2100, step = 10,
                               width = "100%"),
            ),
            column(6,
                   sliderInput("width_statboxplot_microarray_raw", 
                               "Width",
                               min = 800, max = 5000,
                               value = 2100, step = 10,
                               width = "100%"),
            )
          ),
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
    # TAB2: P value and logFC histograms
    #********************************************************************#
    
    observe({
      req(input$comparisons_view_microarray_raw)
      req(rv$top_table)
      
      if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
        
        # Make p-value histogram
        rv$Phistogram <- makePHistogram(P = rv$top_table[[input$comparisons_view_microarray_raw]][,"p-value"],
                                        color = input$histogram_color_microarray_raw,
                                        bins = input$histogram_bins_microarray_raw)
        
        # Plot p-value histogram
        output$Phistogram_microarray_raw <- plotly::renderPlotly({
          return(rv$Phistogram)
        })
        
        #***************************#
        # Modal to P value histogram
        #***************************#
        
        # Download plot
        observe({
          req(input$Phistogram_file_microarray_raw)
          output$realdownload_Phistogram_microarray_raw <- downloadHandler(
            filename = ifelse(input$Phistogram_file_microarray_raw == "HTML", "Phistogram.html",
                              ifelse(input$Phistogram_file_microarray_raw == "PNG", "Phistogram.png",
                                     ifelse(input$Phistogram_file_microarray_raw == "PDF", "Phistogram.pdf",
                                            "Phistogram.tif"))),
            content = function(file){
              
              if (input$Phistogram_file_microarray_raw != "HTML"){
                
                # Make histogram
                p <- makePHistogram(P = rv$top_table[[input$comparisons_view_microarray_raw]][,"p-value"],
                                    color = input$histogram_color_microarray_raw,
                                    bins = input$histogram_bins_microarray_raw,
                                    static = TRUE)
                
                ggplot2::ggsave(plot = p, 
                                filename = file,
                                width = input$width_Phistogram_microarray_raw,
                                height = input$height_Phistogram_microarray_raw,
                                units = "px")
              } else{
                htmlwidgets::saveWidget(rv$Phistogram, 
                                        file)
              }
            }
          )
        })
        
        
        # Make modal
        observeEvent(input$download_Phistogram_microarray_raw, {
          showModal(modalDialog(
            title = NULL,
            easyClose = TRUE,
            size = "m",
            footer = tagList(
              fluidRow(
                column(6,align = "left",
                       shinyWidgets::prettyRadioButtons(
                         inputId = "Phistogram_file_microarray_raw",
                         label = NULL,
                         choices = c("PNG","PDF", "TIF", "HTML"),
                         selected = "PNG",
                         fill = TRUE
                       )
                )
              ),
              fluidRow(
                column(6,
                       conditionalPanel(
                         condition = "input.Phistogram_file_microarray_raw!=`HTML`",
                         sliderInput("height_Phistogram_microarray_raw", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                       )
                ),
                column(6,
                       conditionalPanel(
                         condition = "input.Phistogram_file_microarray_raw!=`HTML`",
                         sliderInput("width_Phistogram_microarray_raw", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                       )
                )
              ),
              
              fluidRow(
                column(12, align = "left",
                       downloadButton('realdownload_Phistogram_microarray_raw', 
                                      'Download')
                )
              )
              
            )
            
          ))
        })
        
        
        # Make logFC histogram
        rv$logFChistogram <- makelogFCHistogram(logFC = rv$top_table[[input$comparisons_view_microarray_raw]][,"log2FC"],
                                                color = input$histogram_color_microarray_raw,
                                                bins = input$histogram_bins_microarray_raw)
        
        # Plot logFC histogram
        output$logFChistogram_microarray_raw <- plotly::renderPlotly({
          return(rv$logFChistogram)
        })
        
        #********************************#
        # Modal to logFC value histogram
        #********************************#
        
        # Download plot
        observe({
          req(input$logFChistogram_file_microarray_raw)
          output$realdownload_logFChistogram_microarray_raw <- downloadHandler(
            filename = ifelse(input$logFChistogram_file_microarray_raw == "HTML", "logFChistogram.html",
                              ifelse(input$logFChistogram_file_microarray_raw == "PNG", "logFChistogram.png",
                                     ifelse(input$logFChistogram_file_microarray_raw == "PDF", "logFChistogram.pdf",
                                            "logFChistogram.tif"))),
            content = function(file){
              
              if (input$logFChistogram_file_microarray_raw != "HTML"){
                
                # Make volcano plot
                p <- makelogFCHistogram(logFC = rv$top_table[[input$comparisons_view_microarray_raw]][,"log2FC"],
                                        color = input$histogram_color_microarray_raw,
                                        bins = input$histogram_bins_microarray_raw,
                                        static = TRUE)
                
                ggplot2::ggsave(plot = p, 
                                filename = file,
                                width = input$width_logFChistogram_microarray_raw,
                                height = input$height_logFChistogram_microarray_raw,
                                units = "px")
              } else{
                htmlwidgets::saveWidget(rv$logFChistogram, 
                                        file)
              }
            }
          )
        })
        
        
        # Make modal
        observeEvent(input$download_logFChistogram_microarray_raw, {
          showModal(modalDialog(
            title = NULL,
            easyClose = TRUE,
            size = "m",
            footer = tagList(
              fluidRow(
                column(6,align = "left",
                       shinyWidgets::prettyRadioButtons(
                         inputId = "logFChistogram_file_microarray_raw",
                         label = NULL,
                         choices = c("PNG","PDF", "TIF", "HTML"),
                         selected = "PNG",
                         fill = TRUE,
                         inline = TRUE
                       )
                )
              ),
              fluidRow(
                column(6,
                       conditionalPanel(
                         condition = "input.logFChistogram_file_microarray_raw!=`HTML`",
                         sliderInput("height_logFChistogram_microarray_raw", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                       )
                ),
                column(6,
                       conditionalPanel(
                         condition = "input.logFChistogram_file_microarray_raw!=`HTML`",
                         sliderInput("width_logFChistogram_microarray_raw", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                       )
                )
              ),
              
              fluidRow(
                column(12, align = "left",
                       downloadButton('realdownload_logFChistogram_microarray_raw', 
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
    
    observeEvent(input$plot_volcano_microarray_raw, {
      req(rv$top_table)
      req(input$rawp_volcano_microarray_raw)
      req(input$p_thres_volcano_microarray_raw)
      req(input$logFC_thres_volcano_microarray_raw)
      req(input$comparisons_view_microarray_raw)
      
      if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
        rv$volcano <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                                  p = input$rawp_volcano_microarray_raw, 
                                  p_threshold = input$p_thres_volcano_microarray_raw, 
                                  logFC_threshold = input$logFC_thres_volcano_microarray_raw,
                                  unchanged_color = input$volcano_unchanged_color_microarray_raw,
                                  down_color = input$volcano_down_color_microarray_raw,
                                  up_color = input$volcano_up_color_microarray_raw
        )
        
        output$volcano_microarray_raw <- plotly::renderPlotly(rv$volcano)
      }
    }, ignoreNULL = FALSE) 
    # ignoreNULL: generate plot even if action button is not pressed
    
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$volcano_file_microarray_raw)
      output$realdownload_volcano_microarray_raw <- downloadHandler(
        filename = ifelse(input$volcano_file_microarray_raw == "HTML", "Volcano.html",
                          ifelse(input$volcano_file_microarray_raw == "PNG", "Volcano.png",
                                 ifelse(input$volcano_file_microarray_raw == "PDF", "Volcano.pdf",
                                        "Volcano.tif"))),
        content = function(file){
          
          if (input$volcano_file_microarray_raw != "HTML"){
            
            # Make volcano plot
            p <- makeVolcano_static(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                                    p = input$rawp_volcano_microarray_raw, 
                                    p_threshold = input$p_thres_volcano_microarray_raw, 
                                    logFC_threshold = input$logFC_thres_volcano_microarray_raw,
                                    unchanged_color = input$volcano_unchanged_color_microarray_raw,
                                    down_color = input$volcano_down_color_microarray_raw,
                                    up_color = input$volcano_up_color_microarray_raw)
            
            ggplot2::ggsave(plot = p, 
                            filename = file,
                            width = input$width_volcano_microarray_raw,
                            height = input$height_volcano_microarray_raw,
                            units = "px")
          } else{
            htmlwidgets::saveWidget(rv$volcano, 
                                    file)
          }
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_volcano_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "volcano_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF", "HTML"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.volcano_file_microarray_raw!=`HTML`",
                     sliderInput("height_volcano_microarray_raw", 
                                 "Height",
                                 min = 800, max = 3000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.volcano_file_microarray_raw!=`HTML`",
                     sliderInput("width_volcano_microarray_raw", 
                                 "Width",
                                 min = 800, max = 4000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_volcano_microarray_raw', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    #********************************************************************#
    # TAB4: MA plot
    #********************************************************************#
    
    observeEvent(input$plot_MA_microarray_raw, {
      req(rv$top_table)
      req(input$rawp_MA_microarray_raw)
      req(input$p_thres_MA_microarray_raw)
      req(input$logFC_thres_MA_microarray_raw)
      req(input$comparisons_view_microarray_raw)
      
      if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
        rv$MA <- makeMAplot(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                            p = input$rawp_MA_microarray_raw, 
                            p_threshold = input$p_thres_MA_microarray_raw, 
                            logFC_threshold = input$logFC_thres_MA_microarray_raw,
                            unchanged_color = input$MA_unchanged_color_microarray_raw,
                            down_color = input$MA_down_color_microarray_raw,
                            up_color = input$MA_up_color_microarray_raw,
                            RNAseq = FALSE)
        
        output$MA_microarray_raw <- plotly::renderPlotly(rv$MA)
      }
    }, ignoreNULL = FALSE) 
    # ignoreNULL: generate plot even if action button is not pressed
    
    #***************************#
    # Modal to download figure
    #***************************#
    
    # Download plot
    observe({
      req(input$MA_file_microarray_raw)
      output$realdownload_MA_microarray_raw <- downloadHandler(
        filename = ifelse(input$MA_file_microarray_raw == "HTML", "MA.html",
                          ifelse(input$MA_file_microarray_raw == "PNG", "MA.png",
                                 ifelse(input$MA_file_microarray_raw == "PDF", "MA.pdf",
                                        "MA.tif"))),
        content = function(file){
          
          if (input$MA_file_microarray_raw != "HTML"){
            
            # Make MA plot
            p <- makeMAplot_static(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                                   p = input$rawp_MA_microarray_raw, 
                                   p_threshold = input$p_thres_MA_microarray_raw, 
                                   logFC_threshold = input$logFC_thres_MA_microarray_raw,
                                   unchanged_color = input$MA_unchanged_color_microarray_raw,
                                   down_color = input$MA_down_color_microarray_raw,
                                   up_color = input$MA_up_color_microarray_raw,
                                   RNAseq =  FALSE)
            
            ggplot2::ggsave(plot = p, 
                            filename = file,
                            width = input$width_MA_microarray_raw,
                            height = input$height_MA_microarray_raw,
                            units = "px")
          } else{
            htmlwidgets::saveWidget(rv$MA, 
                                    file)
          }
        }
      )
    })
    
    
    # Make modal
    observeEvent(input$download_MA_microarray_raw, {
      showModal(modalDialog(
        title = NULL,
        easyClose = TRUE,
        size = "m",
        footer = tagList(
          fluidRow(
            column(6,align = "left",
                   shinyWidgets::prettyRadioButtons(
                     inputId = "MA_file_microarray_raw",
                     label = NULL,
                     choices = c("PNG","PDF", "TIF", "HTML"),
                     selected = "PNG",
                     fill = TRUE,
                     inline = TRUE
                   )
            )
          ),
          fluidRow(
            column(6,
                   conditionalPanel(
                     condition = "input.MA_file_microarray_raw!=`HTML`",
                     sliderInput("height_MA_microarray_raw", 
                                 "Height",
                                 min = 800, max = 3000,
                                 value = 1200, step = 10,
                                 width = "100%")
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = "input.MA_file_microarray_raw!=`HTML`",
                     sliderInput("width_MA_microarray_raw", 
                                 "Width",
                                 min = 800, max = 4000,
                                 value = 1500, step = 10,
                                 width = "100%")
                   )
            )
          ),
          
          fluidRow(
            column(12, align = "left",
                   downloadButton('realdownload_MA_microarray_raw', 
                                  'Download')
            )
          )
          
        )
        
      ))
    })
    
    
    #********************************************************************#
    # TAB5: Summary
    #********************************************************************#
    observeEvent(input$get_summary_microarray_raw, {
      req(rv$top_table)
      req(input$rawp_summary_microarray_raw)
      req(input$p_thres_summary_microarray_raw)
      req(input$logFC_thres_summary_microarray_raw)
      req(input$comparisons_view_microarray_raw)
      
      if (input$comparisons_view_microarray_raw %in% names(rv$top_table)){
        # DEG table
        temp <- matrix(NA,2,3)
        
        if (input$rawp_summary_microarray_raw == "raw"){
          rownames(temp) <- c(paste0("p-value < ", input$p_thres_summary_microarray_raw),
                              paste0("p-value > ", input$p_thres_summary_microarray_raw))
          colnames(temp) <- c(paste0("log2FC < -", input$logFC_thres_summary_microarray_raw),
                              paste0("-",input$logFC_thres_summary_microarray_raw, " < log2FC < ", input$logFC_thres_summary_microarray_raw),
                              paste0("log2FC > ", input$logFC_thres_summary_microarray_raw))
          
          temp[1,1] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` < input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC < -1*input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[1,2] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` < input$p_thres_summary_microarray_raw) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC) < input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[1,3] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` < input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC > input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,1] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` > input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC < -1*input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,2] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` > input$p_thres_summary_microarray_raw) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC) < input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,3] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` > input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC > input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
        } else{
          rownames(temp) <- c(paste0("adj. p-value < ", input$p_thres_summary_microarray_raw),
                              paste0("adj. p-value > ", input$p_thres_summary_microarray_raw))
          colnames(temp) <- c(paste0("log2FC < -", input$logFC_thres_summary_microarray_raw),
                              paste0("-",input$logFC_thres_summary_microarray_raw, " < log2FC < ", input$logFC_thres_summary_microarray_raw),
                              paste0("log2FC > ", input$logFC_thres_summary_microarray_raw))
          
          temp[1,1] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`adj. p-value` < input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC < -1*input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[1,2] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`adj. p-value` < input$p_thres_summary_microarray_raw) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC) < input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[1,3] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`adj. p-value` < input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC > input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,1] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`adj. p-value` > input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC < -1*input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,2] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`p-value` > input$p_thres_summary_microarray_raw) &
                             (abs(rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC) < input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
          temp[2,3] <- sum((rv$top_table[[input$comparisons_view_microarray_raw]]$`adj. p-value` > input$p_thres_summary_microarray_raw) &
                             (rv$top_table[[input$comparisons_view_microarray_raw]]$log2FC > input$logFC_thres_summary_microarray_raw), na.rm = TRUE)
        }
        rv$summaryTable <- temp
        rm(temp)
        
        output$summaryTable_microarray_raw <- DT::renderDataTable({
          return(rv$summaryTable)
        },options = list(pageLength = 10,
                         dom = 't'),
        selection = "none")
        
      }
    }, ignoreNULL = FALSE)
    
    
    # Download button
    output$downloadSummaryTable_microarray_raw <- downloadHandler(
      filename = "SummaryTable.csv",
      content = function(file){
        write.csv(rv$summaryTable, file, quote = FALSE, row.names = TRUE)
      }
    )
    
    #********************************************************************#
    # TAB6: Settings
    #********************************************************************#
    observe({
      req(input$comparisons_view_microarray_raw)
      
      # Print table with settings
      output$statSettings_microarray_raw <- DT::renderDataTable({
        return(rv$statSettings[[input$comparisons_view_microarray_raw]])
      },options = list(pageLength = 10,
                       dom = 't',
                       rownames = FALSE),
      selection = "none")
      
      # Download button
      output$downloadStatSettings_microarray_raw <- downloadHandler(
        filename = "StatisticalAnalysis_Settings.csv",
        content = function(file){
          write.csv(rv$statSettings[[input$comparisons_view_microarray_raw]], file, quote = FALSE, row.names = FALSE)
        }
      )
      
    })
    
    #********************************************************************#
    # Statistical analysis report
    #********************************************************************#
    output$SAreport_microarray_raw <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "SAreport.html",
      content = function(file) {
        shinybusy::show_modal_spinner(text = "Making statistical analysis report...",
                                      color="#0dc5c1")
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "SAreport_microarray_raw.Rmd")
        file.copy("Reports/SAreport_microarray_raw.Rmd", tempReport, overwrite = TRUE)
        
        tempLogo <- file.path(tempdir(), "logo_main.PNG")
        file.copy("www/logo_main.PNG", tempLogo, overwrite = TRUE)
        
        tempHeader <- file.path(tempdir(), "header.html")
        file.copy("www/header.html", tempHeader, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(statSettings = rv$statSettings[[input$comparisons_view_microarray_raw]],
                       topTable = rv$top_table[[input$comparisons_view_microarray_raw]],
                       volcanoTable = rv$summaryTable,
                       dir = tempdir(),
                       ArrayAnalysis_version = ArrayAnalysis_version)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
        shinybusy::remove_modal_spinner()
      }
    )
    
    #--------------------------------------#
    # dynamic UI: Output in different tabs
    #--------------------------------------#
    
    observe({
      req(rv$top_table)
      
      output$UI_boxplotAll_microarray_raw <- renderUI({
        tagList(
          tags$div(class = "dropdown-content",
                   
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
                       h5("The top table includes the output of the statistical analysis. 
                                  Click on the table to explore the data!"),
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
                       tags$div(class = "dropdown",
                                tags$button(class = "dropbtn", icon("cog")),
                                uiOutput("UI_boxplotAll_microarray_raw")
                       ), # EO dropdownButton
                       
                       plotOutput("ExprBoxplot_statistics_microarray_raw")%>% 
                         shinycssloaders::withSpinner(color="#0dc5c1"),
                       # Plot all experimental groups?
                       conditionalPanel(
                         condition = "output.length_statboxplot_microarray_raw > 2",
                         br(),
                         shinyWidgets::materialSwitch(inputId = "boxplotAll_microarray_raw",
                                                      label = "Plot all experimental groups?", 
                                                      value =  TRUE,
                                                      status = "primary"),
                         
                       ),
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
                       # Dropdown Button to adjust the plot settings
                       tags$div(class = "dropdown",
                                tags$button(class = "dropbtn", icon("cog")),
                                tags$div(class = "dropdown-content",
                                         tags$h4("Select color"),
                                         colourpicker::colourInput("histogram_color_microarray_raw", 
                                                                   NULL, 
                                                                   "#d3d3d3"),
                                         br(),
                                         tags$h4("Number of bins"),
                                         numericInput(inputId = "histogram_bins_microarray_raw",
                                                      label = NULL,
                                                      value = 100),
                                         br(),br()
                                )), # EO dropdownButton
                       plotly::plotlyOutput("Phistogram_microarray_raw")%>% 
                         shinycssloaders::withSpinner(color="#0dc5c1"),
                       actionButton("download_Phistogram_microarray_raw", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       hr(),
                       plotly::plotlyOutput("logFChistogram_microarray_raw")%>% 
                         shinycssloaders::withSpinner(color="#0dc5c1"),
                       actionButton("download_logFChistogram_microarray_raw", 
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
                       tags$div(class = "dropdown",
                                tags$button(class = "dropbtn", icon("cog")),
                                tags$div(class = "dropdown-content",
                                         tags$h4("Color of unchanged genes"),
                                         colourpicker::colourInput("volcano_unchanged_color_microarray_raw", 
                                                                   NULL, 
                                                                   "darkgrey"),
                                         tags$h4("Color of downregulated genes"),
                                         colourpicker::colourInput("volcano_down_color_microarray_raw", 
                                                                   NULL, 
                                                                   "blue"),
                                         tags$h4("Color of upregulated genes"),
                                         colourpicker::colourInput("volcano_up_color_microarray_raw", 
                                                                   NULL, 
                                                                   "red"),
                                         br()
                                )), # EO dropdownButton
                       
                       # Volcano plot output
                       plotly::plotlyOutput("volcano_microarray_raw")%>% 
                         withSpinner(color="#0dc5c1"),
                       
                       
                       br(),
                       actionButton("download_volcano_microarray_raw", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                shinyWidgets::prettyRadioButtons(
                                  inputId = "rawp_volcano_microarray_raw",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_volcano_microarray_raw",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "plot_volcano_microarray_raw", 
                                                         label = "Load",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_volcano_microarray_raw",
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
                       tags$div(class = "dropdown",
                                tags$button(class = "dropbtn", icon("cog")),
                                tags$div(class = "dropdown-content",
                                         tags$h4("Color of unchanged genes"),
                                         colourpicker::colourInput("MA_unchanged_color_microarray_raw", 
                                                                   NULL, 
                                                                   "darkgrey"),
                                         tags$h4("Color of downregulated genes"),
                                         colourpicker::colourInput("MA_down_color_microarray_raw", 
                                                                   NULL, 
                                                                   "blue"),
                                         tags$h4("Color of upregulated genes"),
                                         colourpicker::colourInput("MA_up_color_microarray_raw", 
                                                                   NULL, 
                                                                   "red"),
                                         br()
                                )), # EO dropdownButton
                       
                       # MA plot output
                       plotly::plotlyOutput("MA_microarray_raw")%>% 
                         withSpinner(color="#0dc5c1"),
                       br(),
                       actionButton("download_MA_microarray_raw", 
                                    "Download figure",
                                    icon = shiny::icon("download")),
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                shinyWidgets::prettyRadioButtons(
                                  inputId = "rawp_MA_microarray_raw",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_MA_microarray_raw",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "plot_MA_microarray_raw", 
                                                         label = "Load",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_MA_microarray_raw",
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
                       DT::dataTableOutput(outputId = "summaryTable_microarray_raw") %>%
                         withSpinner(color="#0dc5c1"),
                       br(),
                       downloadButton("downloadSummaryTable_microarray_raw",
                                      "Download table"),
                       
                       br(),
                       hr(),
                       fluidRow(
                         column(2,
                                # Raw or adjusted P value?
                                shinyWidgets::prettyRadioButtons(
                                  inputId = "rawp_summary_microarray_raw",
                                  label = "P value", 
                                  choices = 
                                    c("p-value" = "raw", 
                                      "Adj. p-value" = "adj"))
                         ),
                         column(3,
                                #P value threshold
                                numericInput(
                                  inputId = "p_thres_summary_microarray_raw",
                                  label = "P threshold",
                                  value = 0.05,
                                  width = "80%")
                         )
                       ),
                       fluidRow(
                         column(2,
                                br(),
                                # Button to reload the volcano plot
                                shinyWidgets::actionBttn(inputId = "get_summary_microarray_raw", 
                                                         label = "Calculate",
                                                         style = "jelly",
                                                         color = "primary",
                                                         icon = icon("sync"))
                         ),
                         column(3,
                                #logFC threshold
                                numericInput(
                                  inputId = "logFC_thres_summary_microarray_raw",
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
                       h5("To enhance reproducibility, download the overview of chosen statistical analysis settings."),
                       hr(),
                       DT::dataTableOutput(outputId = "statSettings_microarray_raw") %>% 
                         withSpinner(color="#0dc5c1"),
                       br(),
                       downloadButton("downloadStatSettings_microarray_raw", 
                                      "Download table"),
                       downloadButton("downloadSessionInfo_microarray_raw", 
                                      "Session info")
                       
              ) # EO Settings tabPanel
              
              
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
        shinyWidgets::downloadBttn(outputId = "SAreport_microarray_raw",
                                   label = "Get report",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("download")),
        shinyWidgets::actionBttn(inputId = "next_statistics_microarray_raw",
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
  observeEvent(input$next_statistics_microarray_raw,{
    
    # Go to microarray statistics tab
    updateNavbarPage(session, "navbar",
                     selected = "panel_ORA_microarray_raw")
  })
  
  #********************************************************************#
  # Options for gene set analysis
  #********************************************************************#
  
  # Comparisons for which gene set analysis should be performed
  observe({
    req(rv$top_table)
    output$UI_comparisons_view_ORA_microarray_raw <- renderUI({
      selectInput(inputId = "comparisons_view_ORA_microarray_raw",
                  label = NULL,
                  choices = names(rv$top_table),
                  selected = names(rv$top_table)[1],
                  multiple = FALSE)
    })
  })
  
  # Gene ID information
  observe({
    
    # Get all columns of the top table that can possibly contain the gene IDs
    #req(input$comparisons_view_ORA_microarray_raw)
    col_choice <- 1
    if (length(colnames(rv$top_table[[1]])) > 6){
      col_choice <- c(1,7:ncol(rv$top_table[[1]]))
    }
    
    # Several options need to be selected before ORA can be performed
    output$UI_geneID_ORA_microarray_raw <- renderUI({
      tagList(
        
        # Which columns of the top table contains the gene ids?
        selectInput(inputId = "geneID_ORA_microarray_raw",
                    label = "Which column of the top table contains the gene IDs?",
                    choices = colnames(rv$top_table[[1]])[col_choice],
                    selected = colnames(rv$top_table[[1]])[1],
                    multiple = FALSE)
      )
    })
  })
  
  
  observeEvent(input$calculate_ORA_microarray_raw,{
    
    #==========================================================================#
    
    # ORA
    
    #==========================================================================#
    
    if (input$ORA_or_GSEA_microarray_raw == "ORA"){
      # Show modal
      shinybusy::show_modal_spinner(text = "Overrepresentation analysis...",
                                    color="#0dc5c1")
      
      # Get gene set version
      if (input$geneset_ORA_microarray_raw == "WikiPathways"){
        load("Objects/WPversion.RData")
        rv$GeneSetVersion <- WPversion
      } else {
        pkg <- switch(input$organism_ORA_microarray_raw,
                      "Homo sapiens" = "org.Hs.eg.db",
                      "Bos taurus" = "org.Bt.eg.db",
                      "Caenorhabditis elegans" = "org.Ce.eg.db",
                      "Mus musculus" = "org.Mm.eg.db",
                      "Rattus norvegicus" = "org.Rn.eg.db"
        )
        rv$GeneSetVersion <- paste0(pkg, " v", packageVersion(pkg))
      }
      
      # Perform ORA:
      
      # Perform ORA based on logFC/P value threshold(s)
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
        
        
        rv$ORA_settings <- data.frame(
          Option = c("Comparison",
                     "Gene set collection",
                     "Gene ID",
                     "Organism",
                     "Method",
                     "p-value threshold",
                     "logFC threshold"
          ),
          Selected = c(input$comparisons_view_ORA_microarray_raw,
                       paste0(input$geneset_ORA_microarray_raw, " (", rv$GeneSetVersion, ")"),
                       paste0(input$selID_ORA_microarray_raw, " (top table column: ", input$geneID_ORA_microarray_raw, ")"),
                       input$organism_ORA_microarray_raw,
                       paste0("ORA (", input$topNorThres_microarray_raw, "; ", input$updown_ORA_microarray_raw, ")"),
                       paste0(input$p_thres_ORA_microarray_raw, " (", input$rawp_ORA_microarray_raw, ")"),
                       input$logFC_thres_ORA_microarray_raw
          )
        )
        
        # Perform ORA based on top N most significant genes
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
        
        
        rv$ORA_settings <- data.frame(
          Option = c("Comparison",
                     "Gene set collection",
                     "Gene ID",
                     "Organism",
                     "Method",
                     "Top N"
          ),
          Selected = c(input$comparisons_view_ORA_microarray_raw,
                       paste0(input$geneset_ORA_microarray_raw, " (", rv$GeneSetVersion, ")"),
                       paste0(input$selID_ORA_microarray_raw, " (top table column: ", input$geneID_ORA_microarray_raw, ")"),
                       input$organism_ORA_microarray_raw,
                       paste0("ORA (",input$topNorThres_microarray_raw, "; ", input$updown_ORA_microarray_raw, ")"),
                       input$topN_microarray_raw
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
          
          output$ORA_table_microarray_raw <- DT::renderDataTable({
            req(input$geneset_ORA_microarray_raw)
            req(rv$ORA_data)
            output <- rv$ORA_data@result
            
            # Link to wikipathways website if gene set ID is from WikiPathways
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
            
            # Link to QuickGO website if gene set ID is a GO term
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
          output$ORAgene_table_microarray_raw <- DT::renderDataTable({
            req(input$ORA_table_microarray_raw_rows_selected)
            req(rv$ORA_data)
            
            # Make ORA gene table
            output <- make_ORAgene_table(ORA_data = rv$ORA_data,
                                         top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                         geneID_col = input$geneID_ORA_microarray_raw,
                                         sel_row_ORA = input$ORA_table_microarray_raw_rows_selected)
            
            output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
            output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
            output$meanExpr <- round(output$meanExpr,3)
            output$log2FC <- round(output$log2FC,3)
            output$`log2FC SE` <- round(output$`log2FC SE`,3)
            
            return(output)
          }, options = list(pageLength = 6), escape = FALSE)
          
          # Text for gene table
          output$text_ORAgene_table_microarray_raw <- renderText({
            req(rv$ORA_data)
            text <- paste0("<h3><b>Gene table: ",rv$ORA_data@result[input$ORA_table_microarray_raw_rows_selected,"ID"],
                           "</b></h3>")
            return(text)
          })
          
          #--------------------------------------------------------------#
          # ORA barchart
          #--------------------------------------------------------------#
          
          observe({
            req(input$nSets_ORAplot_microarray_raw)
            req(rv$ORA_data)
            rv$ORAplot <- makeORAplot(rv$ORA_data,
                                      nSets = input$nSets_ORAplot_microarray_raw,
                                      color = input$color_ORAplot_microarray_raw)
            
            output$ORAplot_microarray_raw <- plotly::renderPlotly(rv$ORAplot%>% 
                                                                    plotly::layout(height = 500, width = 1000))
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          observe({
            req(input$ORAplot_file_microarray_raw)
            output$realdownload_ORAplot_microarray_raw <- downloadHandler(
              filename = ifelse(input$ORAplot_file_microarray_raw == "HTML", "ORA_barchart.html",
                                ifelse(input$ORAplot_file_microarray_raw == "PNG", "ORA_barchart.png",
                                       ifelse(input$ORAplot_file_microarray_raw == "PDF", "ORA_barchart.pdf",
                                              "ORA_barchart.tif"))),
              content = function(file){
                
                if (input$ORAplot_file_microarray_raw != "HTML"){
                  
                  
                  # Make MA plot
                  p <- makeORAplot(rv$ORA_data,
                                   nSets = input$nSets_ORAplot_microarray_raw,
                                   color = input$color_ORAplot_microarray_raw,
                                   static = TRUE)
                  
                  ggplot2::ggsave(plot = p, 
                                  filename = file,
                                  width = input$width_ORAplot_microarray_raw,
                                  height = input$height_ORAplot_microarray_raw,
                                  units = "px")
                } else{
                  htmlwidgets::saveWidget(rv$ORAplot, 
                                          file)
                }
              }
            )
          })
          
          
          # Make modal
          observeEvent(input$download_ORAplot_microarray_raw, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,align = "left",
                         shinyWidgets::prettyRadioButtons(
                           inputId = "ORAplot_file_microarray_raw",
                           label = NULL,
                           choices = c("PNG","PDF", "TIF", "HTML"),
                           selected = "PNG",
                           fill = TRUE,
                           inline = TRUE
                         )
                  )
                ),
                fluidRow(
                  column(6,
                         conditionalPanel(
                           condition = "input.ORAplot_file_microarray_raw!=`HTML`",
                           sliderInput("height_ORAplot_microarray_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 1200, step = 10,
                                       width = "100%")
                         )
                  ),
                  column(6,
                         conditionalPanel(
                           condition = "input.ORAplot_file_microarray_raw!=`HTML`",
                           sliderInput("width_ORAplot_microarray_raw", 
                                       "Width",
                                       min = 800, max = 4000,
                                       value = 1500, step = 10,
                                       width = "100%")
                         )
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_ORAplot_microarray_raw', 
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
            req(input$layout_ORAnetwork_microarray_raw)
            req(input$nSets_ORAnetwork_microarray_raw)
            req(rv$ORA_data)
            rv$ORAnetwork <- makeORAnetwork(ORA_data = rv$ORA_data,
                                            layout = input$layout_ORAnetwork_microarray_raw,
                                            nSets = input$nSets_ORAnetwork_microarray_raw,
                                            color = input$color_ORAnetwork_microarray_raw)
            
            output$ORAnetwork_microarray_raw <- renderPlot(rv$ORAnetwork,
                                                           height = 500, 
                                                           width = 800)
            
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          observe({
            req(input$ORAnetwork_file_microarray_raw)
            output$realdownload_ORAnetwork_microarray_raw <- downloadHandler(
              filename = ifelse(input$ORAnetwork_file_microarray_raw == "PNG", "ORA_network.png",
                                ifelse(input$ORAnetwork_file_microarray_raw == "PDF", "ORA_network.pdf",
                                       "ORA_network.tif")),
              content = function(file){
                
                ggplot2::ggsave(plot = rv$ORAnetwork, 
                                filename = file,
                                width = input$width_ORAnetwork_microarray_raw*2,
                                height = input$height_ORAnetwork_microarray_raw*2,
                                units = "px")
              }
            )
          })
          
          
          # Make modal
          observeEvent(input$download_ORAnetwork_microarray_raw, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,align = "left",
                         shinyWidgets::prettyRadioButtons(
                           inputId = "ORAnetwork_file_microarray_raw",
                           label = NULL,
                           choices = c("PNG","PDF", "TIF"),
                           selected = "PNG",
                           inline = TRUE,
                           fill = TRUE
                         )
                  )
                ),
                fluidRow(
                  column(6,
                         sliderInput("height_ORAnetwork_microarray_raw", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                         
                  ),
                  column(6,
                         sliderInput("width_ORAnetwork_microarray_raw", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_ORAnetwork_microarray_raw', 
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
            output$ORASettings_microarray_raw <- DT::renderDataTable({
              return(rv$ORA_settings)
            },options = list(pageLength = 10,
                             dom = 't',
                             rownames = FALSE),
            selection = "none")
            
            # Download button
            output$downloadORASettings_microarray_raw <- downloadHandler(
              filename = "ORA_Settings.csv",
              content = function(file){
                write.csv(rv$ORA_settings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
          })
          
          #--------------------------------------------------------------#
          # ORA report
          #--------------------------------------------------------------#
          
          output$ORAreport_microarray_raw <- downloadHandler(
            # For PDF output, change this to "report.pdf"
            filename = "ORAreport.html",
            content = function(file) {
              shinybusy::show_modal_spinner(text = "Making ORA report...",
                                           color="#0dc5c1")
              # Copy the report file to a temporary directory before processing it, in
              # case we don't have write permissions to the current working dir (which
              # can happen when deployed).
              tempReport <- file.path(tempdir(), "ORAreport_microarray_raw.Rmd")
              file.copy("Reports/ORAreport_microarray_raw.Rmd", tempReport, overwrite = TRUE)
              
              tempLogo <- file.path(tempdir(), "logo_main.PNG")
              file.copy("www/logo_main.PNG", tempLogo, overwrite = TRUE)
              
              tempHeader <- file.path(tempdir(), "header.html")
              file.copy("www/header.html", tempHeader, overwrite = TRUE)
              
              # Set up parameters to pass to Rmd document
              params <- list(ORASettings = rv$ORA_settings,
                             ORATable = rv$ORA_data,
                             dir = tempdir(),
                             ArrayAnalysis_version = ArrayAnalysis_version)
              
              # Knit the document, passing in the `params` list, and eval it in a
              # child of the global environment (this isolates the code in the document
              # from the code in this app).
              rmarkdown::render(tempReport, output_file = file,
                                params = params,
                                envir = new.env(parent = globalenv())
              )
              shinybusy::remove_modal_spinner()
            }
          )
          
        } # EO ifelse
      }) # EO observe
      
      
      
      
      
      #********************************************************************#
      # UI: ORA output
      #********************************************************************#
      observe({
        if (!is.null(rv$ORA_data)){
          output$UI_output_ORA_microarray_raw <- renderUI({
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
                         h5("Here you can view the statistics from the ORA. 
                            Click on the table to get the statistics per gene."),
                         hr(),
                         
                         # Statistics table
                         DT::dataTableOutput(outputId = "ORA_table_microarray_raw") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         
                         # Download button
                         downloadButton("download_ORA_table_microarray_raw", 
                                        "Download"),
                         br(),
                         
                         # Title + description of gene table
                         htmlOutput("text_ORAgene_table_microarray_raw"),
                         h5("Here you can view the differential expression results for the selected gene set."),
                         hr(),
                         
                         # Gene table
                         DT::dataTableOutput(outputId = "ORAgene_table_microarray_raw") %>% 
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
                         actionButton("download_ORAplot_microarray_raw", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         tags$div(class = "dropdown",
                                  tags$button(class = "dropbtn", icon("cog")),
                                  tags$div(class = "dropdown-content",
                                           
                                           # Color gradient
                                           tags$h4("Color gradient"),
                                           selectInput(inputId = "color_ORAplot_microarray_raw",
                                                       label = NULL,
                                                       choices = c("Viridis", 
                                                                   "Yellow-red", 
                                                                   "Blues", 
                                                                   "Reds")),
                                           br(),
                                           
                                           # Number of gene sets
                                           tags$h4("# Gene sets"),
                                           sliderInput(
                                             inputId = "nSets_ORAplot_microarray_raw",
                                             label = NULL,
                                             value = 10,
                                             min = 5,
                                             max = 20,
                                             step = 1),
                                           
                                           br(),br()
                                  )), # EO dropdownButton
                         
                         # Interactive plot output
                         plotly::plotlyOutput("ORAplot_microarray_raw") %>% 
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
                         h5("The network visualizes the similarity between the most significant gene sets."),
                         hr(),
                         actionButton("download_ORAnetwork_microarray_raw", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         tags$div(class = "dropdown",
                                  tags$button(class = "dropbtn", icon("cog")),
                                  tags$div(class = "dropdown-content",
                                           
                                           # Color gradient
                                           tags$h4("Color gradient"),
                                           selectInput(inputId = "color_ORAnetwork_microarray_raw",
                                                       label = NULL,
                                                       choices = c("Viridis", 
                                                                   "Yellow-red", 
                                                                   "Blues", 
                                                                   "Reds")),
                                           br(),
                                           
                                           # Network layout
                                           tags$h4("Network layout"),
                                           selectInput(inputId = "layout_ORAnetwork_microarray_raw",
                                                       label = NULL,
                                                       choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                                                                   'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                       selected = 'graphopt',
                                                       multiple = FALSE),
                                           br(),
                                           
                                           # Number of gene sets
                                           tags$h4("# Gene sets"),
                                           sliderInput(
                                             inputId = "nSets_ORAnetwork_microarray_raw",
                                             label = NULL,
                                             value = 10,
                                             min = 5,
                                             max = 20,
                                             step = 1),
                                           br(),br()
                                  )), # EO dropdownButton
                         
                         
                         
                         # Make plot
                         plotOutput("ORAnetwork_microarray_raw") %>% 
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         br(),
                         
                ),
                
                #--------------------------------------------------------------#
                # Settings tab
                #--------------------------------------------------------------#
                tabPanel("Settings overview",
                         icon = icon("fas fa-file"),
                         h3(strong("ORA settings")),
                         h5("To enhance reproducibility, download the overview of chosen ORA settings."),
                         hr(),
                         DT::dataTableOutput(outputId = "ORASettings_microarray_raw") %>% 
                           withSpinner(color="#0dc5c1"),
                         br(),
                         downloadButton("downloadORASettings_microarray_raw", 
                                        "Download table"),
                         downloadButton("downloadSessionInfo_microarray_raw", 
                                        "Session info")
                         
                )
                
              ) # EO tabSetPanel
            ) # EO tagList
          }) # EO renderUI
        }# EO if !is.null(rv$ORA_data)
        
        # Allow user to download ORA report
        output$UI_ORAreport_microarray_raw <- renderUI({
          req(rv$ORA_data)
          tagList(
            shinyWidgets::downloadBttn(outputId = "ORAreport_microarray_raw",
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
    
    if (input$ORA_or_GSEA_microarray_raw == "GSEA"){
      
      # Show modal
      shinybusy::show_modal_spinner(text = "Gene Set Enrichment Analysis...",
                                    color="#0dc5c1")
      
      # Get gene set version
      if (input$geneset_ORA_microarray_raw == "WikiPathways"){
        load("Objects/WPversion.RData")
        rv$GeneSetVersion <- WPversion
      } else {
        pkg <- switch(input$organism_ORA_microarray_raw,
                      "Homo sapiens" = "org.Hs.eg.db",
                      "Bos taurus" = "org.Bt.eg.db",
                      "Caenorhabditis elegans" = "org.Ce.eg.db",
                      "Mus musculus" = "org.Mm.eg.db",
                      "Rattus norvegicus" = "org.Rn.eg.db"
        )
        rv$GeneSetVersion <- paste0(pkg, " v", packageVersion(pkg))
      }
      
      # Perform GSEA:
      rv$GSEA_data <- performGSEA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                  geneset = input$geneset_ORA_microarray_raw,
                                  geneID_col = input$geneID_ORA_microarray_raw,
                                  geneID_type = input$selID_ORA_microarray_raw,
                                  organism = input$organism_ORA_microarray_raw,
                                  rankingVar = input$ranking_GSEA_microarray_raw)
      
      if (input$ranking_GSEA_microarray_raw == "pvalue"){
        rankvar <- "-log p-value"
      }
      if (input$ranking_GSEA_microarray_raw == "signed_pvalue"){
        rankvar <- "-log p-value x sign logFC"
      }
      if (input$ranking_GSEA_microarray_raw == "logFC"){
        rankvar <- "logFC"
      }
      
      rv$GSEA_settings <- data.frame(
        Option = c("Comparison",
                   "Gene set collection",
                   "Gene ID",
                   "Organism",
                   "Method",
                   "Ranking variable"
        ),
        Selected = c(input$comparisons_view_ORA_microarray_raw,
                     paste0(input$geneset_ORA_microarray_raw, " (", rv$GeneSetVersion, ")"),
                     paste0(input$selID_ORA_microarray_raw, " (top table column: ", input$geneID_ORA_microarray_raw, ")"),
                     input$organism_ORA_microarray_raw,
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
            text = "Gene Set Enrichment Analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
            type = "info")
          
          #--------------------------------------------------------------#
          # GSEA statistics table
          #--------------------------------------------------------------#
          
          output$GSEA_table_microarray_raw <- DT::renderDataTable({
            req(input$geneset_ORA_microarray_raw)
            req(rv$GSEA_data)
            output <- rv$GSEA_data@result
            
            # Link to wikipathways website if gene set ID is from WikiPathways
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
            
            # Link to QuickGO website if gene set ID is a GO term
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
          output$download_GSEA_table_microarray_raw <- downloadHandler(
            filename = paste0("GSEATable_",input$comparisons_view_ORA_microarray_raw,"_",input$geneset_ORA_microarray_raw,".csv"),
            content = function(file){
              write.csv(rv$GSEA_data@result, file, quote = FALSE, row.names = FALSE)
            }
          )
          
          # Print statistics of genes in selected Term
          output$GSEAgene_table_microarray_raw <- DT::renderDataTable({
            req(input$GSEA_table_microarray_raw_rows_selected)
            req(rv$GSEA_data)
            
            # Make GSEA gene table
            output <- make_ORAgene_table(ORA_data = rv$GSEA_data,
                                         top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                                         geneID_col = input$geneID_ORA_microarray_raw,
                                         sel_row_ORA = input$GSEA_table_microarray_raw_rows_selected)
            
            output$`p-value` <- format(output$`p-value`, scientific=TRUE, digits = 3)
            output$`adj. p-value` <- format(output$`adj. p-value`, scientific=TRUE, digits = 3)
            output$meanExpr <- round(output$meanExpr,3)
            output$log2FC <- round(output$log2FC,3)
            output$`log2FC SE` <- round(output$`log2FC SE`,3)
            
            return(output)
          }, options = list(pageLength = 6), escape = FALSE)
          
          # Text for gene table
          output$text_GSEAgene_table_microarray_raw <- renderText({
            req(rv$GSEA_data)
            text <- paste0("<h3><b>Gene table: ",rv$GSEA_data@result[input$GSEA_table_microarray_raw_rows_selected,"ID"],
                           "</b></h3>")
            return(text)
          })
          
          #--------------------------------------------------------------#
          # GSEA barchart
          #--------------------------------------------------------------#
          
          observe({
            req(input$nSets_GSEAplot_microarray_raw)
            req(rv$GSEA_data)
            rv$GSEAplot <- makeGSEAplot(rv$GSEA_data,
                                        nSets = input$nSets_GSEAplot_microarray_raw,
                                        color = c(input$lowcol_GSEAplot_microarray_raw,
                                                  input$midcol_GSEAplot_microarray_raw,
                                                  input$highcol_GSEAplot_microarray_raw))
            
            output$GSEAplot_microarray_raw <- plotly::renderPlotly(rv$GSEAplot%>% 
                                                                     plotly::layout(height = 500, width = 1000))
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          observe({
            req(input$GSEAplot_file_microarray_raw)
            output$realdownload_GSEAplot_microarray_raw <- downloadHandler(
              filename = ifelse(input$GSEAplot_file_microarray_raw == "HTML", "GSEA_barchart.html",
                                ifelse(input$GSEAplot_file_microarray_raw == "PNG", "GSEA_barchart.png",
                                       ifelse(input$GSEAplot_file_microarray_raw == "PDF", "GSEA_barchart.pdf",
                                              "GSEA_barchart.tif"))),
              content = function(file){
                
                if (input$GSEAplot_file_microarray_raw != "HTML"){
                  
                  
                  # Make MA plot
                  p <- makeGSEAplot(rv$GSEA_data,
                                    nSets = input$nSets_GSEAplot_microarray_raw,
                                    color = input$color_GSEAplot_microarray_raw,
                                    static = TRUE)
                  
                  ggplot2::ggsave(plot = p, 
                                  filename = file,
                                  width = input$width_GSEAplot_microarray_raw,
                                  height = input$height_GSEAplot_microarray_raw,
                                  units = "px")
                } else{
                  htmlwidgets::saveWidget(rv$GSEAplot, 
                                          file)
                }
              }
            )
          })
          
          
          # Make modal
          observeEvent(input$download_GSEAplot_microarray_raw, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,align = "left",
                         shinyWidgets::prettyRadioButtons(
                           inputId = "GSEAplot_file_microarray_raw",
                           label = NULL,
                           choices = c("PNG","PDF", "TIF", "HTML"),
                           selected = "PNG",
                           fill = TRUE,
                           inline = TRUE
                         )
                  )
                ),
                fluidRow(
                  column(6,
                         conditionalPanel(
                           condition = "input.GSEAplot_file_microarray_raw!=`HTML`",
                           sliderInput("height_GSEAplot_microarray_raw", 
                                       "Height",
                                       min = 800, max = 3000,
                                       value = 1200, step = 10,
                                       width = "100%")
                         )
                  ),
                  column(6,
                         conditionalPanel(
                           condition = "input.GSEAplot_file_microarray_raw!=`HTML`",
                           sliderInput("width_GSEAplot_microarray_raw", 
                                       "Width",
                                       min = 800, max = 4000,
                                       value = 1500, step = 10,
                                       width = "100%")
                         )
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_GSEAplot_microarray_raw', 
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
            req(input$layout_GSEAnetwork_microarray_raw)
            req(input$nSets_GSEAnetwork_microarray_raw)
            req(rv$GSEA_data)
            rv$GSEAnetwork <- makeGSEAnetwork(GSEA_data = rv$GSEA_data,
                                              layout = input$layout_GSEAnetwork_microarray_raw,
                                              nSets = input$nSets_GSEAnetwork_microarray_raw,
                                              color = c(input$lowcol_GSEAnetwork_microarray_raw,
                                                        input$midcol_GSEAnetwork_microarray_raw,
                                                        input$highcol_GSEAnetwork_microarray_raw))
            
            output$GSEAnetwork_microarray_raw <- renderPlot(rv$GSEAnetwork,
                                                            height = 500, 
                                                            width = 800)
            
          })
          
          #***************************#
          # Modal to download figure
          #***************************#
          
          # Download plot
          observe({
            req(input$GSEAnetwork_file_microarray_raw)
            output$realdownload_GSEAnetwork_microarray_raw <- downloadHandler(
              filename = ifelse(input$GSEAnetwork_file_microarray_raw == "PNG", "GSEA_network.png",
                                ifelse(input$GSEAnetwork_file_microarray_raw == "PDF", "GSEA_network.pdf",
                                       "GSEA_network.tif")),
              content = function(file){
                
                ggplot2::ggsave(plot = rv$GSEAnetwork, 
                                filename = file,
                                width = input$width_GSEAnetwork_microarray_raw*2,
                                height = input$height_GSEAnetwork_microarray_raw*2,
                                units = "px")
              }
            )
          })
          
          
          # Make modal
          observeEvent(input$download_GSEAnetwork_microarray_raw, {
            showModal(modalDialog(
              title = NULL,
              easyClose = TRUE,
              size = "m",
              footer = tagList(
                fluidRow(
                  column(6,align = "left",
                         shinyWidgets::prettyRadioButtons(
                           inputId = "GSEAnetwork_file_microarray_raw",
                           label = NULL,
                           choices = c("PNG","PDF", "TIF"),
                           selected = "PNG",
                           fill = TRUE,
                           inline = TRUE
                         )
                  )
                ),
                fluidRow(
                  column(6,
                         sliderInput("height_GSEAnetwork_microarray_raw", 
                                     "Height",
                                     min = 800, max = 2000,
                                     value = 1200, step = 10,
                                     width = "100%")
                         
                  ),
                  column(6,
                         sliderInput("width_GSEAnetwork_microarray_raw", 
                                     "Width",
                                     min = 800, max = 2000,
                                     value = 1500, step = 10,
                                     width = "100%")
                  )
                ),
                
                fluidRow(
                  column(12, align = "left",
                         downloadButton('realdownload_GSEAnetwork_microarray_raw', 
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
            output$GSEASettings_microarray_raw <- DT::renderDataTable({
              return(rv$GSEA_settings)
            },options = list(pageLength = 10,
                             dom = 't',
                             rownames = FALSE),
            selection = "none")
            
            # Download button
            output$downloadGSEASettings_microarray_raw <- downloadHandler(
              filename = "GSEA_Settings.csv",
              content = function(file){
                write.csv(rv$GSEA_settings, file, quote = FALSE, row.names = FALSE)
              }
            )
            
          })
          
          
          #--------------------------------------------------------------#
          # GSEA report
          #--------------------------------------------------------------#
          
          output$GSEAreport_microarray_raw <- downloadHandler(
            # For PDF output, change this to "report.pdf"
            filename = "GSEAreport.html",
            content = function(file) {
              shinybusy::show_modal_spinner(text = "Making GSEA report...",
                                            color="#0dc5c1")
              # Copy the report file to a temporary directory before processing it, in
              # case we don't have write permissions to the current working dir (which
              # can happen when deployed).
              tempReport <- file.path(tempdir(), "GSEAreport_microarray_raw.Rmd")
              file.copy("Reports/GSEAreport_microarray_raw.Rmd", tempReport, overwrite = TRUE)
              
              tempLogo <- file.path(tempdir(), "logo_main.PNG")
              file.copy("www/logo_main.PNG", tempLogo, overwrite = TRUE)
              
              tempHeader <- file.path(tempdir(), "header.html")
              file.copy("www/header.html", tempHeader, overwrite = TRUE)
              
              # Set up parameters to pass to Rmd document
              params <- list(GSEASettings = rv$GSEA_settings,
                             GSEATable = rv$GSEA_data,
                             dir = tempdir(),
                             ArrayAnalysis_version = ArrayAnalysis_version)
              
              # Knit the document, passing in the `params` list, and eval it in a
              # child of the global environment (this isolates the code in the document
              # from the code in this app).
              rmarkdown::render(tempReport, output_file = file,
                                params = params,
                                envir = new.env(parent = globalenv())
              )
              shinybusy::remove_modal_spinner()
            }
          )
          
          
        } # EO ifelse
      }) # EO observe
      
      
      
      #********************************************************************#
      # UI: GSEA output
      #********************************************************************#
      observe({
        if (!is.null(rv$GSEA_data)){
          output$UI_output_GSEA_microarray_raw <- renderUI({
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
                         h5("Here you can view the statistics from the GSEA. 
                            Click on the table to get the statistics per gene."),
                         hr(),
                         
                         # Statistics table
                         DT::dataTableOutput(outputId = "GSEA_table_microarray_raw") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1"),
                         
                         # Download button
                         downloadButton("download_GSEA_table_microarray_raw",
                                        "Download"),
                         br(),
                         
                         # Title + description of gene table
                         htmlOutput("text_GSEAgene_table_microarray_raw"),
                         h5("Here you can view the differential expression results for the selected gene set."),
                         hr(),
                         
                         # Gene table
                         DT::dataTableOutput(outputId = "GSEAgene_table_microarray_raw") %>%
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
                         h5("The bar chart visualizes the results from the Gene Set Enrichment Analysis."),
                         hr(),
                         actionButton("download_GSEAplot_microarray_raw", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         tags$div(class = "dropdown",
                                  tags$button(class = "dropbtn", icon("cog")),
                                  tags$div(class = "dropdown-content",
                                           
                                           # Color gradient
                                           tags$h4("Color gradient"),
                                           colourpicker::colourInput("lowcol_GSEAplot_microarray_raw",
                                                                     NULL,
                                                                     "#000072"),
                                           colourpicker::colourInput("midcol_GSEAplot_microarray_raw",
                                                                     NULL,
                                                                     "#FEE6CE"),
                                           colourpicker::colourInput("highcol_GSEAplot_microarray_raw",
                                                                     NULL,
                                                                     "red"),
                                           br(),
                                           
                                           # Number of gene sets
                                           tags$h4("# Gene sets"),
                                           sliderInput(
                                             inputId = "nSets_GSEAplot_microarray_raw",
                                             label = NULL,
                                             value = 10,
                                             min = 5,
                                             max = 20,
                                             step = 1),
                                           br(),br()
                                  )), # EO dropdownButton
                         
                         # Interactive plot output
                         plotly::plotlyOutput("GSEAplot_microarray_raw") %>%
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
                         actionButton("download_GSEAnetwork_microarray_raw", 
                                      "Download figure",
                                      icon = shiny::icon("download")),
                         br(),
                         br(),
                         
                         # Dropdown Button to adjust the plot settings
                         tags$div(class = "dropdown",
                                  tags$button(class = "dropbtn", icon("cog")),
                                  tags$div(class = "dropdown-content",
                                           
                                           # Color gradient
                                           tags$h4("Color gradient"),
                                           colourpicker::colourInput("lowcol_GSEAnetwork_microarray_raw",
                                                                     NULL,
                                                                     "#000072"),
                                           colourpicker::colourInput("midcol_GSEAnetwork_microarray_raw",
                                                                     NULL,
                                                                     "#FEE6CE"),
                                           colourpicker::colourInput("highcol_GSEAnetwork_microarray_raw",
                                                                     NULL,
                                                                     "red"),
                                           br(),
                                           
                                           # Network layout
                                           tags$h4("Network layout"),
                                           selectInput(inputId = "layout_GSEAnetwork_microarray_raw",
                                                       label = NULL,
                                                       choices = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds',
                                                                   'randomly', 'fr', 'kk', 'drl', 'lgl'),
                                                       selected = 'graphopt',
                                                       multiple = FALSE),
                                           br(),
                                           
                                           # Number of gene sets
                                           tags$h4("# Gene sets"),
                                           sliderInput(
                                             inputId = "nSets_GSEAnetwork_microarray_raw",
                                             label = NULL,
                                             value = 10,
                                             min = 5,
                                             max = 20,
                                             step = 1),
                                           br(),br()
                                  )), # EO dropdownButton
                         
                         # Make plot
                         plotOutput("GSEAnetwork_microarray_raw") %>%
                           shinycssloaders::withSpinner(color="#0dc5c1")
                ),
                
                #--------------------------------------------------------------#
                # Settings tab
                #--------------------------------------------------------------#
                tabPanel("Settings overview",
                         icon = icon("fas fa-file"),
                         h3(strong("GSEA settings")),
                         h5("To enhance reproducibility, download the overview of chosen GSEA settings."),
                         hr(),
                         DT::dataTableOutput(outputId = "GSEASettings_microarray_raw") %>% 
                           withSpinner(color="#0dc5c1"),
                         br(),
                         downloadButton("downloadGSEASettings_microarray_raw", 
                                        "Download table"),
                         downloadButton("downloadSessionInfo_microarray_raw", 
                                        "Session info")
                         
                )
                
              ) # EO tabSetPanel
            ) # EO tagList
          }) # EO renderUI
        }# EO if !is.null(rv$GSEA_data)
        
        
        # Allow user to download GSEA report
        output$UI_GSEAreport_microarray_raw <- renderUI({
          req(rv$GSEA_data)
          tagList(
            shinyWidgets::downloadBttn(outputId = "GSEAreport_microarray_raw",
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