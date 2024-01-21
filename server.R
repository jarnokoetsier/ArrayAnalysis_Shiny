#==============================================================================#
# Name: server.R
# Description: server of the ArrayAnalysis Shiny app
#==============================================================================#

# Increase connection size to allow for bigger data uploads
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

# Start server
server <- function(input, output, session){
  
  ##############################################################################
  
  # Preparations
  
  ##############################################################################
  
  # Set options for data upload
  options(shiny.maxRequestSize=125*1024^10)
  
  # Example meta data file
  output$downloadmeta_example <- downloadHandler(
    filename = "MetaData_example.csv",
    content = function(file){
      write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
    }
  )
  # Example meta data file
  output$downloadmeta_example_norm <- downloadHandler(
    filename = "MetaData_example.csv",
    content = function(file){
      write.csv(exampleMeta, file, quote = FALSE, row.names = FALSE)
    }
  )
  # Make list for reactive values
  rv <- reactiveValues()
  
  # Hide tabs (microarray - raw)
  hideTab("navbar", target = "panel_upload_microarray_raw")
  hideTab("navbar", target = "panel_preprocessing_microarray_raw")
  hideTab("navbar", target = "panel_statistics_microarray_raw" )
  hideTab("navbar", target = "panel_ORA_microarray_raw" )
  
  # Hide tabs (microarray - processed/normalized)
  hideTab("navbar", target = "panel_upload_microarray_norm")
  hideTab("navbar", target = "panel_preprocessing_microarray_norm")
  hideTab("navbar", target = "panel_statistics_microarray_norm" )
  hideTab("navbar", target = "panel_ORA_microarray_norm" )
  
  observe({
    
    ############################################################################
    
    # Home page
    
    ############################################################################
    observeEvent(input$startAnalysis,{
    # Analyse microarray data or RNA-seq data
    microarray_or_rnaseq <- reactive({
      req(input$microarray_or_rnaseq)
      return(input$microarray_or_rnaseq)
    })
    
    # Analyse raw or pre-processed data
    raw_or_norm <- reactive({
      req(input$raw_or_norm)
      return(input$raw_or_norm)
    })
    
    if (microarray_or_rnaseq() == "Documentation"){
      observeEvent(input$startAnalysis,{
      updateNavbarPage(session, "navbar",
                       selected = "documentation")
      })
    }
    
    
    ############################################################################
    
    # RNA-seq Analysis (WORK IN PROGRESS)
    
    ############################################################################
    
    if (microarray_or_rnaseq() == "RNA-Seq"){
      
      
      #************************************************************************#
      # Raw RNA-seq data
      #************************************************************************#
      if (raw_or_norm() == "Raw data"){
        
        #======================================================================#
        # Data Upload
        #======================================================================#
        
        # Go to data upload tab
        observeEvent(input$startAnalysis,{
          # Show RNA-seq upload tab
          showTab("navbar", target = "panel_upload_rnaseq_raw")
          
          # Go to RNA-seq tab
          updateNavbarPage(session, "navbar",
                           selected = "panel_upload_rnaseq_raw")
        })
        
        # Select (raw) expression matrix
        
        # Select meta data
        
        # Get overlap in samples
        
        
        #======================================================================#
        # Data Pre-processing
        #======================================================================#
        
        
        #======================================================================#
        # QC
        #======================================================================#
        
        
        #======================================================================#
        # Statistical analysis
        #======================================================================#
        
        
      } # raw or norm
      
      
      
      #************************************************************************#
      # Preprocessed RNA-seq data
      #************************************************************************#
      if (raw_or_norm() == "Processed data"){
        
        #======================================================================#
        # Data Upload
        #======================================================================#
        
        # Go to data upload tab
        observeEvent(input$startAnalysis,{
          # Show RNA-seq upload tab
          showTab("navbar", target = "panel_upload_rnaseq_norm")
          
          # Go to RNA-seq tab
          updateNavbarPage(session, "navbar",
                           selected = "panel_upload_rnaseq_norm")
        })
        
        # Select (preprocessed) expression matrix
        
        # Select meta data
        
        # Get overlap in samples
        
        
        #======================================================================#
        # Data Pre-processing
        #======================================================================#
        
        
        #======================================================================#
        # QC
        #======================================================================#
        
        
        #======================================================================#
        # Statistical analysis
        #======================================================================#
        
        
        
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
          
          # Remove the other microarray (raw) tabs
          hideTab("navbar", target = "panel_upload_microarray_norm")
          hideTab("navbar", target = "panel_preprocessing_microarray_norm")
          hideTab("navbar", target = "panel_statistics_microarray_norm" )
          hideTab("navbar", target = "panel_ORA_microarray_norm" )
          
          # Go to microarray (raw) upload tab
          updateNavbarPage(session, "navbar",
                           selected = "panel_upload_microarray_raw")
        })
        
        #----------------------------------------------------------------------#
        # Upload files
        #----------------------------------------------------------------------#
        
        # Observe "upload" input
        observeEvent(input$upload_microarray_raw,{
          
          # Show loading modal
          show_modal_spinner(text = "Reading data...",
                             color="#0dc5c1")
          
          # Select (raw) CEL files
          rv$celfiles <- getCELs(zippath = input$uploadCEL_microarray_raw)
          
          # Get metadata
          req(rv$celfiles) # CEL files are required for the metadata table
          if(length(rv$celfiles)>0){
            
            # Metadata in tsv or csv format
            if (input$MetaFileType==".tsv/.csv file"){
              rv$metaData <- getMetaData(path = input$uploadMeta_microarray_raw_tsv$datapath,
                                         celfiles = rv$celfiles,
                                         filetype = input$MetaFileType)
            }
            
            # Metadata in Series Matrix File format
            if (input$MetaFileType=="Series Matrix File"){
              rv$metaData <- getMetaData(path = input$uploadMeta_microarray_raw_smf$datapath,
                                         celfiles = rv$celfiles,
                                         filetype = input$MetaFileType)
            }
          } else {
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No expression data",
              type = "error")
          }
          
          # Read raw expression data
          req(rv$metaData)
          if(nrow(rv$metaData)>0){
            
            # Read data
            rv$celfiles_fil <- rv$celfiles[str_remove(basename(rv$celfiles),"\\.CEL.*") %in% rownames(rv$metaData)]
            rv$gxData <- readCELs(rv$celfiles_fil)
            
            # Check whether all expression samples have meta data available
            if (nrow(rv$metaData) < length(rv$celfiles)){
              sendSweetAlert(
                session = session,
                title = "Warning!",
                text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                type = "warning")
            }
          
            
            
            #------------------------------------------------------------------#
            # Outputs
            #------------------------------------------------------------------#
            
            # Set tables to zero
            output$exprTable_microarray_raw <- renderDataTable(NULL)
            output$metaTable_microarray_raw <- renderDataTable(NULL)
            
            # Print expression table
            output$exprTable_microarray_raw <- renderDataTable({
              req(rv$gxData)
              output <- head(exprs(rv$gxData),6)
              colnames(output) <- str_remove(colnames(output),"\\.CEL.*")
              return(output)
              
            }, options = list(pageLength = 6))
            
            # Print meta table
            output$metaTable_microarray_raw <- renderDataTable({
              req(rv$metaData)
              return(rv$metaData)
            }, options = list(pageLength = 6))
            
            # Render UI for main tab
            output$UI_upload_microarray_raw <- renderUI({
              tagList(
                tabsetPanel(
                  tabPanel("Expression matrix",
                           h3(strong("Expression matrix")),
                           h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "exprTable_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1")),
                  tabPanel("Meta data",                  # Meta table
                           h3(strong("Meta data")),
                           h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "metaTable_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1")),
                )
              )
            })
            
            # Allow user to go to next tab
            output$next_upload_microarray_raw <- renderUI({
              req(rv$metaData)
              req(rv$gxData)
              
              # Remove modal
              remove_modal_spinner()
              
              # Show RNA-seq upload tab
              showTab("navbar", target = "panel_preprocessing_microarray_raw")
              
              # Show message
              if (nrow(rv$metaData) >= length(rv$celfiles)){
                sendSweetAlert(
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
                actionBttn(inputId = "next_upload_microarray_raw",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          } else {
            # No common samples between meta data and expression data
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No common samples",
              type = "error")
            
            # Remove loading modal
            remove_modal_spinner()
          }
          
        }) # EO observeEvent
        
        #----------------------------------------------------------------------#
        # Run Example
        #----------------------------------------------------------------------#
        
        # Observe "upload" input
        observeEvent(input$example_microarray_raw,{
          
          # Show loading modal
          show_modal_spinner(text = "Reading data...",
                             color="#0dc5c1")
          
          # Select (raw) CEL files
          rv$celfiles <- getCELs(zippath = "Data/GSE6955_RAW.zip",
                                 shiny_upload = FALSE)
          
          # Get metadata
          req(rv$celfiles) # CEL files are required for the metadata table
          if(length(rv$celfiles)>0){
            
            # Metadata in tsv or csv format
            if (input$MetaFileType==".tsv/.csv file"){
              rv$metaData <- getMetaData(path = "Data/metaData_GSE6955.csv",
                                         celfiles = rv$celfiles,
                                         filetype = ".tsv/.csv file")
            }
          } else {
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No expression data",
              type = "error")
          }
          
          # Read raw expression data
          req(rv$metaData)
          if(nrow(rv$metaData)>0){
            
            # Read data
            rv$celfiles_fil <- rv$celfiles[str_remove(basename(rv$celfiles),"\\.CEL.*") %in% rownames(rv$metaData)]
            rv$gxData <- readCELs(rv$celfiles_fil)
            
            # Check whether all expression samples have meta data available
            if (nrow(rv$metaData) < length(rv$celfiles)){
              sendSweetAlert(
                session = session,
                title = "Warning!",
                text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                type = "warning")
            }
            
            
            
            #------------------------------------------------------------------#
            # Outputs
            #------------------------------------------------------------------#
            
            # Set tables to zero
            output$exprTable_microarray_raw <- renderDataTable(NULL)
            output$metaTable_microarray_raw <- renderDataTable(NULL)
            
            # Print expression table
            output$exprTable_microarray_raw <- renderDataTable({
              req(rv$gxData)
              output <- head(exprs(rv$gxData),6)
              colnames(output) <- str_remove(colnames(output),"\\.CEL.*")
              return(output)
              
            }, options = list(pageLength = 6))
            
            # Print meta table
            output$metaTable_microarray_raw <- renderDataTable({
              req(rv$metaData)
              return(rv$metaData)
            }, options = list(pageLength = 6))
            
            # Render UI for main tab
            output$UI_upload_microarray_raw <- renderUI({
              tagList(
                tabsetPanel(
                  tabPanel("Expression matrix",
                           h3(strong("Expression matrix")),
                           h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "exprTable_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1")),
                  tabPanel("Meta data",                  # Meta table
                           h3(strong("Meta data")),
                           h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "metaTable_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1")),
                )
              )
            })
            
            # Allow user to go to next tab
            output$next_upload_microarray_raw <- renderUI({
              req(rv$metaData)
              req(rv$gxData)
              
              # Remove modal
              remove_modal_spinner()
              
              # Show RNA-seq upload tab
              showTab("navbar", target = "panel_preprocessing_microarray_raw")
              
              # Show message
              if (nrow(rv$metaData) >= length(rv$celfiles)){
                sendSweetAlert(
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
                actionBttn(inputId = "next_upload_microarray_raw",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          } else {
            # No common samples between meta data and expression data
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No common samples",
              type = "error")
            
            # Remove loading modal
            remove_modal_spinner()
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
        
        
        # 1. Select outliers
        output$UI_outlier_microarray_raw <- renderUI({
          
          # If outliers are checked, show possible sample to select as outliers
          if(!input$outlier){
            samples <- rownames(rv$metaData)
            pickerInput(inputId = "select_outliers",
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
        output$UI_groupselect_microarray_raw <- renderUI({
            pickerInput(inputId = "groupselect",
                        label = NULL,
                        choices = colnames(rv$metaData),
                        selected = autoGroup(rv$metaData),
                        multiple = TRUE)
        })
        
        # print the experimental levels
        output$experimentallevels <- renderText({
          req(input$groupselect)
          
          if(length(input$groupselect) > 1){
            # If more than one metadata column is selected: concatenate  the vectors
            experimentFactor <- apply(rv$metaData[,input$groupselect], 1, paste, collapse = "_" )
          } else{
            # If only one metadata column is selected: nothing to do!
            experimentFactor <- rv$metaData[,input$groupselect]
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
        
        # Get organism
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
        

        # 3. Pre-process the raw data
        observeEvent(input$start_preprocessing_microarray_raw,{
          
          rv$annotations <- input$annotations_microarray 
          if (rv$annotations == "Custom annotations"){
            rv$ProbeAnnotation <- input$CDFtype_microarray 
            rv$Organism <- input$species_microarray_raw 
          } else{
            rv$ProbeAnnotation <- NULL
          }
          
          # Show modal
          show_modal_spinner(text = "Pre-processing data...",
                             color="#0dc5c1")
          
          # Select outlier
          if (!isTRUE(input$outier)){
            rv$outlier <- input$select_outliers
          } else{
            rv$outlier <- NULL
          }
          
          # Remove outlier
          rv$celfiles_sel <- rv$celfiles_fil
          if (!is.null(rv$outlier)){
            for (i in 1:length(rv$outlier)){
              rv$celfiles_sel <- rv$celfiles_sel[!str_detect(rv$celfiles_sel, rv$outlier[i])]
            }
          }
          
          # Filter metadata and expression data (samples in correct order)
          rv$gxData_fil <- readCELs(rv$celfiles_sel)
          rv$metaData_fil <- rv$metaData[str_remove(colnames(exprs(rv$gxData_fil)),"\\.CEL.*"),]
          
          # Experiment factor
          if(length(input$groupselect) > 1){
            rv$experimentFactor <- factor(apply(rv$metaData_fil[,input$groupselect], 1, paste, collapse = "_" ))
          } else{
            rv$experimentFactor <- factor(rv$metaData_fil[,input$groupselect])
          }
          
          # Normalization
          if (class(rv$gxData_fil) == "GeneFeatureSet"){
            sendSweetAlert(
              session = session,
              title = "Warning",
              text = "The data is from the Gene ST 2.x series of arrays. 
              Currently only RMA normalization without custom probeset annotation 
              is supported. So, the pre-processing will be performed accordingly.",
              type = "info")
            rv$annotations <- "No annotations"
            rv$ProbeAnnotation <- NULL 
          }
          rv$normData <- microarrayNormalization(gxData = rv$gxData_fil,
                                                 experimentFactor = rv$experimentFactor,
                                                 normMeth = input$normMeth_microarray,
                                                 CDFtype = input$CDFtype_microarray,
                                                 species = input$species_microarray_raw,
                                                 annotations = rv$annotations,
                                                 perGroup_name = input$perGroup_microarray,
                                                 annot_file_datapath = input$annot_file_microarray$datapath)
        
          normMatrix <- exprs(rv$normData)
          
          if (rv$annotations == "Custom annotations"){
            id_names <- as.character(str_remove(rownames(normMatrix),"_.*"))
            normMatrix <- normMatrix[!duplicated(id_names),]
            rownames(normMatrix) <- id_names[!duplicated(id_names)]
          }
          rv$normMatrix <- normMatrix
          rm(normMatrix)
          
          
        #======================================================================#
        # QC
        #======================================================================#
          
          #********************************************************************#
          # Expression values
          #********************************************************************#
          
          # Print expression table
          output$exprTable_microarray_norm <- renderDataTable({
            
            # Remove modal
            remove_modal_spinner()
            
            # Show message
            sendSweetAlert(
              session = session,
              title = "Info",
              text = "The data has been pre-processed. Please check the different 
              QC plots on this page to assess the pre-processing quality.",
              type = "info")
            
            # Show microarray statistics tab
            showTab("navbar", target = "panel_statistics_microarray_raw")
            
            output <- rv$normMatrix
            
            if (!is.null(rv$ProbeAnnotation)){
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
          
          # Boxplot of single gene
          output$ExprBoxplot_microarray <- renderPlotly({
            req(input$exprTable_microarray_norm_rows_selected)
            
            # Set colors
            myPalette <- colorsByFactor(rv$experimentFactor)
            legendColors <- myPalette$legendColors
            
            # Prepare expression data frame
            plotExpr <- data.frame(
              logExpr = rv$normMatrix[input$exprTable_microarray_norm_rows_selected,],
              Grouping = rv$experimentFactor
            )
            
            # make interactive plot
            plotExpr %>%
              plot_ly(x = ~Grouping,y = ~as.numeric(logExpr),
                      color = ~Grouping, colors = legendColors, type = "box") %>%
              layout(xaxis = list(title = " "),
                     yaxis = list(title = 'log intensity'),
                     legend = list(title=list(text='Group')),
                     showlegend = FALSE)
            
          })
          
          #********************************************************************#
          # Boxplots
          #********************************************************************#
          
          # Boxplots of all genes together
          output$boxplots_microarray_norm <- renderImage({
            getBoxplots(experimentFactor = rv$experimentFactor,
                        normData = rv$normData)
          },deleteFile = TRUE)
          
          #********************************************************************#
          # Densityplots
          #********************************************************************#
          
          # Densityplots of all genes together
          output$densityplots_microarray_norm <- renderPlotly({
            getDensityplots(experimentFactor = rv$experimentFactor,
                            normMatrix = rv$normMatrix)
            
          })
          
          #********************************************************************#
          # Heatmap
          #********************************************************************#
          
          output$heatmap_microarray_norm  <- renderPlotly({
                myPalette <- colorsByFactor(rv$experimentFactor)

                # Plot colors
                plotColors <- myPalette$plotColors

                # Legend colors
                legendColors <- myPalette$legendColors

                # Set cluster options
                clusterOption1 <- input$clusteroption1_microarray_raw
                clusterOption2 <- input$clusteroption2_microarray_raw

                #note: for computing array correlation, euclidean would not make sense
                #only use euclidean distance to compute the similarity of the correlation
                #vectors for the arrays
                COpt1 <- "pearson"
                if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
                crp <- cor(rv$normMatrix, use="complete.obs", method=COpt1)

                switch(tolower(clusterOption1),
                       "pearson" = {
                         my.dist <- function(x) cor.dist(x, abs=FALSE)
                       },
                       "spearman" = {
                         my.dist <- function(x) spearman.dist(x, abs=FALSE)
                       },
                       "euclidean" = {
                         my.dist <- function(x) euc(x)
                       }
                )

                # Perform clustering
                my.hclust <- function(d) hclust(d, method=clusterOption2)
                names(legendColors) <- levels(rv$experimentFactor)
                sidecolors <- data.frame(rv$experimentFactor)
                colnames(sidecolors) <- "Experimental group"

                if (input$heatmaptheme_microarray_raw == "Default"){
                 gradient = viridis(n = 256, alpha = 1, begin = 0, end = 1,
                                    option = "viridis")
                }

                if (input$heatmaptheme_microarray_raw == "Yellow-red"){
                  gradient = heat.colors(100)
                }

                if (input$heatmaptheme_microarray_raw == "Dark"){
                  gradient = RColorBrewer::brewer.pal(8, "Dark2")
                }

                # Make heatmap
                p <- heatmaply(crp, plot_method = "plotly", distfun = my.dist,
                          hclustfun = my.hclust, symm = TRUE, seriate = "mean",
                          titleX = FALSE, titleY = FALSE, key.title = NULL,
                          show_dendrogram = c(TRUE, FALSE), col_side_colors = sidecolors,
                          col_side_palette = legendColors, column_text_angle = 90,
                          colors = gradient)
                
                return(p)

             })
          
          #********************************************************************#
          # PCA
          #********************************************************************#
          
          #Perform PCA
          rv$PCA_data <- prcomp(t(rv$normMatrix[apply(rv$normMatrix, 1, var) != 0,]),
                                retx = TRUE, 
                                center = TRUE,
                                scale.=TRUE)
          
          
          # Make PCA plot
          output$PCA_microarray_raw <- renderPlotly({
            
            if(length(input$colorFactor_microarray_raw) > 1){
              colorFactor <- factor(apply(rv$metaData_fil[,input$colorFactor_microarray_raw], 1, paste, collapse = "_" ))
            } else{
              colorFactor <- factor(rv$metaData_fil[,input$colorFactor_microarray_raw])
            }
            print(head(colorFactor))
            
            plot_PCA(PC_data = rv$PCA_data, 
                     colorFactor = colorFactor, 
                     xpc = as.numeric(str_remove(input$xpca_microarray_raw,"PC")), 
                     ypc = as.numeric(str_remove(input$ypca_microarray_raw,"PC")), 
                     zpc = ifelse(input$xyz_microarray_raw,as.numeric(str_remove(input$zpca_microarray_raw,"PC")),3), 
                     xyz = input$xyz_microarray_raw)
            
          })
          #********************************************************************#
          # UI
          #********************************************************************#
          
          # Render UI for main tab
          output$UI_QC_microarray_norm <- renderUI({
            tagList(
              tabsetPanel(
                tabPanel("Expression values",
                         icon = icon("fas fa-mouse-pointer"),
                         h3(strong("Normalized expression values")),
                         h5("Here you can view the normalized log intensity (expression) values."),
                         hr(),
                         dataTableOutput(outputId = "exprTable_microarray_norm") %>% 
                           withSpinner(color="#0dc5c1"),
                         downloadButton("downloadNormalizedData_microarray_raw", 
                                        "Download"),
                         plotlyOutput("ExprBoxplot_microarray")%>% 
                           withSpinner(color="#0dc5c1")
                ),
                tabPanel("Boxplots",
                         icon = icon("fas fa-file"),
                         h3(strong("Boxplots")),
                         h5("These are boxplots of the normalized expression expression values."),
                         hr(),
                         plotOutput(outputId = "boxplots_microarray_norm") %>% 
                           withSpinner(color="#0dc5c1")
                ),
                
                tabPanel("Density plots",
                         icon = icon("fas fa-mouse-pointer"),
                         h3(strong("Density plots")),
                         h5("These are density plots of the normalized expression expression values."),
                         hr(),
                         plotlyOutput(outputId = "densityplots_microarray_norm") %>% 
                           withSpinner(color="#0dc5c1")
                ),
                
                tabPanel("Correlation plot", 
                         icon = icon("fas fa-mouse-pointer"),
                         
                         br(),
                         
                         fluidRow(
                           
                           #Distance method
                           column(4,
                                  pickerInput(
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
                                  pickerInput(
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
                           
                           selectInput(inputId = 'heatmaptheme_microarray_raw',
                                       label = NULL,
                                       choices = c("Default", 
                                                   "Yellow-red", 
                                                   "Dark")),
                           
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"), width = "300px",
                           
                           tooltip = tooltipOptions(
                             title = "Click to change theme!")
                           
                         ),
                         
                         plotlyOutput("heatmap_microarray_norm", 
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
                                  pickerInput(inputId = "colorFactor_microarray_raw",
                                              label = "Color by:",
                                              choices = colnames(rv$metaData_fil),
                                              selected = autoGroup(rv$metaData_fil),
                                              multiple = TRUE)
                                  ),
                           column(3,
                                  br(),
                                  # 3D plot?
                                  materialSwitch(
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
                           plotlyOutput("PCA_microarray_raw")%>% 
                             withSpinner(color="#0dc5c1")
                         ) 

                         ) # PCA tabpanel
              ) # tabsetpanel
            ) # taglist
          }) # renderUI
          
          # Allow user to go to next tab
          output$UI_next_preprocessing_microarray_raw <- renderUI({
            req(rv$normData)
            tagList(
              hr(),
              h2(strong("Continue your analysis")),
              actionBttn(inputId = "next_preprocessing_microarray_raw",
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
        observeEvent(input$next_preprocessing_microarray_raw,{
          
          # Go to microarray statistics tab
          updateNavbarPage(session, "navbar",
                           selected = "panel_statistics_microarray_raw")
        })
        
        # SIDEPANEL:
        
        # Select experimental group
        output$UI_expFactor_microarray_raw <- renderUI({
          tagList(
          pickerInput(inputId = "expFactor_microarray_raw",
                      label = "Experiment factor:",
                      choices = colnames(rv$metaData_fil),
                      selected = autoGroup(rv$metaData_fil),
                      multiple = TRUE)
          )
        })
        
        # print the experimental levels
        output$expFactor_levels_micorarray_raw <- renderText({
          req(input$expFactor_microarray_raw)
          
          if(length(input$expFactor_microarray_raw) > 1){
            experimentFactor <- apply(rv$metaData_fil[,input$expFactor_microarray_raw], 1, paste, collapse = "_" )
          } else{
            experimentFactor <- rv$metaData_fil[,input$expFactor_microarray_raw]
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
        
        # Select continuous covariates
        output$UI_covGroups_num_microarray_raw <- renderUI({
          req(input$expFactor_microarray_raw)
          tagList(
            pickerInput(inputId = "covGroups_num_microarray_raw",
                        label = "Continuous covariates (e.g., age):",
                        choices = setdiff(colnames(rv$metaData_fil),
                                          input$expFactor_microarray_raw),
                        selected = NULL,
                        multiple = TRUE)
          )
        })
        
        # Select discrete covariates
        output$UI_covGroups_char_microarray_raw <- renderUI({
          req(input$expFactor_microarray_raw)
          tagList(
            pickerInput(inputId = "covGroups_char_microarray_raw",
                        label = "Discrete covariates (e.g., sex):",
                        choices = setdiff(colnames(rv$metaData_fil),
                                          input$expFactor_microarray_raw),
                        selected = NULL,
                        multiple = TRUE)
          )
        })
        
        # Select comparisons
        output$UI_comparisons_microarray_raw <- renderUI({
          req(input$expFactor_microarray_raw)
          tagList(
            multiInput(
              inputId = "comparisons_microarray_raw",
              label = "Comparisons:", 
              choices = makeComparisons(make.names(unique(rv$metaData[,input$expFactor_microarray_raw]))),
              selected = makeComparisons(make.names(unique(rv$metaData[,input$expFactor_microarray_raw])))[1]
            )
          )
        })
        
        
        # Select comparisons
        output$UI_biomart_dataset_microarray_raw <- renderUI({
          req(input$addAnnotation_microarray_raw)
          pickerInput(inputId = "biomart_dataset_microarray_raw",
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
            
            pickerInput(inputId = "biomart_filter_microarray_raw",
                        label = "Probeset ID",
                        choices = filterList[[input$biomart_dataset_microarray_raw]],
                        selected = selFilter(rv$ProbeAnnotation),
                        multiple = FALSE),
            
            pickerInput(inputId = "biomart_attributes_microarray_raw",
                        label = "Output",
                        choices = c("ensembl_gene_id",
                                    "entrezgene_id",
                                    "gene_name"),
                        selected = "gene_name",
                        multiple = TRUE)
          )

        })
        
      

        #======================================================================#
        # Output of statistical analysis
        #======================================================================#
        observeEvent(input$calculate_statistics_microarray_raw,{
          show_modal_spinner(text = "Statistical analysis...",
                             color="#0dc5c1")
          
          # Calculate statistics
          if (isTRUE(input$addAnnotation_microarray_raw)){
            rv$top_table <- getStatistics(normMatrix = rv$normMatrix, 
                                          metaData = rv$metaData_fil, 
                                          expFactor = input$expFactor_microarray_raw, 
                                          covGroups_num = input$covGroups_num_microarray_raw,
                                          covGroups_char = input$covGroups_char_microarray_raw,
                                          comparisons = input$comparisons_microarray_raw,
                                          addAnnotation = input$addAnnotation_microarray_raw,
                                          biomart_dataset = input$biomart_dataset_microarray_raw,
                                          biomart_attributes = unique(c(input$biomart_filter_microarray_raw,
                                                                        input$biomart_attributes_microarray_raw)),
                                          biomart_filters = input$biomart_filter_microarray_raw)
            
          } else{
            rv$top_table <- getStatistics(normMatrix = rv$normMatrix, 
                                          metaData = rv$metaData_fil, 
                                          expFactor = input$expFactor_microarray_raw, 
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
            #req(rv$top_table)
            
            if (!is.null(rv$top_table)){
              # Remove modal
              remove_modal_spinner()
              
              # Show comparisons
              output$UI_comparisons_view_microarray_raw <- renderUI({
                pickerInput(inputId = "comparisons_view_microarray_raw",
                            label = "Select comparison:",
                            choices = names(rv$top_table),
                            selected = names(rv$top_table)[1],
                            multiple = FALSE)
              })
              
              # Show message
              sendSweetAlert(
                session = session,
                title = "Info",
                text = "Statistical analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
                type = "info")
              
              # Show microarray statistics tab
              showTab("navbar", target = "panel_ORA_microarray_raw")
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
          output$top_table_microarray_raw <- renderDataTable({
            req(input$comparisons_view_microarray_raw)
            
            output <- rv$top_table[[input$comparisons_view_microarray_raw]]
            
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
          
          # Boxplot of single gene
          output$ExprBoxplot_statistics_microarray_raw <- renderPlotly({
            
            req(input$top_table_microarray_raw_rows_selected)
            req(input$comparisons_view_microarray_raw)
            
            # Set colors
            myPalette <- colorsByFactor(rv$experimentFactor)
            legendColors <- myPalette$legendColors
            
            # Prepare expression data frame
            gene <- rv$top_table[[input$comparisons_view_microarray_raw]]$GeneID[input$top_table_microarray_raw_rows_selected]
            plotExpr <- data.frame(
              logExpr = rv$normMatrix[as.character(rownames(rv$normMatrix)) %in% as.character(gene),],
              Grouping = factor(rv$metaData_fil[,input$expFactor_microarray_raw])
            )
            
            # make interactive plot
            plotExpr %>%
              plot_ly(x = ~Grouping,y = ~as.numeric(logExpr),
                      color = ~Grouping, colors = legendColors, type = "box") %>%
              layout(xaxis = list(title = " "),
                     yaxis = list(title = 'log intensity'),
                     legend = list(title=list(text='Group')),
                     showlegend = FALSE)
            
          })
          
          #********************************************************************#
          # Histograms
          #********************************************************************#
          output$Phistogram_microarray_raw <- renderPlotly({
            req(rv$top_table)
            p <- makePHistogram(rv$top_table[[input$comparisons_view_microarray_raw]][,"p-value"])
            return(p)
          })
          
          output$logFChistogram_microarray_raw <- renderPlotly({
            req(rv$top_table)
            p <- makelogFCHistogram(rv$top_table[[input$comparisons_view_microarray_raw]][,"log2FC"])
            return(p)
          })
          
          #********************************************************************#
          # Volcano plot
          #********************************************************************#
          observeEvent(input$plot_volcano_microarray_raw, {
            req(rv$top_table)
            req(input$rawp_volcano_microarray_raw)
            req(input$p_thres_volcano_microarray_raw)
            req(input$logFC_thres_volcano_microarray_raw)
            p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_raw]], 
                             p = input$rawp_volcano_microarray_raw, 
                             p_threshold = input$p_thres_volcano_microarray_raw, 
                             logFC_threshold = input$logFC_thres_volcano_microarray_raw)
            
            output$volcano_microarray_raw <- renderPlotly(p)
          }, ignoreNULL = FALSE)
          
          
          # UI: Output in different tabs
          
          observe({
            output$UI_output_statistics_microarray_raw <- renderUI({
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
                           dataTableOutput(outputId = "top_table_microarray_raw") %>% 
                             withSpinner(color="#0dc5c1"),
                           downloadButton("download_top_table_microarray_raw", 
                                          "Download"),
                           plotlyOutput("ExprBoxplot_statistics_microarray_raw")%>% 
                             withSpinner(color="#0dc5c1")
                  ),
                  #********************************************************************#
                  # histogram tab
                  #********************************************************************#
                  
                  tabPanel("Histograms",
                           icon = icon("fas fa-mouse-pointer"),
                           br(),
                           plotlyOutput("Phistogram_microarray_raw")%>% 
                             withSpinner(color="#0dc5c1"),
                           hr(),
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
                                    #P value threshold
                                    numericInput(
                                      inputId = "p_thres_volcano_microarray_raw",
                                      label = "P threshold",
                                      value = 0.05)
                                    ),
                             column(3,
                                    #logFC threshold
                                    numericInput(
                                      inputId = "logFC_thres_volcano_microarray_raw",
                                      label = "logFC threshold",
                                      value = 1)
                                    )
                           ),
                           hr(),
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
            pickerInput(inputId = "comparisons_view_ORA_microarray_raw",
                        label = NULL,
                        choices = names(rv$top_table),
                        selected = names(rv$top_table)[1],
                        multiple = FALSE)
          })
        })
        
        # Which columns contain gene ids
        observe({
          req(input$comparisons_view_ORA_microarray_raw)
          col_choice <- 1
          if (length(colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])) > 6){
            col_choice <- c(1,7:ncol(rv$top_table[[input$comparisons_view_ORA_microarray_raw]]))
          }
          output$UI_geneID_ORA_microarray_raw <- renderUI({
            tagList(
              selectInput(inputId = "organism_ORA_microarray_raw",
                          label = "Organism",
                          choices = c("Bos taurus",
                                      "Caenorhabditis elegans",
                                      "Homo sapiens",
                                      "Mus musculus", 
                                      "Rattus norvegicus"),
                          selected = rv$Organism),
              
              pickerInput(inputId = "geneID_ORA_microarray_raw",
                          label = "Which column of the top table contains the gene IDs?",
                          choices = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])[col_choice],
                          selected = colnames(rv$top_table[[input$comparisons_view_ORA_microarray_raw]])[1],
                          multiple = FALSE),
              
              pickerInput(inputId = "selID_ORA_microarray_raw",
                          label = "Which gene ID to use?",
                          choices = c("ensembl_gene_id" = "ENSEMBL", 
                                      "entrezgene_id" = "ENTREZID", 
                                      "gene_name" = "SYMBOL"),
                          selected = selID(rv$ProbeAnnotation),
                          multiple = FALSE),
            )
          })
        })
        
        # Show modal
        observeEvent(input$calculate_ORA_microarray_raw,{
          show_modal_spinner(text = "Overrepresentation analysis...",
                             color="#0dc5c1")
          
          # Perform ORA
          rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_raw]],
                               geneset = input$geneset_ORA_microarray_raw,
                            geneID_col = input$geneID_ORA_microarray_raw,
                            geneID_type = input$selID_ORA_microarray_raw,
                            organism = input$organism_ORA_microarray_raw,
                            rawadj = input$rawp_ORA_microarray_raw,
                            updown = input$updown_ORA_microarray_raw,
                            p_thres = input$p_thres_ORA_microarray_raw,
                            logFC_thres = input$logFC_thres_ORA_microarray_raw)
          
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
                
                if (input$geneset_ORA_microarray_raw != "WikiPathways"){
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
          # Show microarray upload tab
          showTab("navbar", target = "panel_upload_microarray_norm")
          
          hideTab("navbar", target = "panel_upload_microarray_raw")
          hideTab("navbar", target = "panel_preprocessing_microarray_raw")
          hideTab("navbar", target = "panel_statistics_microarray_raw" )
          hideTab("navbar", target = "panel_ORA_microarray_raw" )
          
          # Go to microarray tab
          updateNavbarPage(session, "navbar",
                           selected = "panel_upload_microarray_norm")
        })
        
        #----------------------------------------------------------------------#
        # Upload files
        #----------------------------------------------------------------------#
        
        observeEvent(input$upload_microarray_norm,{
          
          # Show modal
          show_modal_spinner(text = "Reading data...",
                             color="#0dc5c1")
          
          # Read expression data
          rv$gxData <- getGEO(filename=input$uploadExprData_microarray_norm_smf$datapath)
         
           # Guess the organism
          rv$Organism <- getOrganism(gxData = rv$gxData)
          
          # Get metadata
          if (input$MetaFileType_norm==".tsv/.csv file"){
            rv$metaData <- getMetaData(path = input$uploadMeta_microarray_norm_tsv$datapath,
                                       celfiles = colnames(exprs(rv$gxData)),
                                       filetype = input$MetaFileType_norm)
          }
          if (input$MetaFileType_norm=="Series Matrix File"){
            rv$metaData <- getMetaData(path = input$uploadMeta_microarray_norm_smf$datapath,
                                       celfiles =  colnames(exprs(rv$gxData)),
                                       filetype = input$MetaFileType_norm)
          }

          # Read raw expression data
          req(rv$metaData)
          if(nrow(rv$metaData)>0){
            
            # check if some samples are removed
            if (nrow(rv$metaData) != ncol(exprs(rv$gxData))){
              sendSweetAlert(
                session = session,
                title = "Warning!",
                text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                type = "warning")
            }
            
            # Filter expression data for samples with metadata
            rv$gxData <- ExpressionSet(assayData = exprs(rv$gxData)[,rownames(rv$metaData)])
            
            
            #------------------------------------------------------------------#
            # Outputs
            #------------------------------------------------------------------#
            
            # Set tables to zero
            output$exprTable_microarray_norm <- renderDataTable(NULL)
            output$metaTable_microarray_norm <- renderDataTable(NULL)
            
            # Print expression table
            output$exprTable_microarray_norm <- renderDataTable({
              req(rv$gxData)
              output <- head(exprs(rv$gxData),6)
              colnames(output) <- str_remove(colnames(output),"\\.CEL.*")
              return(output)
              
            }, options = list(pageLength = 6))
            
            # Print meta table
            output$metaTable_microarray_norm <- renderDataTable({
              req(rv$metaData)
              return(rv$metaData)
            }, options = list(pageLength = 6))
            
            # Render UI for main tab
            output$UI_upload_microarray_norm <- renderUI({
              tagList(
                tabsetPanel(
                  tabPanel("Expression matrix",
                           h3(strong("Expression matrix")),
                           h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "exprTable_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1")),
                  tabPanel("Meta data",                  # Meta table
                           h3(strong("Meta data")),
                           h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "metaTable_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1")),
                )
              )
            })
            
            # Allow user to go to next tab
            output$next_upload_microarray_norm <- renderUI({
              req(rv$metaData)
              req(rv$gxData)
              
              # Remove modal
              remove_modal_spinner()
              
              # Show RNA-seq upload tab
              showTab("navbar", target = "panel_preprocessing_microarray_norm")
              
              # Show message
              if (nrow(rv$metaData) >= ncol(exprs(rv$gxData))){
                sendSweetAlert(
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
                actionBttn(inputId = "next_upload_microarray_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          } else {
            # No common samples between meta data and expression data
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No common samples",
              type = "error")
            
            remove_modal_spinner()
          }
          
        }) # Observe event
   
        #----------------------------------------------------------------------#
        # Run Example
        #----------------------------------------------------------------------#
        
        observeEvent(input$example_microarray_norm,{
          
          # Show modal
          show_modal_spinner(text = "Reading data...",
                             color="#0dc5c1")
          
          # Read expression data
          rv$gxData <- getGEO(filename="Data/GSE6955_series_matrix.txt.gz")
          
          # Guess the organism
          rv$Organism <- getOrganism(gxData = rv$gxData)
          
          # Get metadata
          rv$metaData <- getMetaData(path = "Data/metaData_GSE6955.csv",
                                     celfiles = colnames(exprs(rv$gxData)),
                                     filetype = ".tsv/.csv file")
         

          # Read raw expression data
          req(rv$metaData)
          if(nrow(rv$metaData)>0){
            
            # check if some samples are removed
            if (nrow(rv$metaData) != ncol(exprs(rv$gxData))){
              sendSweetAlert(
                session = session,
                title = "Warning!",
                text = "One or more samples in the expression data file do not have meta
                  data available. These samples are excluded from the analysis.",
                type = "warning")
            }
            
            # Filter expression data for samples with metadata
            rv$gxData <- ExpressionSet(assayData = exprs(rv$gxData)[,rownames(rv$metaData)])
            
            
            #------------------------------------------------------------------#
            # Outputs
            #------------------------------------------------------------------#
            
            # Set tables to zero
            output$exprTable_microarray_norm <- renderDataTable(NULL)
            output$metaTable_microarray_norm <- renderDataTable(NULL)
            
            # Print expression table
            output$exprTable_microarray_norm <- renderDataTable({
              req(rv$gxData)
              output <- head(exprs(rv$gxData),6)
              colnames(output) <- str_remove(colnames(output),"\\.CEL.*")
              return(output)
              
            }, options = list(pageLength = 6))
            
            # Print meta table
            output$metaTable_microarray_norm <- renderDataTable({
              req(rv$metaData)
              return(rv$metaData)
            }, options = list(pageLength = 6))
            
            # Render UI for main tab
            output$UI_upload_microarray_norm <- renderUI({
              tagList(
                tabsetPanel(
                  tabPanel("Expression matrix",
                           h3(strong("Expression matrix")),
                           h5("These are the first six entries of the expression matrix. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "exprTable_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1")),
                  tabPanel("Meta data",                  # Meta table
                           h3(strong("Meta data")),
                           h5("This is a preview of the meta data. Please check if the data 
               has been correctly imported."),
                           hr(),
                           dataTableOutput(outputId = "metaTable_microarray_norm") %>% 
                             withSpinner(color="#0dc5c1")),
                )
              )
            })
            
            # Allow user to go to next tab
            output$next_upload_microarray_norm <- renderUI({
              req(rv$metaData)
              req(rv$gxData)
              
              # Remove modal
              remove_modal_spinner()
              
              # Show RNA-seq upload tab
              showTab("navbar", target = "panel_preprocessing_microarray_norm")
              
              # Show message
              if (nrow(rv$metaData) >= ncol(exprs(rv$gxData))){
                sendSweetAlert(
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
                actionBttn(inputId = "next_upload_microarray_norm",
                           label = "Next",
                           style = "jelly",
                           color = "danger",
                           icon = icon("arrow-right"))
              )
            })
            
          } else {
            # No common samples between meta data and expression data
            sendSweetAlert(
              session = session,
              title = "Error!!",
              text = "No common samples",
              type = "error")
            
            remove_modal_spinner()
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
            
            pickerInput(inputId = "select_outliers_norm",
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
          pickerInput(inputId = "groupselect_norm",
                      label = NULL,
                      choices = colnames(rv$metaData),
                      selected = autoGroup(rv$metaData),
                      multiple = TRUE)
        })
        
        # print the experimental levels
        output$experimentallevels_norm <- renderText({
          req(input$groupselect_norm)
          
          if(length(input$groupselect) > 1){
            experimentFactor <- apply(rv$metaData[,input$groupselect_norm], 1, paste, collapse = "_" )
          } else{
            experimentFactor <- rv$metaData[,input$groupselect_norm]
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
        observeEvent(input$start_preprocessing_microarray_norm,{

          # Show modal
          show_modal_spinner(text = "Pre-processing data...",
                             color="#0dc5c1")
          
          # Select outlier
          if (!isTRUE(input$outier)){
            rv$outlier <- input$select_outliers_norm
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
          rv$metaData_fil <- rv$metaData[str_remove(colnames(exprs(rv$gxData_fil)),"\\.CEL.*"),]
          
          # Experiment factor
          if(length(input$groupselect) > 1){
            rv$experimentFactor <- factor(apply(rv$metaData_fil[,input$groupselect_norm], 1, paste, collapse = "_" ))
          } else{
            rv$experimentFactor <- factor(rv$metaData_fil[,input$groupselect_norm])
          }
          
          # Normalization
          rv$normData <- microarrayNormalization_processed(gxData = rv$gxData_fil,
                                                           experimentFactor = rv$experimentFactor_norm,
                                                           transMeth = input$transformation_microarray_norm,
                                                           normMeth = input$normMeth_microarray_norm,
                                                           perGroup_name = input$perGroup_microarray_norm)
          
          rv$normMatrix <- exprs(rv$normData)

          #======================================================================#
          # QC
          #======================================================================#
          
          #********************************************************************#
          # Expression values
          #********************************************************************#
          
          # Print expression table
          output$exprTable_microarray_norm_norm <- renderDataTable({
            
            # Remove modal
            remove_modal_spinner()
            
            # Show message
            sendSweetAlert(
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
          
          # Boxplot of single gene
          output$ExprBoxplot_microarray_norm <- renderPlotly({
            req(input$exprTable_microarray_norm_norm_rows_selected)
            
            # Set colors
            myPalette <- colorsByFactor(rv$experimentFactor)
            legendColors <- myPalette$legendColors
            
            # Prepare expression data frame
            plotExpr <- data.frame(
              logExpr = rv$normMatrix[input$exprTable_microarray_norm_norm_rows_selected,],
              Grouping = rv$experimentFactor
            )
            
            # make interactive plot
            plotExpr %>%
              plot_ly(x = ~Grouping,y = ~as.numeric(logExpr),
                      color = ~Grouping, colors = legendColors, type = "box") %>%
              layout(xaxis = list(title = " "),
                     yaxis = list(title = 'log intensity'),
                     legend = list(title=list(text='Group')),
                     showlegend = FALSE)
            
          })
          
          #********************************************************************#
          # Boxplots
          #********************************************************************#
          
          # Boxplots of all genes together
          output$boxplots_microarray_norm_norm <- renderImage({
            getBoxplots(experimentFactor = rv$experimentFactor,
                        normData = rv$normData)
          },deleteFile = TRUE)
          
          #********************************************************************#
          # Densityplots
          #********************************************************************#
          
          # Densityplots of all genes together
          output$densityplots_microarray_norm_norm <- renderPlotly({
            getDensityplots(experimentFactor = rv$experimentFactor,
                            normMatrix = rv$normMatrix)
            
          })
          
          #********************************************************************#
          # Heatmap
          #********************************************************************#
          
          output$heatmap_microarray_norm_norm  <- renderPlotly({
            myPalette <- colorsByFactor(rv$experimentFactor)
            
            # Plot colors
            plotColors <- myPalette$plotColors
            
            # Legend colors
            legendColors <- myPalette$legendColors
            
            # Set cluster options
            clusterOption1 <- input$clusteroption1_microarray_norm
            clusterOption2 <- input$clusteroption2_microarray_norm
            
            #note: for computing array correlation, euclidean would not make sense
            #only use euclidean distance to compute the similarity of the correlation
            #vectors for the arrays
            COpt1 <- "pearson"
            if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
            crp <- cor(rv$normMatrix, use="complete.obs", method=COpt1)
            
            switch(tolower(clusterOption1),
                   "pearson" = {
                     my.dist <- function(x) cor.dist(x, abs=FALSE)
                   },
                   "spearman" = {
                     my.dist <- function(x) spearman.dist(x, abs=FALSE)
                   },
                   "euclidean" = {
                     my.dist <- function(x) euc(x)
                   }
            )
            
            # Perform clustering
            my.hclust <- function(d) hclust(d, method=clusterOption2)
            names(legendColors) <- levels(rv$experimentFactor)
            sidecolors <- data.frame(rv$experimentFactor)
            colnames(sidecolors) <- "Experimental group"
            
            if (input$heatmaptheme_microarray_norm == "Default"){
              gradient = viridis(n = 256, alpha = 1, begin = 0, end = 1,
                                 option = "viridis")
            }
            
            if (input$heatmaptheme_microarray_norm == "Yellow-red"){
              gradient = heat.colors(100)
            }
            
            if (input$heatmaptheme_microarray_norm == "Dark"){
              gradient = RColorBrewer::brewer.pal(8, "Dark2")
            }
            
            # Make heatmap
            p <- heatmaply(crp, plot_method = "plotly", distfun = my.dist,
                           hclustfun = my.hclust, symm = TRUE, seriate = "mean",
                           titleX = FALSE, titleY = FALSE, key.title = NULL,
                           show_dendrogram = c(TRUE, FALSE), col_side_colors = sidecolors,
                           col_side_palette = legendColors, column_text_angle = 90,
                           colors = gradient)
            
            return(p)
            
          })
          
          #********************************************************************#
          # PCA
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
            print(head(colorFactor))
            
            plot_PCA(PC_data = rv$PCA_data, 
                     colorFactor = colorFactor, 
                     xpc = as.numeric(str_remove(input$xpca_microarray_norm,"PC")), 
                     ypc = as.numeric(str_remove(input$ypca_microarray_norm,"PC")), 
                     zpc = ifelse(input$xyz_microarray_norm,as.numeric(str_remove(input$zpca_microarray_norm,"PC")),3), 
                     xyz = input$xyz_microarray_norm)
            
          })
          #********************************************************************#
          # UI
          #********************************************************************#
          
          # Render UI for main tab
          output$UI_QC_microarray_norm_norm <- renderUI({
            tagList(
              tabsetPanel(
                tabPanel("Expression values",
                         icon = icon("fas fa-mouse-pointer"),
                         h3(strong("Normalized expression values")),
                         h5("Here you can view the normalized log intensity (expression) values."),
                         hr(),
                         dataTableOutput(outputId = "exprTable_microarray_norm_norm") %>% 
                           withSpinner(color="#0dc5c1"),
                         downloadButton("downloadNormalizedData_microarray_norm", 
                                        "Download"),
                         plotlyOutput("ExprBoxplot_microarray_norm")%>% 
                           withSpinner(color="#0dc5c1")
                ),
                tabPanel("Boxplots",
                         icon = icon("fas fa-file"),
                         h3(strong("Boxplots")),
                         h5("These are boxplots of the normalized expression expression values."),
                         hr(),
                         plotOutput(outputId = "boxplots_microarray_norm_norm") %>% 
                           withSpinner(color="#0dc5c1")
                ),
                
                tabPanel("Density plots",
                         icon = icon("fas fa-mouse-pointer"),
                         h3(strong("Density plots")),
                         h5("These are density plots of the normalized expression expression values."),
                         hr(),
                         plotlyOutput(outputId = "densityplots_microarray_norm_norm") %>% 
                           withSpinner(color="#0dc5c1")
                ),
                
                tabPanel("Correlation plot", 
                         icon = icon("fas fa-mouse-pointer"),
                         
                         br(),
                         
                         fluidRow(
                           
                           #Distance method
                           column(4,
                                  pickerInput(
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
                                  pickerInput(
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
                           
                           selectInput(inputId = 'heatmaptheme_microarray_norm',
                                       label = NULL,
                                       choices = c("Default", 
                                                   "Yellow-red", 
                                                   "Dark")),
                           
                           circle = TRUE, status = "info",
                           icon = icon("fas fa-cog"), width = "300px",
                           
                           tooltip = tooltipOptions(
                             title = "Click to change theme!")
                           
                         ),
                         
                         plotlyOutput("heatmap_microarray_norm_norm", 
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
                                  pickerInput(inputId = "colorFactor_microarray_norm",
                                              label = "Color by:",
                                              choices = colnames(rv$metaData_fil),
                                              selected = autoGroup(rv$metaData_fil),
                                              multiple = TRUE)
                           ),
                           column(3,
                                  br(),
                                  # 3D plot?
                                  materialSwitch(
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
                                    condition = "input.xyz_microarray_raw==true",
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
                           plotlyOutput("PCA_microarray_norm")%>% 
                             withSpinner(color="#0dc5c1")
                         ) 
                         
                ) # PCA tabpanel
              ) # tabsetpanel
            ) # taglist
          }) # renderUI
          
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
        
        # Select experimental group
        output$UI_expFactor_microarray_norm <- renderUI({
          tagList(
            pickerInput(inputId = "expFactor_microarray_norm",
                        label = "Experiment factor:",
                        choices = colnames(rv$metaData_fil),
                        selected = autoGroup(rv$metaData_fil),
                        multiple = TRUE)
          )
        })
        
        # print the experimental levels
        output$expFactor_levels_micorarray_norm <- renderText({
          req(input$expFactor_microarray_norm)
          
          if(length(input$expFactor_microarray_norm) > 1){
            experimentFactor <- apply(rv$metaData_fil[,input$expFactor_microarray_norm], 1, paste, collapse = "_" )
          } else{
            experimentFactor <- rv$metaData_fil[,input$expFactor_microarray_norm]
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
        
        # Select continuous covariates
        output$UI_covGroups_num_microarray_norm <- renderUI({
          req(input$expFactor_microarray_norm)
          tagList(
            pickerInput(inputId = "covGroups_num_microarray_norm",
                        label = "Continuous covariates (e.g., age):",
                        choices = setdiff(colnames(rv$metaData_fil),
                                          input$expFactor_microarray_norm),
                        selected = NULL,
                        multiple = TRUE)
          )
        })
        
        # Select discrete covariates
        output$UI_covGroups_char_microarray_norm <- renderUI({
          req(input$expFactor_microarray_norm)
          tagList(
            pickerInput(inputId = "covGroups_char_microarray_norm",
                        label = "Discrete covariates (e.g., sex):",
                        choices = setdiff(colnames(rv$metaData_fil),
                                          input$expFactor_microarray_norm),
                        selected = NULL,
                        multiple = TRUE)
          )
        })
        
        # Select comparisons
        output$UI_comparisons_microarray_norm <- renderUI({
          req(input$expFactor_microarray_norm)
          tagList(
            multiInput(
              inputId = "comparisons_microarray_norm",
              label = "Comparisons:", 
              choices = makeComparisons(make.names(unique(rv$metaData[,input$expFactor_microarray_norm]))),
              selected = makeComparisons(make.names(unique(rv$metaData[,input$expFactor_microarray_norm])))[1]
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
                        choices = c("ensembl_gene_id",
                                    "entrezgene_id",
                                    "gene_name"),
                        selected = "gene_name",
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
            rv$top_table <- getStatistics(normMatrix = rv$normMatrix, 
                                          metaData = rv$metaData_fil, 
                                          expFactor = input$expFactor_microarray_norm, 
                                          covGroups_num = input$covGroups_num_microarray_norm,
                                          covGroups_char = input$covGroups_char_microarray_norm,
                                          comparisons = input$comparisons_microarray_norm,
                                          addAnnotation = input$addAnnotation_microarray_norm,
                                          biomart_dataset = input$biomart_dataset_microarray_norm,
                                          biomart_attributes = unique(c(input$biomart_filter_microarray_norm,
                                                                        input$biomart_attributes_microarray_norm)),
                                          biomart_filters = input$biomart_filter_microarray_norm)
            
          } else{
            rv$top_table <- getStatistics(normMatrix = rv$normMatrix, 
                                          metaData = rv$metaData_fil, 
                                          expFactor = input$expFactor_microarray_norm, 
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
            #req(rv$top_table)
            
            if (!is.null(rv$top_table)){
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
                text = "Statistical analysis has been performed. You can download 
              the results as well as view them in interactive plots.",
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
          
          # Boxplot of single gene
          output$ExprBoxplot_statistics_microarray_norm <- renderPlotly({
            
            req(input$top_table_microarray_norm_rows_selected)
            req(input$comparisons_view_microarray_norm)
            
            # Set colors
            myPalette <- colorsByFactor(rv$experimentFactor)
            legendColors <- myPalette$legendColors
            
            # Prepare expression data frame
            gene <- rv$top_table[[input$comparisons_view_microarray_norm]]$GeneID[input$top_table_microarray_norm_rows_selected]
            plotExpr <- data.frame(
              logExpr = rv$normMatrix[as.character(rownames(rv$normMatrix)) %in% as.character(gene),],
              Grouping = factor(rv$metaData_fil[,input$expFactor_microarray_norm])
            )
            
            # make interactive plot
            plotExpr %>%
              plot_ly(x = ~Grouping,y = ~as.numeric(logExpr),
                      color = ~Grouping, colors = legendColors, type = "box") %>%
              layout(xaxis = list(title = " "),
                     yaxis = list(title = 'log intensity'),
                     legend = list(title=list(text='Group')),
                     showlegend = FALSE)
            
          })
          
          #********************************************************************#
          # Histograms
          #********************************************************************#
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
          
          #********************************************************************#
          # Volcano plot
          #********************************************************************#
          observeEvent(input$plot_volcano_microarray_norm, {
            req(rv$top_table)
            req(input$rawp_volcano_microarray_norm)
            req(input$p_thres_volcano_microarray_norm)
            req(input$logFC_thres_volcano_microarray_norm)
            p <- makeVolcano(top_table = rv$top_table[[input$comparisons_view_microarray_norm]], 
                             p = input$rawp_volcano_microarray_norm, 
                             p_threshold = input$p_thres_volcano_microarray_norm, 
                             logFC_threshold = input$logFC_thres_volcano_microarray_norm)
            
            output$volcano_microarray_norm <- renderPlotly(p)
          }, ignoreNULL = FALSE)
          
          
          # UI: Output in different tabs
          
          observe({
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
                                          "Download"),
                           plotlyOutput("ExprBoxplot_statistics_microarray_norm")%>% 
                             withSpinner(color="#0dc5c1")
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
                          choices = c("ensembl_gene_id" = "ENSEMBL", 
                                      "entrezgene_id" = "ENTREZID", 
                                      "gene_name" = "SYMBOL"),
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
          rv$ORA_data <- ORA(top_table = rv$top_table[[input$comparisons_view_ORA_microarray_norm]],
                             geneset = input$geneset_ORA_microarray_norm,
                             geneID_col = input$geneID_ORA_microarray_norm,
                             geneID_type = input$selID_ORA_microarray_norm,
                             organism = input$organism_ORA_microarray_norm,
                             rawadj = input$rawp_ORA_microarray_norm,
                             updown = input$updown_ORA_microarray_norm,
                             p_thres = input$p_thres_ORA_microarray_norm,
                             logFC_thres = input$logFC_thres_ORA_microarray_norm)
          
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
                
                if (input$geneset_ORA_microarray_norm != "WikiPathways"){
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