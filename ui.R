#==============================================================================#
# Name: ui.R
# Description: User interface of the ArrayAnalysis Shiny app
#==============================================================================#

# Start UI
ui <- tagList(
  
  # Set the style of the UI
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                            .navbar-nav > li:nth-child(18) {
                           float: right;
                           }
                           .my_style_1{ 
                             background-image: url(background.jpg);
                           }
                           
                           .my_style_1 { margin-top: -20px; }
                           
                           .my_style_1 { width: 100%; }
                           
                           .container-fluid { padding-left: 0; padding-right: 0; }
                           
                           .my_style_1 { position: absolute; left: 0; }
                           
                           "))),
  
  
  fluidPage(
    
    # This allow for pop-up messages to show up
    shinyWidgets::useSweetAlert(),
    
    # This allows for insertion of information boxes in the app
    prompter::use_prompt(),
    
    # Make a page with a navigation bar
    navbarPage(title = "ArrayAnalysis", id = "navbar",
               
               ###################################################################
               
               #  Home panel                        
               
               ###################################################################
               
               tabPanel("Home", 
                        value = "home_panel", 
                        icon = icon("fas fa-home"), 
                        class = "my_style_1",
                        
                        # Set spacing between top of page and text box
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        
                        # Make text box
                        fluidRow(
                          
                          column(6, offset = 3, 
                                 align = "center", 
                                 style = "background-color:#FFFFFF;",
                                 
                                 # Line break
                                 br(),
                                 
                                 # ArrayAnalysis logo
                                 img(src = "logo.png", width = "100%"),
                                 
                                 # Welcome message
                                 h1(strong(span(style = "color:#000000", 
                                                "Welcome to ArrayAnalysis!"))),
                                 h5(span(style = "color:#000000", 
                                         "Do you want to analyze microarray or RNA-Seq data?")),
                                 
                                 # Line break
                                 br(),
                                 
                                 # Analyse microarray or RNA-seq data?
                                 radioGroupButtons(
                                   inputId = "microarray_or_rnaseq",
                                   label = NULL,
                                   choices = c("Microarray", 
                                               "RNA-Seq",
                                               "Documentation"),
                                   status = "danger",
                                   selected = "Microarray"
                                 ),
                                 
                                 # Analyse raw or processed data?
                                 prettyRadioButtons(
                                   inputId = "raw_or_norm",
                                   label = NULL, 
                                   choices = c("Raw data", "Processed data"),
                                   inline = TRUE, 
                                   status = "danger",
                                   fill = TRUE),
                                 
                                 # Line break
                                 br(),
                                 
                                 # Action button: start analysis by clicking
                                 actionBttn(inputId = "startAnalysis",
                                            label = "Start Analysis",
                                            style = "simple",
                                            color = "primary",
                                            icon = icon("arrow-right")),
                                 
                                 # Line breaks
                                 br(),
                                 br(),
                                 br()
                                 
                          ) # EO column
                        ), # EO fluidRow
                        
                        # Set spacing between bottom of page and text box
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br()
               ), # EO home panel
               
               ###################################################################
               
               #  Microarray: raw
               
               ###################################################################
               
               #*****************************************************************#
               # Upload microarray data
               #*****************************************************************#
               
               tabPanel("Upload", value = "panel_upload_microarray_raw", 
                        icon = icon("fas fa-upload"),
                        
                        # Side panel
                        sidebarPanel(
                          
                          #Title
                          h2(strong("Data upload")),
                          h5("Before you can run the analysis workflow, you first 
                           need to upload the expression data as well as the meta data."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as an .zip folder
                           containing all '.CEL.gz' files. The file names should match with
                           the sample IDs in the meta data."),
                          fileInput(inputId = "uploadCEL_microarray_raw",
                                    label = NULL,
                                    accept = ".zip",
                                    placeholder = "Select .zip data file"),
                          
                          h4(strong("2. Upload meta data")),
                          h5("The meta data includes relevant information (e.g., diagnostic group)
                           about the samples. You can upload the meta data as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example', 'here'),
                             "for an example .csv meta data file. A Series Matrix File
                           can be downloaded from the GEO website."),
                          
                          prettyRadioButtons(inputId = "MetaFileType", 
                                             label = NULL, 
                                             choices = c(".tsv/.csv file",
                                                         "Series Matrix File"),
                                             inline = TRUE,
                                             fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_microarray_raw_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file"),
                            
                          ),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_microarray_raw_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File"),
                            
                          ),
                          
                          # Confirm upload
                          actionBttn(inputId = "upload_microarray_raw",
                                     label = "Read data",
                                     style = "simple",
                                     color = "primary",
                                     icon = icon("fas fa-upload")),
                          actionBttn(inputId = "example_microarray_raw",
                                     label = "Run example",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("cloud-arrow-up")),
                          
                          
                          uiOutput("next_upload_microarray_raw")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          uiOutput("UI_upload_microarray_raw")
                        )
                        
               ), # End of upload microarray raw tab
               
               #*****************************************************************#
               # Pre-processing of microarray data
               #*****************************************************************#
               tabPanel("Pre-processing", value = "panel_preprocessing_microarray_raw", icon = icon("sync"),
                        
                        sidebarPanel(
                          
                          #******************************************************#
                          #   Information
                          #******************************************************#
                          
                          h2(strong("Pre-processing")),
                          
                          h5("In this pre-processing step, you can remove samples 
                           (e.g., outliers), perform normalization, 
                           and choose your desired probeset annotation."),
                          
                          hr(),
                          
                          #******************************************************#
                          #   Options
                          #******************************************************#
                          
                          #Remove outliers
                          h4(strong(tags$span(
                            "1. Remove samples",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Uncheck the box to exclude samples from the analysis.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          awesomeCheckbox(inputId = "outlier",
                                          label = "Keep all samples", 
                                          value = TRUE,
                                          status = "danger"),
                          
                          uiOutput("UI_outlier_microarray_raw"),
                          br(),
                          
                          # Select experimental group
                          h4(strong(tags$span(
                            "2. Select experimental group",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_microarray_raw"),
                          
                          htmlOutput("experimentallevels"),
                          br(),
                          
                          #Normalization
                          h4(strong(tags$span(
                            "3. Normalization",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all arrays or per experimental group.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          radioGroupButtons(inputId = "normMeth_microarray", 
                                            label = NULL, 
                                            choices = c("RMA","GCRMA","PLIER"), #"None"
                                            status = "danger"),
                          
                          prettyRadioButtons(inputId = "perGroup_microarray", 
                                             label = NULL, 
                                             choices = c("Use all arrays",
                                                         "Per experimental group"),
                                             inline = TRUE, 
                                             status = "danger",
                                             fill = TRUE),
                          
                          br(),
                          
                          #Annotation
                          h4(strong(tags$span(
                            "4. Annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the probeset annotation. 
                                       ENTREZG or ENSG custom annotations are recommended.
                                       If no annotation is selected, the affy annotation will be used.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          selectInput(inputId = "annotations_microarray",
                                      label = NULL,
                                      choices = c("No annotations",
                                                  "Custom annotations",
                                                  "Upload annotation file"),
                                      selected = "Custom annotations"),
                          
                          conditionalPanel(
                            condition = "input.annotations_microarray=='Custom annotations'",
                            
                            uiOutput("UI_species_microarray_raw"),
                            
                            selectInput(inputId = "CDFtype_microarray",
                                        label = "Annotation format",
                                        choices = c("ENTREZG","REFSEQ","ENSG",
                                                    "ENSE","ENST","VEGAG","VEGAE",
                                                    "VEGAT","TAIRG","TAIRT","UG",
                                                    "MIRBASEF","MIRBASEG"))
                          ),
                          conditionalPanel(
                            condition = "input.annotations_microarray=='Upload annotation file'",
                            
                            fileInput(inputId = "annot_file_microarray",
                                      label = "Upload annotation file",
                                      multiple = FALSE)
                          ),
                          
                          br(),
                          
                          #Pre-processing
                          h4(strong(tags$span(
                            "5. Pre-processing",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to perform the pre-processing!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          actionBttn(inputId = "start_preprocessing_microarray_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_preprocessing_microarray_raw")
                          
                        ),
                        
                        #********************************************************#
                        #   Output panel
                        #********************************************************#
                        mainPanel(
                          uiOutput("UI_QC_microarray_raw")
                          
                          
                        ) #Main panel
                        
               ), #Tab panel
               
               #*****************************************************************#
               # Statistical analysis microarray data
               #*****************************************************************#
               tabPanel("Statistical analysis", value = "panel_statistics_microarray_raw", 
                        icon = icon("fas fa-chart-bar"),
                        
                        sidebarPanel(
                          h2(strong("Statistical analysis")),
                          
                          h5("In the statistical analysis step, 
                        you can select which groups to compare to each other and 
                           which covariates to add to the statistical model."),
                          
                          hr(),
                          h4(strong(tags$span(
                            "1. Make comparisons",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_expFactor_microarray_raw"),
                          uiOutput("UI_comparisons_microarray_raw"),
                          htmlOutput("expFactor_levels_micorarray_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariated",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Adjust for variables like age and sex 
                                       by adding them as covariates to the model.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_covGroups_num_microarray_raw"),
                          uiOutput("UI_covGroups_char_microarray_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "3. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Add gene annotations from biomaRt to 
                                       the output.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          awesomeCheckbox(inputId = "addAnnotation_microarray_raw",
                                          label = "Add gene annotations from biomaRt", 
                                          value = FALSE,
                                          status = "danger"),
                          conditionalPanel(
                            condition = "input.addAnnotation_microarray_raw==true",
                            uiOutput("UI_biomart_dataset_microarray_raw"),
                            uiOutput("UI_addAnnotations_microarray_raw")
                          ),
                          br(),
                          
                          h4(strong(tags$span(
                            "4. Perform statistical analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the statistical analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_statistics_microarray_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_statistics_microarray_raw")
                          
                        ),
                        mainPanel(
                          uiOutput("UI_comparisons_view_microarray_raw"),
                          br(),
                          uiOutput("UI_output_statistics_microarray_raw")
                        )
               ), # Tab panel
               
               #*****************************************************************#
               # ORA microarray data
               #*****************************************************************#
               tabPanel("ORA", value = "panel_ORA_microarray_raw", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Overrepresentation analysis")),
                          
                          h5("In the gene overrepresentation analysis (ORA) step, 
                        you can find processes that are enriched by the 
                           differentially expressed genes (DEGs)."),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Choose from the comparisons for which 
                                       the statistical analysis was performed in the previous step.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_microarray_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select geneset collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "A geneset is a collection of genes that are association 
                                       with a specific biological process (GO-BP), molecular function (GO-MF),
                                       cellular component (GO-CC), or biological pathway (WikiPathways).", 
                                         position = "right",
                                         size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_microarray_raw",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways")),
                          
                          br(),
                          h4(strong(tags$span(
                            "3. Select differentially expressed genes",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Tip: look at the volcano plot in the previous 
                                       step (statistical analysis) to find the optimal P value and 
                                       logFC thresholds.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          # Up/down regulated genes only
                          prettyRadioButtons(
                            inputId = "updown_ORA_microarray_raw",
                            label = "Perform ORA on ...", 
                            choices = 
                              c("Upregulated genes only", 
                                "Downregulated genes only",
                                "Both"),
                            selected = "Both",
                            status = "danger",
                            fill = TRUE),
                          
                          # top N or logFC/P value threshold
                          radioGroupButtons(inputId = "topNorThres_microarray_raw",
                                            label = "Select DEGs based on ...", 
                                            choices = c("Threshold", 
                                                        "Top N"),
                                            status = "danger"),
                          
                          # Select DEGs based on logFC/P value threshold
                          conditionalPanel(
                            condition = "input.topNorThres_microarray_raw=='Threshold'",
                            
                            #P value threshold
                            numericInput(
                              inputId = "p_thres_ORA_microarray_raw",
                              label = "P threshold",
                              value = 0.05),
                            
                            # Raw or adjusted P value?
                            prettyRadioButtons(inputId = "rawp_ORA_microarray_raw", 
                                               label = NULL, 
                                               choices = c("Raw P value" = "raw", 
                                                           "Adjusted P value" = "adj"),
                                               inline = TRUE, 
                                               status = "danger",
                                               fill = TRUE),
                            
                            #logFC threshold
                            numericInput(
                              inputId = "logFC_thres_ORA_microarray_raw",
                              label = "logFC threshold",
                              value = 0),
                          ),
                          
                          # Select top N most significant genes
                          conditionalPanel(
                            condition = "input.topNorThres_microarray_raw=='Top N'",
                            numericInput(
                              inputId = "topN_microarray_raw",
                              label = "Top N most significant genes",
                              value = 100),
                          ),

                          br(),
                          
                          # Select gene identifier
                          h4(strong(tags$span(
                            "4. Select gene identifier",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "We need to know which gene identifiers are 
                                       used, so we can link the genes to their correct genesets.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_geneID_ORA_microarray_raw"),
                          br(),
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform ORA",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the overrepresentation analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_ORA_microarray_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync"))
                          
                        ),
                        mainPanel(
                          uiOutput("UI_output_ORA_microarray_raw")
                        )
               ), # Tab panel
               
               ###################################################################
               
               #  Microarray: processed
               
               ###################################################################
               
               #*****************************************************************#
               # Upload microarray data
               #*****************************************************************#
               
               tabPanel("Upload", value = "panel_upload_microarray_norm", 
                        icon = icon("fas fa-upload"),
                        
                        # Side panel
                        sidebarPanel(
                          
                          #Title
                          h2(strong("Data upload")),
                          h5("Before you can run the analysis workflow, you first 
                           need to upload the expression data as well as the meta data."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a Series 
                           Matrix File."),
                          fileInput(inputId = "uploadExprData_microarray_norm_smf",
                                    label = NULL,
                                    accept = ".txt.gz",
                                    placeholder = "Select .txt.gz Series Matrix File"),
                          
                          h4(strong("2. Upload meta data")),
                          h5("The meta data includes relevant information (e.g., diagnostic group)
                           about the samples. You can upload the meta data as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_norm', 'here'),
                             "for an example .csv meta data file. A Series Matrix File
                           can be downloaded from the GEO website."),
                          
                          prettyRadioButtons(inputId = "MetaFileType_norm", 
                                             label = NULL, 
                                             choices = c(".tsv/.csv file",
                                                         "Series Matrix File"),
                                             inline = TRUE,
                                             fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_norm=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_microarray_norm_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file")
                            
                          ),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_norm=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_microarray_norm_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File")
                            
                          ),
                          
                          # Confirm upload
                          actionBttn(inputId = "upload_microarray_norm",
                                     label = "Read data",
                                     style = "simple",
                                     color = "primary",
                                     icon = icon("fas fa-upload")),
                          actionBttn(inputId = "example_microarray_norm",
                                     label = "Run example",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("cloud-arrow-up")),
                          
                          
                          uiOutput("next_upload_microarray_norm")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          uiOutput("UI_upload_microarray_norm")
                        )
                        
               ), # End of upload microarray norm tab
               
               #*****************************************************************#
               # Pre-processing of microarray data
               #*****************************************************************#
               tabPanel("Pre-processing", value = "panel_preprocessing_microarray_norm", icon = icon("sync"),
                        
                        sidebarPanel(
                          
                          #******************************************************#
                          #   Information
                          #******************************************************#
                          
                          h2(strong("Pre-processing")),
                          
                          h5("In this pre-processing step, you can remove samples 
                           (e.g., outliers), perform normalization, 
                           and choose your desired probeset annotation."),
                          
                          hr(),
                          
                          #******************************************************#
                          #   Options
                          #******************************************************#
                          
                          #Remove outliers
                          h4(strong(tags$span(
                            "1. Remove samples",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Uncheck the box to exclude samples from the analysis.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          awesomeCheckbox(inputId = "outlier_norm",
                                          label = "Keep all samples", 
                                          value = TRUE,
                                          status = "danger"),
                          
                          uiOutput("UI_outlier_microarray_norm"),
                          br(),
                          
                          # Select experimental group
                          h4(strong(tags$span(
                            "2. Select experimental group",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_microarray_norm"),
                          
                          htmlOutput("experimentallevels_norm"),
                          br(),
                          
                          #Transformation
                          h4(strong(tags$span(
                            "3. Transformation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the transformation method. 
                                       The aim of data transformation is to make the 
                                       data more normally distributed.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_transformation_microarray_norm"),
                          br(),
                          #Normalization
                          h4(strong(tags$span(
                            "4. Normalization",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all arrays or per experimental group.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          radioGroupButtons(inputId = "normMeth_microarray_norm", 
                                            label = NULL, 
                                            choices = c("Quantile","None"),
                                            selected = "None",
                                            status = "danger"),
                          
                          prettyRadioButtons(inputId = "perGroup_microarray_norm", 
                                             label = NULL, 
                                             choices = c("Use all arrays",
                                                         "Per experimental group"),
                                             inline = TRUE, 
                                             status = "danger",
                                             fill = TRUE),
                          
                          br(),
                          
                          #Pre-processing
                          h4(strong(tags$span(
                            "5. Pre-processing",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to perform the pre-processing!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          actionBttn(inputId = "start_preprocessing_microarray_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_preprocessing_microarray_norm")
                          
                        ),
                        
                        #********************************************************#
                        #   Output panel
                        #********************************************************#
                        mainPanel(
                          uiOutput("UI_QC_microarray_norm")
                          
                          
                        ) #Main panel
                        
               ), #Tab panel
               
               #*****************************************************************#
               # Statistical analysis microarray data
               #*****************************************************************#
               tabPanel("Statistical analysis", value = "panel_statistics_microarray_norm", 
                        icon = icon("fas fa-chart-bar"),
                        
                        sidebarPanel(
                          h2(strong("Statistical analysis")),
                          
                          h5("In the statistical analysis step, 
                        you can select which groups to compare to each other and 
                           which covariates to add to the statistical model."),
                          
                          hr(),
                          h4(strong(tags$span(
                            "1. Make comparisons",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_expFactor_microarray_norm"),
                          uiOutput("UI_comparisons_microarray_norm"),
                          htmlOutput("expFactor_levels_micorarray_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariated",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Adjust for variables like age and sex 
                                       by adding them as covariates to the model.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_covGroups_num_microarray_norm"),
                          uiOutput("UI_covGroups_char_microarray_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "3. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Add gene annotations from biomaRt to 
                                       the output.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          awesomeCheckbox(inputId = "addAnnotation_microarray_norm",
                                          label = "Add gene annotations from biomaRt", 
                                          value = FALSE,
                                          status = "danger"),
                          conditionalPanel(
                            condition = "input.addAnnotation_microarray_norm==true",
                            uiOutput("UI_biomart_dataset_microarray_norm"),
                            uiOutput("UI_addAnnotations_microarray_norm")
                          ),
                          br(),
                          
                          h4(strong(tags$span(
                            "4. Perform statistical analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the statistical analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_statistics_microarray_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_statistics_microarray_norm")
                          
                        ),
                        mainPanel(
                          uiOutput("UI_comparisons_view_microarray_norm"),
                          br(),
                          uiOutput("UI_output_statistics_microarray_norm")
                        )
               ), # Tab panel
               
               #*****************************************************************#
               # ORA microarray data
               #*****************************************************************#
               tabPanel("ORA", value = "panel_ORA_microarray_norm", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Overrepresentation analysis")),
                          
                          h5("In the gene overrepresentation analysis (ORA) step, 
                        you can find processes that are enriched by the 
                           differentially expressed genes."),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Choose from the comparisons for which 
                                       the statistical analysis was performed in the previous step.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_microarray_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select geneset collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "A geneset is a collection of genes that are association 
                                       with a specific biological process (GO-BP), molecular function (GO-MF),
                                       cellular component (GO-CC), or biological pathway (WikiPathways).", 
                                         position = "right",
                                         size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_microarray_norm",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways")),
                          
                          br(),
                          
                          h4(strong(tags$span(
                            "3. Select differentially expressed genes",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Tip: look at the volcano plot in the previous 
                                       step (statistical analysis) to find the optimal P value and 
                                       logFC thresholds.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          # Up/down regulated genes only
                          prettyRadioButtons(
                            inputId = "updown_ORA_microarray_norm",
                            label = "Perform ORA on ...", 
                            choices = 
                              c("Upregulated genes only", 
                                "Downregulated genes only",
                                "Both"),
                            selected = "Both",
                            status = "danger",
                            fill = TRUE),
                          
                          # top N or logFC/P value threshold
                          radioGroupButtons(inputId = "topNorThres_microarray_norm",
                                            label = "Select DEGs based on ...", 
                                            choices = c("Threshold", 
                                                        "Top N"),
                                            status = "danger"),
                          
                          # Select DEGs based on logFC/P value threshold
                          conditionalPanel(
                            condition = "input.topNorThres_microarray_norm=='Threshold'",
                            
                            #P value threshold
                            numericInput(
                              inputId = "p_thres_ORA_microarray_norm",
                              label = "P threshold",
                              value = 0.05),
                            
                            # Raw or adjusted P value?
                            prettyRadioButtons(inputId = "rawp_ORA_microarray_norm", 
                                               label = NULL, 
                                               choices = c("Raw P value" = "raw", 
                                                           "Adjusted P value" = "adj"),
                                               inline = TRUE, 
                                               status = "danger",
                                               fill = TRUE),
                            
                            #logFC threshold
                            numericInput(
                              inputId = "logFC_thres_ORA_microarray_morm",
                              label = "logFC threshold",
                              value = 0),
                          ),
                          
                          # Select top N most significant genes
                          conditionalPanel(
                            condition = "input.topNorThres_microarray_norm=='Top N'",
                            numericInput(
                              inputId = "topN_microarray_norm",
                              label = "Top N most significant genes",
                              value = 100),
                          ),
                          
                          br(),
                          
                          # Select gene identifier
                          h4(strong(tags$span(
                            "4. Select gene identifier",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "We need to know which gene identifiers are 
                                       used, so we can link the genes to their correct genesets.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_geneID_ORA_microarray_norm"),
                          br(),
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform ORA",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the overrepresentation analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_ORA_microarray_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync"))
                          
                        ),
                        mainPanel(
                          uiOutput("UI_output_ORA_microarray_norm")
                        )
               ), # Tab panel
               
               ###################################################################
               
               #  RNA-seq: raw
               
               ###################################################################
               
               #*****************************************************************#
               # Upload RNA-seq data
               #*****************************************************************#
               
               tabPanel("Upload", value = "panel_upload_rnaseq_raw", 
                        icon = icon("fas fa-upload"),
                        
                        # Side panel
                        sidebarPanel(
                          
                          #Title
                          h2(strong("Data upload")),
                          h5("Before you can run the analysis workflow, you first 
                           need to upload the expression data as well as the meta data."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a .tsv/.csv file."),
                          fileInput(inputId = "uploadExprData_rnaseq_raw",
                                    label = NULL,
                                    accept = c(".tsv",".csv"),
                                    placeholder = "Select .tsv/.csv file"),
                          
                          h4(strong("2. Upload meta data")),
                          h5("The meta data includes relevant information (e.g., diagnostic group)
                           about the samples. You can upload the meta data as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_rnaseq_raw', 'here'),
                             "for an example .csv meta data file. A Series Matrix File
                           can be downloaded from the GEO website."),
                          
                          prettyRadioButtons(inputId = "MetaFileType_rnaseq_raw", 
                                             label = NULL, 
                                             choices = c(".tsv/.csv file",
                                                         "Series Matrix File"),
                                             inline = TRUE,
                                             fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_rnaseq_raw=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_rnaseq_raw_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file")
                            
                          ),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_rnaseq_raw=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_rnaseq_raw_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File")
                            
                          ),
                          
                          # Confirm upload
                          actionBttn(inputId = "upload_rnaseq_raw",
                                     label = "Read data",
                                     style = "simple",
                                     color = "primary",
                                     icon = icon("fas fa-upload")),
                          actionBttn(inputId = "example_rnaseq_raw",
                                     label = "Run example",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("cloud-arrow-up")),
                          
                          
                          uiOutput("next_upload_rnaseq_raw")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          uiOutput("UI_upload_rnaseq_raw")
                        )
                        
               ), # End of upload RNA-seq (raw) tab
               
               #*****************************************************************#
               # Pre-processing of RNA-seq data
               #*****************************************************************#
               tabPanel("Pre-processing", value = "panel_preprocessing_rnaseq_raw", icon = icon("sync"),
                        
                        sidebarPanel(
                          
                          #******************************************************#
                          #   Information
                          #******************************************************#
                          
                          h2(strong("Pre-processing")),
                          
                          h5("In this pre-processing step, you can remove samples 
                           (e.g., outliers), perform normalization, 
                           and choose your desired probeset annotation."),
                          
                          hr(),
                          
                          #******************************************************#
                          #   Options
                          #******************************************************#
                          
                          #Remove outliers
                          h4(strong(tags$span(
                            "1. Remove samples",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Uncheck the box to exclude samples from the analysis.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          awesomeCheckbox(inputId = "outlier_rnaseq_raw",
                                          label = "Keep all samples", 
                                          value = TRUE,
                                          status = "danger"),
                          
                          uiOutput("UI_outlier_rnaseq_raw"),
                          br(),
                          
                          # Select experimental group
                          h4(strong(tags$span(
                            "2. Select experimental group",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_rnaseq_raw"),
                          htmlOutput("experimentallevels_rnaseq_raw"),
                          br(),
                          
                          # Filter threshold
                          h4(strong(tags$span(
                            "3. Filtering",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the minimum number of  
                                         counts that is required for n = size of the smallest 
                                         experimental group", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          numericInput(
                            inputId = "countfilter_rnaseq_raw",
                            label = "Minimum number of counts in smallest group size",
                            value = 10),
                          br(),
                          
                          #Pre-processing
                          h4(strong(tags$span(
                            "4. Pre-processing",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to perform the pre-processing!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          actionBttn(inputId = "start_preprocessing_rnaseq_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_preprocessing_rnaseq_raw")
                          
                        ),
                        
                        #********************************************************#
                        #   Output panel
                        #********************************************************#
                        mainPanel(
                          uiOutput("UI_QC_rnaseq_raw")
                        ) #Main panel
                        
               ), #Tab panel
               
               #*****************************************************************#
               # Statistical analysis RNA-seq data
               #*****************************************************************#
               tabPanel("Statistical analysis", value = "panel_statistics_rnaseq_raw", 
                        icon = icon("fas fa-chart-bar"),
                        
                        sidebarPanel(
                          h2(strong("Statistical analysis")),
                          
                          h5("In the statistical analysis step, 
                        you can select which groups to compare to each other and 
                           which covariates to add to the statistical model."),
                          
                          hr(),
                          h4(strong(tags$span(
                            "1. Make comparisons",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_expFactor_rnaseq_raw"),
                          uiOutput("UI_comparisons_rnaseq_raw"),
                          htmlOutput("expFactor_levels_rnaseq_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariated",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Adjust for variables like age and sex 
                                       by adding them as covariates to the model.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_covGroups_num_rnaseq_raw"),
                          uiOutput("UI_covGroups_char_rnaseq_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "3. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Add gene annotations from biomaRt to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          awesomeCheckbox(inputId = "addAnnotation_rnaseq_raw",
                                          label = "Add gene annotations from biomaRt", 
                                          value = FALSE,
                                          status = "danger"),
                          conditionalPanel(
                            condition = "input.addAnnotation_rnaseq_raw==true",
                            uiOutput("UI_biomart_dataset_rnaseq_raw"),
                            uiOutput("UI_addAnnotations_rnaseq_raw")
                          ),
                          br(),
                          
                          h4(strong(tags$span(
                            "4. Perform statistical analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the statistical analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_statistics_rnaseq_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_statistics_rnaseq_raw")
                          
                        ),
                        mainPanel(
                          uiOutput("UI_comparisons_view_rnaseq_raw"),
                          br(),
                          uiOutput("UI_output_statistics_rnaseq_raw")
                        )
               ), # Tab panel
               
               #*****************************************************************#
               # ORA RNA-seq data
               #*****************************************************************#
               tabPanel("ORA", value = "panel_ORA_rnaseq_raw", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Overrepresentation analysis")),
                          
                          h5("In the gene overrepresentation analysis (ORA) step, 
                        you can find processes that are enriched by the 
                           differentially expressed genes."),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Choose from the comparisons for which 
                                       the statistical analysis was performed in the previous step.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_rnaseq_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select geneset collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "A geneset is a collection of genes that are association 
                                       with a specific biological process (GO-BP), molecular function (GO-MF),
                                       cellular component (GO-CC), or biological pathway (WikiPathways).", 
                                         position = "right",
                                         size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_rnaseq_raw",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways")),
                          
                          br(),
                          h4(strong(tags$span(
                            "3. Select differentially expressed genes",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Tip: look at the volcano plot in the previous 
                                       step (statistical analysis) to find the optimal P value and 
                                       logFC thresholds.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          # Up/down regulated genes only
                          prettyRadioButtons(
                            inputId = "updown_ORA_rnaseq_raw",
                            label = "Perform ORA on ...", 
                            choices = 
                              c("Upregulated genes only", 
                                "Downregulated genes only",
                                "Both"),
                            selected = "Both",
                            status = "danger",
                            fill = TRUE),
                          
                          # top N or logFC/P value threshold
                          radioGroupButtons(inputId = "topNorThres_rnaseq_raw",
                                            label = "Select DEGs based on ...", 
                                            choices = c("Threshold", 
                                                        "Top N"),
                                            status = "danger"),
                          
                          # Select DEGs based on logFC/P value threshold
                          conditionalPanel(
                            condition = "input.topNorThres_rnaseq_raw=='Threshold'",
                            
                            #P value threshold
                            numericInput(
                              inputId = "p_thres_ORA_rnaseq_raw",
                              label = "P threshold",
                              value = 0.05),
                            
                            # Raw or adjusted P value?
                            prettyRadioButtons(inputId = "rawp_ORA_rnaseq_raw", 
                                               label = NULL, 
                                               choices = c("Raw P value" = "raw", 
                                                           "Adjusted P value" = "adj"),
                                               inline = TRUE, 
                                               status = "danger",
                                               fill = TRUE),
                            
                            #logFC threshold
                            numericInput(
                              inputId = "logFC_thres_ORA_rnaseq_raw",
                              label = "logFC threshold",
                              value = 0),
                          ),
                          
                          # Select top N most significant genes
                          conditionalPanel(
                            condition = "input.topNorThres_rnaseq_raw=='Top N'",
                            numericInput(
                              inputId = "topN_rnaseq_raw",
                              label = "Top N most significant genes",
                              value = 100),
                          ),
                          
                          br(),
                          
                          # Select gene identifier
                          h4(strong(tags$span(
                            "4. Select gene identifier",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "We need to know which gene identifiers are 
                                       used, so we can link the genes to their correct genesets.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_geneID_ORA_rnaseq_raw"),
                          br(),
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform ORA",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the overrepresentation analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          actionBttn(inputId = "calculate_ORA_rnaseq_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync"))
                          
                        ),
                        mainPanel(
                          uiOutput("UI_output_ORA_rnaseq_raw")
                        )
               ), # Tab panel
               
               ###################################################################
               
               #  RNA-seq: processed
               
               ###################################################################
               
               #*****************************************************************#
               # Upload RNA-seq data
               #*****************************************************************#
               
               tabPanel("Upload", value = "panel_upload_rnaseq_norm", 
                        icon = icon("fas fa-upload"),
                        
                        # Side panel
                        sidebarPanel(
                          
                          #Title
                          h2(strong("Data upload")),
                          h5("Before you can run the analysis workflow, you first 
                           need to upload the expression data as well as the meta data."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a .tsv/.csv file."),
                          fileInput(inputId = "uploadExprData_rnaseq_norm",
                                    label = NULL,
                                    accept = c(".tsv",".csv"),
                                    placeholder = "Select .tsv/.csv file"),
                          
                          h4(strong("2. Upload meta data")),
                          h5("The meta data includes relevant information (e.g., diagnostic group)
                           about the samples. You can upload the meta data as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_rnaseq_norm', 'here'),
                             "for an example .csv meta data file. A Series Matrix File
                           can be downloaded from the GEO website."),
                          
                          prettyRadioButtons(inputId = "MetaFileType_rnaseq_norm", 
                                             label = NULL, 
                                             choices = c(".tsv/.csv file",
                                                         "Series Matrix File"),
                                             inline = TRUE,
                                             fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_rnaseq_norm=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_rnaseq_norm_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file")
                            
                          ),
                          
                          conditionalPanel(
                            condition = "input.MetaFileType_rnaseq_norm=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_rnaseq_norm_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File")
                            
                          ),
                          
                          # Confirm upload
                          actionBttn(inputId = "upload_rnaseq_norm",
                                     label = "Read data",
                                     style = "simple",
                                     color = "primary",
                                     icon = icon("fas fa-upload")),
                          actionBttn(inputId = "example_rnaseq_norm",
                                     label = "Run example",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("cloud-arrow-up")),
                          
                          
                          uiOutput("next_upload_rnaseq_norm")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          uiOutput("UI_upload_rnaseq_norm")
                        )
                        
               ), # End of upload RNA-seq (raw) tab
               
               #*****************************************************************#
               # Pre-processing of RNA-seq data
               #*****************************************************************#
               tabPanel("Pre-processing", value = "panel_preprocessing_rnaseq_norm", icon = icon("sync"),
                        sidebarPanel(
                          
                          #******************************************************#
                          #   Information
                          #******************************************************#
                          
                          h2(strong("Pre-processing")),
                          
                          h5("In this pre-processing step, you can remove samples 
                           (e.g., outliers), perform normalization, 
                           and choose your desired probeset annotation."),
                          
                          hr(),
                          
                          #******************************************************#
                          #   Options
                          #******************************************************#
                          
                          #Remove outliers
                          h4(strong(tags$span(
                            "1. Remove samples",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Uncheck the box to exclude samples from the analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          awesomeCheckbox(inputId = "outlier_rnaseq_norm",
                                          label = "Keep all samples", 
                                          value = TRUE,
                                          status = "danger"),
                          
                          uiOutput("UI_outlier_rnaseq_norm"),
                          br(),
                          
                          # Select experimental group
                          h4(strong(tags$span(
                            "2. Select experimental group",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_rnaseq_norm"),
                          
                          htmlOutput("experimentallevels_rnaseq_norm"),
                          br(),
                          
                          #Transformation
                          h4(strong(tags$span(
                            "3. Transformation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the transformation method. 
                                       The aim of data transformation is to make the 
                                       data more normally distributed.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          uiOutput("UI_transformation_rnaseq_norm"),
                          br(),
                          #Normalization
                          h4(strong(tags$span(
                            "4. Normalization",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all arrays or per experimental group.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          radioGroupButtons(inputId = "normMeth_rnaseq_norm", 
                                            label = NULL, 
                                            choices = c("Quantile","None"),
                                            selected = "None",
                                            status = "danger"),
                          
                          prettyRadioButtons(inputId = "perGroup_rnaseq_norm", 
                                             label = NULL, 
                                             choices = c("Use all arrays",
                                                         "Per experimental group"),
                                             inline = TRUE, 
                                             status = "danger",
                                             fill = TRUE),
                          
                          br(),
                          
                          #Pre-processing
                          h4(strong(tags$span(
                            "5. Pre-processing",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to perform the pre-processing!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          
                          
                          actionBttn(inputId = "start_preprocessing_rnaseq_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_preprocessing_rnaseq_norm")
                          
                        ),
                        
                        #********************************************************#
                        #   Output panel
                        #********************************************************#
                        mainPanel(
                          uiOutput("UI_QC_rnaseq_norm")
                        ) #Main panel
                        
               ), #Tab panel
               
               #*****************************************************************#
               # Statistical analysis RNA-seq data
               #*****************************************************************#
               tabPanel("Statistical analysis", value = "panel_statistics_rnaseq_norm", 
                        icon = icon("fas fa-chart-bar"),
                        
                        sidebarPanel(
                          h2(strong("Statistical analysis")),
                          
                          h5("In the statistical analysis step, 
                        you can select which groups to compare to each other and 
                           which covariates to add to the statistical model."),
                          
                          hr(),
                          h4(strong(tags$span(
                            "1. Make comparisons",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_expFactor_rnaseq_norm"),
                          uiOutput("UI_comparisons_rnaseq_norm"),
                          htmlOutput("expFactor_levels_rnaseq_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariated",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Adjust for variables like age and sex 
                                       by adding them as covariates to the model.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_covGroups_num_rnaseq_norm"),
                          uiOutput("UI_covGroups_char_rnaseq_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "3. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Add gene annotations from biomaRt to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          awesomeCheckbox(inputId = "addAnnotation_rnaseq_norm",
                                          label = "Add gene annotations from biomaRt", 
                                          value = FALSE,
                                          status = "danger"),
                          conditionalPanel(
                            condition = "input.addAnnotation_rnaseq_norm==true",
                            uiOutput("UI_biomart_dataset_rnaseq_norm"),
                            uiOutput("UI_addAnnotations_rnaseq_norm")
                          ),
                          br(),
                          
                          h4(strong(tags$span(
                            "4. Perform statistical analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Click to start the statistical analysis!", 
                                         position = "right",
                                         size = "large")
                          ))),
                          actionBttn(inputId = "calculate_statistics_rnaseq_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          
                          br(),
                          
                          br(),
                          
                          #Go forward
                          uiOutput("UI_next_statistics_rnaseq_norm")
                          
                        ),
                        mainPanel(
                          uiOutput("UI_comparisons_view_rnaseq_norm"),
                          br(),
                          uiOutput("UI_output_statistics_rnaseq_norm")
                        )
               ), # Tab panel
               
               #*****************************************************************#
               # ORA RNA-seq data
               #*****************************************************************#
               tabPanel("ORA", value = "panel_ORA_rnaseq_norm", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Overrepresentation analysis")),
                          
                          h5("In the gene overrepresentation analysis (ORA) step, 
                        you can find processes that are enriched by the 
                           differentially expressed genes."),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Choose from the comparisons for which 
                                       the statistical analysis was performed in the previous step.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_rnaseq_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select geneset collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "A geneset is a collection of genes that are association 
                                       with a specific biological process (GO-BP), molecular function (GO-MF),
                                       cellular component (GO-CC), or biological pathway (WikiPathways).", 
                                         position = "right",
                                         size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_rnaseq_norm",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways")),
                          
                          br(),
                          h4(strong(tags$span(
                            "3. Select differentially expressed genes",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              add_prompt(message = "Tip: look at the volcano plot in the previous 
                                       step (statistical analysis) to find the optimal P value and 
                                       logFC thresholds.", 
                                         position = "right",
                                         size = "large")
                          ))),
                          # Up/down regulated genes only
                          prettyRadioButtons(
                            inputId = "updown_ORA_rnaseq_norm",
                            label = "Perform ORA on ...", 
                            choices = 
                              c("Upregulated genes only", 
                                "Downregulated genes only",
                                "Both"),
                            selected = "Both",
                            status = "danger",
                            fill = TRUE),
                          
                          # top N or logFC/P value threshold
                          radioGroupButtons(inputId = "topNorThres_rnaseq_norm",
                                            label = "Select DEGs based on ...", 
                                            choices = c("Threshold", 
                                                        "Top N"),
                                            status = "danger"),
                          
                          # Select DEGs based on logFC/P value threshold
                          conditionalPanel(
                            condition = "input.topNorThres_rnaseq_norm=='Threshold'",
                            
                            #P value threshold
                            numericInput(
                              inputId = "p_thres_ORA_rnaseq_norm",
                              label = "P threshold",
                              value = 0.05),
                            
                            # Raw or adjusted P value?
                            prettyRadioButtons(inputId = "rawp_ORA_rnaseq_norm", 
                                               label = NULL, 
                                               choices = c("Raw P value" = "raw", 
                                                           "Adjusted P value" = "adj"),
                                               inline = TRUE, 
                                               status = "danger",
                                               fill = TRUE),
                            
                            #logFC threshold
                            numericInput(
                              inputId = "logFC_thres_ORA_rnaseq_norm",
                              label = "logFC threshold",
                              value = 0),
                          ),
                          
                          # Select top N most significant genes
                          conditionalPanel(
                            condition = "input.topNorThres_rnaseq_norm=='Top N'",
                            numericInput(
                              inputId = "topN_rnaseq_norm",
                              label = "Top N most significant genes",
                              value = 100),
                          ),
                          
                          br(),
                          
                          # Select gene identifier
                          h4(strong(tags$span(
                            "4. Select gene identifier",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "We need to know which gene identifiers are 
                                       used, so we can link the genes to their correct genesets.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_geneID_ORA_rnaseq_norm"),
                          br(),
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform ORA",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the overrepresentation analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          actionBttn(inputId = "calculate_ORA_rnaseq_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync"))
                          
                        ),
                        mainPanel(
                          uiOutput("UI_output_ORA_rnaseq_norm")
                        )
               ), # Tab panel
               
               ###################################################################
               
               #  Documentation
               
               ###################################################################
               tabPanel("Documentation", value = "documentation", 
                        icon = icon("question-circle")
               ) # Tab panel
               
               
               ###################################################################
               
               #  END
               
               ###################################################################
    ) # EO navbarPage
  ) # EO fluidPage
) # EO tagList

