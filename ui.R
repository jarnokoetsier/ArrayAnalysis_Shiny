#==============================================================================#
# Name: ui.R
# Description: User interface of the ArrayAnalysis Shiny app
#==============================================================================#

# Start UI
ui <- tagList(
  
  # Set the style of the UI
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none;
                           width: 100%;
                           margin: 0;
                           padding: 0;
                           }
                            .navbar-nav > li:nth-child(18) {
                           float: right;
                            }
                           
                           .my_style_1{ 
                            background-image: url(background.jpg);
                            background-size: cover;
                            background-repeat: no-repeat;
                            background-attachment: fixed;
                            background-position: center center;
                            min-height: 100vh;
                            width: 100%;
                            margin-top: -20px; 
                            margin-right: 0; 
                            padding: 0;
                            display: flex;
                            align-items: center;
                           }
                           
                           .container-fluid { 
                           padding-left: 0; 
                           padding-right: 0;
                           }
                           "))
  ),
  
  tags$head(tags$style(HTML("
    .pretty input:checked~.state.p-success label:after, .pretty.p-toggle .state.p-success label:after {
    background-color: black!important;
}"))),
  
  tags$head(tags$style(HTML("
    .dropbtn {
       background-color: #D9D9D9;
      color: white;
      padding: 0;
      display: flex;
      font-size: 28px;
      border: none;
      border-radius: 50%;
      cursor: pointer;
      width: 50px;
      height: 50px;
      align-items: center;
      justify-content: center;
      line-height = 1;
    }

    .dropdown {
      position: relative;
      display: inline-block;
    }

    .dropdown-content {
      display: none;
      position: absolute;
      background-color: #f9f9f9;
      min-width: 250px;
      box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
      padding: 12px;
      z-index: 1;
      border-radius: 4px;
    }

    .dropdown:hover .dropdown-content {
      display: block;
    }
  "))),
  
  fluidPage(
    
    # This allow for pop-up messages to show up
    shinyWidgets::useSweetAlert(),
    
    # This allows for insertion of information boxes in the app
    prompter::use_prompt(),
    
    # Make a page with a navigation bar
    navbarPage(title = div(img(src="logo_navbar.PNG", 
                               style="margin-top: -14px;
                               padding-right:10px;
                               padding-bottom:10px",
                               height = 60)),
               windowTitle = "ArrayAnalysis",
               id = "navbar",
               
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
                                 style = "background-color:rgba(255, 255, 255, 0.95); border-radius: 15px;",
                                 
                                 # Line break
                                 br(),
                                 
                                 # ArrayAnalysis logo
                                 img(src = "logo_main.PNG", width = "100%"),
                                 h5(ArrayAnalysis_version, align = "center"),
                                 
                                 # Welcome message
                                 h1(strong(span(style = "color:#000000", 
                                                "Welcome to ArrayAnalysis!"))),
                                 h4(span(style = "color:#000000", 
                                         "Do you want to analyze microarray or RNA-Seq data?")),
                                 
                                 # Line break
                                 br(),
                                 
                                 # Analyse microarray or RNA-seq data?
                                 selectInput(
                                   inputId = "microarray_or_rnaseq",
                                   label = NULL,
                                   choices = c("RNA-Seq analysis" = "RNA-Seq",
                                               "Microarray analysis" = "Microarray"),
                                   selected = "RNA-Seq",
                                 ),
                                 
                                 # shinyWidgets::radioGroupButtons(
                                 #   inputId = "microarray_or_rnaseq",
                                 #   label = NULL,
                                 #   choices = c("Microarray analysis" = "Microarray",
                                 #               "RNA-Seq analysis" = "RNA-Seq"),
                                 #   status = "info",
                                 #   selected = "RNA-Seq"
                                 # ),
                                 
                                 conditionalPanel(
                                   condition = "input.microarray_or_rnaseq=='RNA-Seq'",
                                   # Analyse raw or processed data?
                                   shinyWidgets::prettyRadioButtons(
                                     inputId = "raw_or_norm_rnaseq",
                                     label = NULL,
                                     choices = c("Raw counts" = "Raw data",
                                                 "Processed counts" = "Processed data"),
                                     inline = TRUE,
                                     status = "info",
                                     fill = TRUE),
                                 ),
                                 
                                 
                                 conditionalPanel(
                                   condition = "input.microarray_or_rnaseq=='Microarray'",
                                   # Analyse raw or processed data?
                                   shinyWidgets::prettyRadioButtons(
                                     inputId = "raw_or_norm_microarray",
                                     label = NULL,
                                     choices = c("CEL files" = "Raw data",
                                                 "Processed intensities" = "Processed data"),
                                     inline = TRUE,
                                     status = "info",
                                     fill = TRUE),
                                 ),
                                 
                                 # Line break
                                 br(),
                                 
                                 shinyWidgets::actionBttn(inputId = "advancedSettings",
                                                          label = NULL,
                                                          style = "simple",
                                                          color = "warning",
                                                          icon = icon("circle-question")),
                                 
                                 # Action button: start analysis by clicking
                                 shinyWidgets::actionBttn(inputId = "startAnalysis",
                                                          label = "Start Analysis",
                                                          style = "simple",
                                                          color = "warning",
                                                          icon = icon("arrow-right")),
                                 
                                 # Line breaks
                                 br(),
                                 br(),
                                 uiOutput("versionMessage"),
                                 br()
                                 
                          ) # EO column
                        ), # EO fluidRow
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
                          h5("Before you can run the analysis, you first 
                           need to upload the expression data and metadata."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as an ", 
                             em(".zip"), " folder containing all ", em(".CEL / .CEL.gz"), 
                             " files. The file names should match with
                           the sample IDs in the metadata table.",
                             "Click ", downloadLink('downloadexpr_example_microarray_raw', 
                                                    'here'),
                             "for an example zip file."),
                          
                          # File input for upload expression data
                          fileInput(inputId = "uploadCEL_microarray_raw",
                                    label = NULL,
                                    accept = ".zip",
                                    placeholder = "Select .zip data file"),
                          
                          h4(strong("2. Upload metadata")),
                          h5("The metadata includes relevant information
                           about the samples (e.g., genotype or experimental group). 
                           You can upload the metadata as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_microarray_raw', 
                                                    'here'),
                             "for an example .csv metadata file. A Series Matrix File
                           can be downloaded from the",
                             a("GEO website.", 
                               href = "https://www.ncbi.nlm.nih.gov/geo/",
                               target="_blank")),
                          
                          # Select .csv/.tsv or Series Matrix File
                          shinyWidgets::prettyRadioButtons(inputId = "MetaFileType_microarray_raw", 
                                                           label = NULL, 
                                                           choices = c(".tsv/.csv file",
                                                                       "Series Matrix File"),
                                                           inline = TRUE,
                                                           fill = TRUE),
                          
                          # .tsv/.csv file input for upload metadata
                          conditionalPanel(
                            condition = "input.MetaFileType_microarray_raw=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_microarray_raw_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file"),
                            
                          ),
                          
                          # Series Matrix File input for upload metadata
                          conditionalPanel(
                            condition = "input.MetaFileType_microarray_raw=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_microarray_raw_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File"),
                            
                          ),
                          
                          # Confirm upload
                          shinyWidgets::actionBttn(inputId = "upload_microarray_raw",
                                                   label = "Read data",
                                                   style = "simple",
                                                   color = "primary",
                                                   icon = icon("fas fa-upload")),
                          
                          # Run example
                          shinyWidgets::actionBttn(inputId = "example_microarray_raw",
                                                   label = "Run example",
                                                   style = "simple",
                                                   color = "warning",
                                                   icon = icon("cloud-arrow-up")),
                          
                          # Button to go to next tab
                          uiOutput("next_upload_microarray_raw")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          
                          # Expression matrix and metadata table
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
                              prompter::add_prompt(message = "Uncheck the box to exclude samples from the analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          shinyWidgets::awesomeCheckbox(inputId = "outlier_microarray_raw",
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
                              prompter::add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare to each other (e.g., disease vs healthy). 
                                                   The selected groups will only be used in the pre-processing in case
                                                   normalization 'Per experimental group' is selected. The selected groups
                                                   will also be used for the visualization of the pre-processing quality.", 
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
                              prompter::add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all arrays or per experimental group.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          # shinyWidgets::radioGroupButtons(
                          #   inputId = "normMeth_microarray_raw",
                          #   label = NULL,
                          #   choices = c("RMA","GCRMA","PLIER"), #"None"
                          #   status = "danger"),
                          
                          selectInput(
                            inputId = "normMeth_microarray_raw",
                            label = NULL,
                            choices = c("RMA","GCRMA","PLIER")),

                          shinyWidgets::prettyRadioButtons(
                            inputId = "perGroup_microarray_raw", 
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
                              prompter::add_prompt(message = "Select the probeset annotation. 
                                       ENTREZG or ENSG custom annotations are recommended.
                                       If no annotation is selected, the standard affy annotations 
                                       will be used by default.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          # selectInput(inputId = "annotations_microarray_raw",
                          #             label = NULL,
                          #             choices = c("No annotations",
                          #                         "Custom annotations",
                          #                         "Upload annotation file"),
                          #             selected = "Custom annotations"),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "annotations_microarray_raw", 
                            label = NULL, 
                            choices = c("Custom" = "Custom annotations",
                                        "Upload" = "Upload annotation file",
                                        "None" = "No annotations"),
                            selected = "Custom annotations",
                            inline = TRUE, 
                            status = "danger",
                            fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.annotations_microarray_raw=='Custom annotations'",
                            
                            uiOutput("UI_species_microarray_raw"),
                            
                            selectInput(inputId = "CDFtype_microarray_raw",
                                        label = "Annotation format",
                                        choices = c("ENTREZG","REFSEQ","ENSG",
                                                    "ENSE","ENST","VEGAG","VEGAE",
                                                    "VEGAT","TAIRG","TAIRT","UG",
                                                    "MIRBASEF","MIRBASEG"))
                          ),
                          conditionalPanel(
                            condition = "input.annotations_microarray_raw=='Upload annotation file'",
                            
                            fileInput(inputId = "annot_file_microarray_raw",
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
                              prompter::add_prompt(message = "Click to perform the pre-processing!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          shinyWidgets::actionBttn(inputId = "start_preprocessing_microarray_raw",
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
                        
                        #******************************************************#
                        #   Information
                        #******************************************************#
                        
                        sidebarPanel(
                          h2(strong("Statistical analysis")),
                          h5("In the statistical analysis step, 
                        you can find differentially expressed genes by selecting 
                        which groups to compare to each other and 
                           which covariates to add to the statistical model."),
                          hr(),
                          
                          #******************************************************#
                          #   Options
                          #******************************************************#
                          
                          # Make comparisons
                          h4(strong(tags$span(
                            "1. Make comparisons",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other. Comparisons are commonly defined 
                                                   as `Disease - Control` or `Treatment - Control`.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          # Select experimental factor
                          uiOutput("UI_expFactor_microarray_raw"),
                          
                          # Select comparison(s)
                          uiOutput("UI_comparisons_microarray_raw"),
                          htmlOutput("expFactor_levels_micorarray_raw"),
                          br(),
                          
                          # Select covariates
                          h4(strong(tags$span(
                            "2. Add covariates",
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
                          
                          # Continuous/numeric covariates
                          uiOutput("UI_covGroups_num_microarray_raw"),
                          
                          # Descrete/character covariates
                          uiOutput("UI_covGroups_char_microarray_raw"),
                          br(),
                          
                          # Add gene annotations to top table using biomaRt
                          h4(strong(tags$span(
                            "3. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Add gene annotations to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::materialSwitch(inputId = "addAnnotation_microarray_raw",
                                                       label = "Add gene annotations", 
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
                              prompter::add_prompt(message = "Click to start the statistical analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_statistics_microarray_raw",
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
               # Gene set analysis
               #*****************************************************************#
               tabPanel("Gene set analysis", value = "panel_ORA_microarray_raw", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Gene set analysis")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "ORA_or_GSEA_microarray_raw",
                            label = NULL, 
                            choices = c("ORA", "GSEA"),
                            selected = "ORA",
                            status = "success",
                            fill = TRUE,
                            inline = TRUE),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw == 'ORA'",
                            h5("With",strong("Overrepresentation Analysis (ORA),"),"dysregulated 
                          processes and pathways can be idenified. These processes/pathways 
                          are identified by testing whether their genes are overrepresented among the (most) significant genes."),
                          ),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw == 'GSEA'",
                            h5("With", strong("Gene Set Enrichment Analysis (GSEA),"), "you can find dysregulated 
                          processes and pathways. These processes/pathways are identified by testing whether 
                               their genes show concordant changes in the data."),
                          ),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Choose for which comparison you want to perform gene set analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_microarray_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select gene set collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "A gene set collection contains different sets of genes that are linked to 
                                       biological processes (GO-BP), molecular functions (GO-MF),
                                       cellular components (GO-CC), or biological pathways (WikiPathways and KEGG).", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_microarray_raw",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways",
                                                  "KEGG")),
                          
                          br(),
                          
                          # ORA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw=='ORA'",
                            h4(strong(tags$span(
                              "3. Select genes",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select which genes are used in the analysis.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Up/down regulated genes only
                            shinyWidgets::prettyRadioButtons(
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
                            shinyWidgets::prettyRadioButtons(
                              inputId = "topNorThres_microarray_raw",
                              label = "Select genes based on ...", 
                              choices = c("Threshold", 
                                          "Top N"),
                              status = "danger",
                              selected = "Top N",
                              fill = TRUE,
                              inline = TRUE),
                            
                            # Select genes based on logFC/P value threshold
                            conditionalPanel(
                              condition = "input.topNorThres_microarray_raw=='Threshold'",
                              
                              #P value threshold
                              numericInput(
                                inputId = "p_thres_ORA_microarray_raw",
                                label = "P threshold",
                                value = 0.05),
                              
                              # Raw or adjusted P value?
                              shinyWidgets::prettyRadioButtons(
                                inputId = "rawp_ORA_microarray_raw", 
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
                            )
                            
                          ),
                          
                          # GSEA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw=='GSEA'",
                            h4(strong(tags$span(
                              "3. Select ranking variable",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select the variable on which the GSEA should be based.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Ranking variable
                            shinyWidgets::prettyRadioButtons(
                              inputId = "ranking_GSEA_microarray_raw",
                              label = NULL, 
                              choices = 
                                c(" logFC" = "logFC",
                                  "-log p-value" = "pvalue",
                                  "-log p-value x sign logFC" = "signed_pvalue"),
                              selected = "logFC",
                              status = "danger",
                              fill = TRUE),
                            
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
                              prompter::add_prompt(message = "To link genes to gene sets, 
                              we need to know which gene identifiers are used.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          # Select organism
                          selectInput(inputId = "organism_ORA_microarray_raw",
                                      label = tags$span(
                                        "Organism", 
                                        tags$span(
                                          icon(
                                            name = "question-circle",
                                          ) 
                                        ) |>
                                          prompter::add_prompt(message = "Select the organism. 
                                               This information is needed to match the gene IDs.", 
                                                               position = "right",
                                                               size = "large")
                                      ),
                                      choices = c("Bos taurus",
                                                  "Caenorhabditis elegans",
                                                  "Homo sapiens",
                                                  "Mus musculus", 
                                                  "Rattus norvegicus"),
                                      selected = "Homo sapiens"),
                          
                          # Which columns of the top table contains the gene ids?
                          uiOutput("UI_geneID_ORA_microarray_raw"),
                          
                          br(),
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the gene set analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_ORA_microarray_raw",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          br(),
                          hr(),
                          uiOutput("UI_ORAreport_microarray_raw"),
                          uiOutput("UI_GSEAreport_microarray_raw")
                          
                        ),
                        mainPanel(
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw=='ORA'",
                            uiOutput("UI_output_ORA_microarray_raw")
                          ),
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_raw=='GSEA'",
                            uiOutput("UI_output_GSEA_microarray_raw")
                          )
                          
                          
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
                          h5("Before you can run the analysis, you first 
                           need to upload the expression data and metadata."),
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a .tsv/.csv file or as a Series 
                           Matrix File.", "Click ", 
                             downloadLink('downloadexpr_example_microarray_norm', 
                                                    'here'),
                             "for an example expression matrix in csv format.",
                             "A Series Matrix File can be downloaded from the ",
                             a("GEO website.", 
                               href = "https://www.ncbi.nlm.nih.gov/geo/",
                               target="_blank")),
                          
                          # .tsv/.csv or Series Matrix File
                          shinyWidgets::prettyRadioButtons(
                            inputId = "ExprDataFileType_microarray_norm", 
                            label = NULL, 
                            choices = c(".tsv/.csv file",
                                        "Series Matrix File"),
                            inline = TRUE,
                            fill = TRUE),
                          
                          conditionalPanel(
                            condition = "input.ExprDataFileType_microarray_norm=='.tsv/.csv file'",
                            
                            # File input for .tsv/.csv file
                            fileInput(inputId = "uploadExprData_microarray_norm_tsv",
                                      label = NULL,
                                      accept = c(".tsv", ".csv"),
                                      placeholder = "Select .tsv or .csv data file"),
                            
                          ),
                          
                          conditionalPanel(
                            condition = "input.ExprDataFileType_microarray_norm=='Series Matrix File'",
                            
                            # File input for Series Matrix File
                            fileInput(inputId = "uploadExprData_microarray_norm_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File"),
                          ),
                          
                          
                          h4(strong("2. Upload metadata")),
                          h5("The metadata includes relevant information
                           about the samples (e.g., genotype or experimental group). 
                           You can upload the metadata as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_microarray_norm', 
                                                    'here'),
                             "for an example .csv metadata file. A Series Matrix File
                           can be downloaded from the",
                             a("GEO website.", 
                               href = "https://www.ncbi.nlm.nih.gov/geo/",
                               target="_blank")),
                          
                          # .tsv/.csv or Series Matrix File
                          shinyWidgets::prettyRadioButtons(
                            inputId = "MetaFileType_microarray_norm", 
                            label = NULL, 
                            choices = c(".tsv/.csv file",
                                        "Series Matrix File"),
                            inline = TRUE,
                            fill = TRUE),
                          
                          # File input for .tsv/.csv file
                          conditionalPanel(
                            condition = "input.MetaFileType_microarray_norm=='.tsv/.csv file'",
                            fileInput(inputId = "uploadMeta_microarray_norm_tsv",
                                      label = NULL,
                                      accept = c(".tsv",".csv"),
                                      placeholder = "Select .tsv or .csv data file")
                            
                          ),
                          
                          # File input for Series Matrix File
                          conditionalPanel(
                            condition = "input.MetaFileType_microarray_norm=='Series Matrix File'",
                            fileInput(inputId = "uploadMeta_microarray_norm_smf",
                                      label = NULL,
                                      accept = ".txt.gz",
                                      placeholder = "Select .txt.gz Series Matrix File")
                            
                          ),
                          
                          # Confirm upload
                          shinyWidgets::actionBttn(inputId = "upload_microarray_norm",
                                                   label = "Read data",
                                                   style = "simple",
                                                   color = "primary",
                                                   icon = icon("fas fa-upload")),
                          
                          # Run example
                          shinyWidgets::actionBttn(inputId = "example_microarray_norm",
                                                   label = "Run example",
                                                   style = "simple",
                                                   color = "warning",
                                                   icon = icon("cloud-arrow-up")),
                          
                          # Button to go to next tab
                          uiOutput("next_upload_microarray_norm")
                          
                        ), # End of side panel
                        
                        # Main panel
                        mainPanel(
                          
                          # Preview of expression matrix and metadata table
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
                           (e.g., outliers) and perform transformation and normalization."),
                          
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
                          
                          shinyWidgets::awesomeCheckbox(inputId = "outlier_norm",
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
                              prompter::add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_microarray_norm"),
                          
                          htmlOutput("experimentallevels_microarray_norm"),
                          br(),
                          
                          #Transformation
                          h4(strong(tags$span(
                            "3. Transformation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Select the transformation method. 
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
                              prompter::add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all arrays or per experimental group.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          selectInput(
                            inputId = "normMeth_microarray_norm", 
                            label = NULL, 
                            choices = c("Quantile normalization" = "Quantile",
                                        "Continue without normalization" = "None"),
                            selected = "None"),
                          
                          conditionalPanel(
                            condition = "input.normMeth_microarray_norm == `Quantile`",
                            shinyWidgets::prettyRadioButtons(
                              inputId = "perGroup_microarray_norm", 
                              label = NULL, 
                              choices = c("Use all arrays",
                                          "Per experimental group"),
                              inline = TRUE, 
                              status = "danger",
                              fill = TRUE)
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
                              prompter::add_prompt(message = "Click to perform the pre-processing!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          shinyWidgets::actionBttn(inputId = "start_preprocessing_microarray_norm",
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
                              prompter::add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other. Comparisons are commonly defined 
                                                   as `Disease - Control` or `Treatment - Control`.", 
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
                              prompter::add_prompt(message = "Adjust for variables like age and sex 
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
                              prompter::add_prompt(message = "Add gene annotations to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::materialSwitch(inputId = "addAnnotation_microarray_norm",
                                                       label = "Add gene annotations", 
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
                              prompter::add_prompt(message = "Click to start the statistical analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_statistics_microarray_norm",
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
               # Gene set analysis
               #*****************************************************************#
               tabPanel("Gene set analysis", value = "panel_ORA_microarray_norm", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Gene set analysis")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "ORA_or_GSEA_microarray_norm",
                            label = NULL, 
                            choices = c("ORA", "GSEA"),
                            selected = "ORA",
                            status = "success",
                            fill = TRUE,
                            inline = TRUE),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm == 'ORA'",
                            h5("With",strong("Overrepresentation Analysis (ORA),"),"dysregulated 
                          processes and pathways can be idenified. These processes/pathways 
                          are identified by testing whether their genes are overrepresented among the (most) significant genes."),
                          ),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm == 'GSEA'",
                            h5("With", strong("Gene Set Enrichment Analysis (GSEA),"), "you can find dysregulated 
                          processes and pathways. These processes/pathways are identified by testing whether 
                               their genes show concordant changes in the data."),
                          ),
                        
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Choose for which comparison you want to perform gene set analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_microarray_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select gene set collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "A gene set collection contains different sets of genes that are linked to 
                                       biological processes (GO-BP), molecular functions (GO-MF),
                                       cellular components (GO-CC), or biological pathways (WikiPathways and KEGG).", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_microarray_norm",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways",
                                                  "KEGG")),
                          
                          br(),
                          
                          # ORA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm=='ORA'",
                            h4(strong(tags$span(
                              "3. Select genes",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select which genes are used in the analysis.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Up/down regulated genes only
                            shinyWidgets::prettyRadioButtons(
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
                            shinyWidgets::prettyRadioButtons(
                              inputId = "topNorThres_microarray_norm",
                              label = "Select genes based on ...", 
                              choices = c("Threshold", 
                                          "Top N"),
                              status = "danger",
                              selected = "Top N",
                              fill = TRUE,
                              inline = TRUE),
                            
                            # Select genes based on logFC/P value threshold
                            conditionalPanel(
                              condition = "input.topNorThres_microarray_norm=='Threshold'",
                              
                              #P value threshold
                              numericInput(
                                inputId = "p_thres_ORA_microarray_norm",
                                label = "P threshold",
                                value = 0.05),
                              
                              # Raw or adjusted P value?
                              shinyWidgets::prettyRadioButtons(inputId = "rawp_ORA_microarray_norm", 
                                                               label = NULL, 
                                                               choices = c("Raw P value" = "raw", 
                                                                           "Adjusted P value" = "adj"),
                                                               inline = TRUE, 
                                                               status = "danger",
                                                               fill = TRUE),
                              
                              #logFC threshold
                              numericInput(
                                inputId = "logFC_thres_ORA_microarray_norm",
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
                            )
                            
                          ),
                          
                          # GSEA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm=='GSEA'",
                            h4(strong(tags$span(
                              "3. Select ranking variable",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select the variable on which the GSEA should be based.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Ranking variable
                            shinyWidgets::prettyRadioButtons(
                              inputId = "ranking_GSEA_microarray_norm",
                              label = NULL, 
                              choices = 
                                c(" logFC" = "logFC",
                                  "-log p-value" = "pvalue",
                                  "-log p-value x sign logFC" = "signed_pvalue"),
                              selected = "logFC",
                              status = "danger",
                              fill = TRUE),
                            
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
                              prompter::add_prompt(message = "To link genes to gene sets, 
                              we need to know which gene identifiers are used.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          # Select organism
                          selectInput(inputId = "organism_ORA_microarray_norm",
                                      label = tags$span(
                                        "Organism", 
                                        tags$span(
                                          icon(
                                            name = "question-circle",
                                          ) 
                                        ) |>
                                          prompter::add_prompt(message = "Select the organism. 
                                               This information is needed to match the gene IDs.", 
                                                               position = "right",
                                                               size = "large")
                                      ),
                                      choices = c("Bos taurus",
                                                  "Caenorhabditis elegans",
                                                  "Homo sapiens",
                                                  "Mus musculus", 
                                                  "Rattus norvegicus"),
                                      selected = "Homo sapiens"),
                          
                          # Which columns of the top table contains the gene ids?
                          uiOutput("UI_geneID_ORA_microarray_norm"),
                          
                          br(),
                          
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the gene set analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_ORA_microarray_norm",
                                     label = "Calculate",
                                     style = "simple",
                                     color = "warning",
                                     icon = icon("sync")),
                          br(),
                          hr(),
                          uiOutput("UI_ORAreport_microarray_norm"),
                          uiOutput("UI_GSEAreport_microarray_norm")
                          
                        ),
                        mainPanel(
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm=='ORA'",
                            uiOutput("UI_output_ORA_microarray_norm")
                          ),
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_microarray_norm=='GSEA'",
                            uiOutput("UI_output_GSEA_microarray_norm")
                          )
                          
                          
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
                          h5("Before you can run the analysis, you first 
                           need to upload the expression data and metadata."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a .tsv/.csv file. 
                             The expression data is a matrix with the genes in the rows  
                             and the samples in the columns. ",
                             "Click ", downloadLink('downloadexpr_example_rnaseq_raw', 
                                                    'here'),
                             "for an example expression matrix."),
                          fileInput(inputId = "uploadExprData_rnaseq_raw",
                                    label = NULL,
                                    accept = c(".tsv",".csv"),
                                    placeholder = "Select .tsv/.csv file"),
                          
                          h4(strong("2. Upload metadata")),
                          h5("The metadata includes relevant information
                           about the samples (e.g., genotype or experimental group). 
                           You can upload the metadata as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_rnaseq_raw', 
                                                    'here'),
                             "for an example .csv metadata file. A Series Matrix File
                           can be downloaded from the",
                             a("GEO website.", 
                               href = "https://www.ncbi.nlm.nih.gov/geo/",
                               target="_blank")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "MetaFileType_rnaseq_raw", 
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
                          shinyWidgets::actionBttn(inputId = "upload_rnaseq_raw",
                                     label = "Read data",
                                     style = "simple",
                                     color = "primary",
                                     icon = icon("fas fa-upload")),
                          shinyWidgets::actionBttn(inputId = "example_rnaseq_raw",
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
                           (e.g., outliers) and perform gene filtering and normalization."),
                          
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
                          
                          shinyWidgets::awesomeCheckbox(inputId = "outlier_rnaseq_raw",
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
                              prompter::add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare to each other 
                                       in the statistical analysis (e.g., disease or treatment status). ", 
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
                              prompter::add_prompt(message = "Select the filtering threshold.
                              Genes with at least this number of (raw) counts 
                              in the smallest group size will be kept for the statistical analysis. 
                                         A value of 10 is often recommended for bulk RNA-seq.", 
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
                              prompter::add_prompt(message = "Click to perform pre-processing!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          shinyWidgets::actionBttn(inputId = "start_preprocessing_rnaseq_raw",
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
                              prompter::add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other. Comparisons are commonly defined 
                                                   as `Disease - Control` or `Treatment - Control`.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_rnaseq_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariates",
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
                            "3. Perform logFC shrinkage",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "If selected, imprecise logFCs are shrunk to 0 with the apeglm method. 
                                                   logFC shrinkage is recommended by DESeq2.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::materialSwitch(inputId = "shrinkage_rnaseq_raw",
                                                       label = NULL, 
                                                       value = TRUE,
                                                       status = "danger"),
                          br(),
                          
                          h4(strong(tags$span(
                            "4. Add gene annotation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Add gene annotations 
                              (e.g., Ensembl gene IDs, Entrez gene IDs, and/or HGNC symbols) to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::materialSwitch(inputId = "addAnnotation_rnaseq_raw",
                                                       label = NULL, 
                                                       value = FALSE,
                                                       status = "danger"),
                          conditionalPanel(
                            condition = "input.addAnnotation_rnaseq_raw==true",
                            uiOutput("UI_biomart_dataset_rnaseq_raw"),
                            uiOutput("UI_addAnnotations_rnaseq_raw")
                          ),
                          br(),
                          
                          h4(strong(tags$span(
                            "5. Perform statistical analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the statistical analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_statistics_rnaseq_raw",
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
               # Gene set analysis RNA-seq data
               #*****************************************************************#
               tabPanel("Gene set analysis", value = "panel_ORA_rnaseq_raw", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Gene set analysis")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "ORA_or_GSEA_rnaseq_raw",
                            label = NULL, 
                            choices = c("ORA", "GSEA"),
                            selected = "ORA",
                            status = "success",
                            fill = TRUE,
                            inline = TRUE),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw == 'ORA'",
                            h5("With",strong("Overrepresentation Analysis (ORA),"),"dysregulated 
                          processes and pathways can be idenified. These processes/pathways 
                          are identified by testing whether their genes are overrepresented among the (most) significant genes."),
                          ),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw == 'GSEA'",
                            h5("With", strong("Gene Set Enrichment Analysis (GSEA),"), "you can find dysregulated 
                          processes and pathways. These processes/pathways are identified by testing whether 
                               their genes show concordant changes in the data."),
                          ),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Choose for which comparison you want to perform gene set analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_rnaseq_raw"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select gene set collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "A gene set collection contains different sets of genes that are linked to 
                                       biological processes (GO-BP), molecular functions (GO-MF),
                                       cellular components (GO-CC), or biological pathways (WikiPathways and KEGG).", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_rnaseq_raw",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways",
                                                  "KEGG")),
                          
                          br(),
                          
                          # ORA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw=='ORA'",
                            h4(strong(tags$span(
                              "3. Select genes",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select which genes are used in the analysis.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Up/down regulated genes only
                            shinyWidgets::prettyRadioButtons(
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
                            shinyWidgets::prettyRadioButtons(
                              inputId = "topNorThres_rnaseq_raw",
                              label = "Select genes based on ...", 
                              choices = c("Threshold", 
                                          "Top N"),
                              status = "danger",
                              selected = "Top N",
                              fill = TRUE,
                              inline = TRUE),
                            
                            # Select genes based on logFC/P value threshold
                            conditionalPanel(
                              condition = "input.topNorThres_rnaseq_raw=='Threshold'",
                              
                              #P value threshold
                              numericInput(
                                inputId = "p_thres_ORA_rnaseq_raw",
                                label = "P threshold",
                                value = 0.05),
                              
                              # Raw or adjusted P value?
                              shinyWidgets::prettyRadioButtons(
                                inputId = "rawp_ORA_rnaseq_raw", 
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
                            )
                            
                          ),
                          
                          # GSEA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw=='GSEA'",
                            h4(strong(tags$span(
                              "3. Select ranking variable",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select the variable on which the GSEA should be based.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Ranking variable
                            shinyWidgets::prettyRadioButtons(
                              inputId = "ranking_GSEA_rnaseq_raw",
                              label = NULL, 
                              choices = 
                                c(" logFC" = "logFC",
                                  "-log p-value" = "pvalue",
                                  "-log p-value x sign logFC" = "signed_pvalue"),
                              selected = "logFC",
                              status = "danger",
                              fill = TRUE),
                            
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
                              prompter::add_prompt(message = "To link genes to gene sets, 
                              we need to know which gene identifiers are used.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          # Select organism
                          selectInput(inputId = "organism_ORA_rnaseq_raw",
                                      label = tags$span(
                                        "Organism", 
                                        tags$span(
                                          icon(
                                            name = "question-circle",
                                          ) 
                                        ) |>
                                          prompter::add_prompt(message = "Select the organism. 
                                               This information is needed to match the gene IDs.", 
                                                               position = "right",
                                                               size = "large")
                                      ),
                                      choices = c("Bos taurus",
                                                  "Caenorhabditis elegans",
                                                  "Homo sapiens",
                                                  "Mus musculus", 
                                                  "Rattus norvegicus"),
                                      selected = "Homo sapiens"),
                          
                          # Which columns of the top table contains the gene ids?
                          uiOutput("UI_geneID_ORA_rnaseq_raw"),
                          
                          br(),
                          
                          # Calculate!
                          h4(strong(tags$span(
                            "5. Perform analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the gene set analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_ORA_rnaseq_raw",
                                                   label = "Calculate",
                                                   style = "simple",
                                                   color = "warning",
                                                   icon = icon("sync")),
                          br(),
                          hr(),
                          uiOutput("UI_ORAreport_rnaseq_raw"),
                          uiOutput("UI_GSEAreport_rnaseq_raw")
                          
                        ),
                        mainPanel(
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw=='ORA'",
                            uiOutput("UI_output_ORA_rnaseq_raw")
                          ),
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_raw=='GSEA'",
                            uiOutput("UI_output_GSEA_rnaseq_raw")
                          )
                          
                          
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
                          h5("Before you can run the analysis, you first 
                           need to upload the expression data and metadata."),
                          
                          hr(),
                          
                          h4(strong("1. Upload expression data")),
                          h5("The expression data should be supplied as a .tsv/.csv file. 
                             The expression data is a matrix with the genes in the rows  
                             and the samples in the columns. ",
                             "Click ", downloadLink('downloadexpr_example_rnaseq_norm', 
                                                    'here'),
                             "for an example expression matrix."),
                          fileInput(inputId = "uploadExprData_rnaseq_norm",
                                    label = NULL,
                                    accept = c(".tsv",".csv"),
                                    placeholder = "Select .tsv/.csv file"),
                          
                          h4(strong("2. Upload metadata")),
                          h5("The metadata includes relevant information
                           about the samples (e.g., genotype or experimental group). 
                           You can upload the metadata as a 
                           .csv/.tsv file or upload a Series Matrix file.", 
                             "Click ", downloadLink('downloadmeta_example_rnaseq_norm', 
                                                    'here'),
                             "for an example .csv metadata file. A Series Matrix File
                           can be downloaded from the",
                             a("GEO website.", 
                               href = "https://www.ncbi.nlm.nih.gov/geo/",
                               target="_blank")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "MetaFileType_rnaseq_norm", 
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
                          shinyWidgets::actionBttn(inputId = "upload_rnaseq_norm",
                                                   label = "Read data",
                                                   style = "simple",
                                                   color = "primary",
                                                   icon = icon("fas fa-upload")),
                          shinyWidgets::actionBttn(inputId = "example_rnaseq_norm",
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
                           (e.g., outliers) and perform gene filtering and normalization."),
                          
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
                          
                          shinyWidgets::awesomeCheckbox(inputId = "outlier_rnaseq_norm",
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
                              prompter::add_prompt(message = "The experimental groups are the 
                                       groups that you would like to compare, like disease status groups.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          uiOutput("UI_groupselect_rnaseq_norm"),
                          
                          htmlOutput("experimentallevels_rnaseq_norm"),
                          br(),
                          
                          # Transformation
                          h4(strong(tags$span(
                            "3. Transformation",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Select the transformation method. 
                                       The aim of data transformation is to make the 
                                       data more normally distributed.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          uiOutput("UI_transformation_rnaseq_norm"),
                          br(),
                          
                          
                          # Filtering
                          h4(strong(tags$span(
                            "4. Filtering",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Select the minimum count required 
                                         for at least s samples, where s is the smallest 
                                         experimental group size", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          uiOutput("UI_filtering_rnaseq_norm"),
                          br(),
                          
                          
                          #Normalization
                          h4(strong(tags$span(
                            "4. Normalization",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Select the normalization method 
                                       and whether you would like to do the normalization
                                       on all samples or per experimental group.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          selectInput(
                            inputId = "normMeth_rnaseq_norm", 
                            label = NULL, 
                            choices = c("Quantile normalization" = "Quantile",
                                        "Continue without normalization" = "None"),
                            selected = "None"),
                          
                          conditionalPanel(
                            condition = "input.normMeth_rnaseq_norm == `Quantile`",
                            shinyWidgets::prettyRadioButtons(
                              inputId = "perGroup_rnaseq_norm", 
                              label = NULL, 
                              choices = c("Use all samples",
                                          "Per experimental group"),
                              inline = TRUE, 
                              status = "danger",
                              fill = TRUE)
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
                              prompter::add_prompt(message = "Click to perform the pre-processing!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          
                          shinyWidgets::actionBttn(inputId = "start_preprocessing_rnaseq_norm",
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
                              prompter::add_prompt(message = "Select which experimental groups you want 
                                       to compare to each other. Comparisons are commonly defined 
                                                   as `Disease - Control` or `Treatment - Control`.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_rnaseq_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Add covariates",
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
                              prompter::add_prompt(message = "Add gene annotations to 
                                       the output.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::materialSwitch(inputId = "addAnnotation_rnaseq_norm",
                                                       label = "Add gene annotations", 
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
                              prompter::add_prompt(message = "Click to start the statistical analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_statistics_rnaseq_norm",
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
               # Gene set analysis
               #*****************************************************************#
               tabPanel("Gene set analysis", value = "panel_ORA_rnaseq_norm", 
                        icon = icon("fas fa-list"),
                        
                        sidebarPanel(
                          h2(strong("Gene set analysis")),
                          
                          shinyWidgets::prettyRadioButtons(
                            inputId = "ORA_or_GSEA_rnaseq_norm",
                            label = NULL, 
                            choices = c("ORA", "GSEA"),
                            selected = "ORA",
                            status = "success",
                            fill = TRUE,
                            inline = TRUE),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm == 'ORA'",
                            h5("With",strong("Overrepresentation Analysis (ORA),"),"dysregulated 
                          processes and pathways can be idenified. These processes/pathways 
                          are identified by testing whether their genes are overrepresented among the (most) significant genes."),
                          ),
                          
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm == 'GSEA'",
                            h5("With", strong("Gene Set Enrichment Analysis (GSEA),"), "you can find dysregulated 
                          processes and pathways. These processes/pathways are identified by testing whether 
                               their genes show concordant changes in the data."),
                          ),
                          
                          hr(),
                          
                          h4(strong(tags$span(
                            "1. Select comparison",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Choose for which comparison you want to perform gene set analysis.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          uiOutput("UI_comparisons_view_ORA_rnaseq_norm"),
                          br(),
                          
                          h4(strong(tags$span(
                            "2. Select gene set collection",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "A gene set collection contains different sets of genes that are linked to 
                                       biological processes (GO-BP), molecular functions (GO-MF),
                                       cellular components (GO-CC), or biological pathways (WikiPathways and KEGG).", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          selectInput(inputId = "geneset_ORA_rnaseq_norm",
                                      label = NULL,
                                      choices = c("GO-BP",
                                                  "GO-MF",
                                                  "GO-CC",
                                                  "WikiPathways",
                                                  "KEGG")),
                          
                          br(),
                          
                          # ORA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm=='ORA'",
                            h4(strong(tags$span(
                              "3. Select genes",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select which genes are used in the analysis.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Up/down regulated genes only
                            shinyWidgets::prettyRadioButtons(
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
                            shinyWidgets::prettyRadioButtons(
                              inputId = "topNorThres_rnaseq_norm",
                              label = "Select genes based on ...", 
                              choices = c("Threshold", 
                                          "Top N"),
                              status = "danger",
                              selected = "Top N",
                              fill = TRUE,
                              inline = TRUE),
                            
                            # Select genes based on logFC/P value threshold
                            conditionalPanel(
                              condition = "input.topNorThres_rnaseq_norm=='Threshold'",
                              
                              #P value threshold
                              numericInput(
                                inputId = "p_thres_ORA_rnaseq_norm",
                                label = "P threshold",
                                value = 0.05),
                              
                              # Raw or adjusted P value?
                              shinyWidgets::prettyRadioButtons(
                                inputId = "rawp_ORA_rnaseq_norm", 
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
                            )
                            
                          ),
                          
                          # GSEA options
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm=='GSEA'",
                            h4(strong(tags$span(
                              "3. Select ranking variable",
                              tags$span(
                                icon(
                                  name = "question-circle",
                                ) 
                              ) |>
                                prompter::add_prompt(message = "Select the variable on which the GSEA should be based.", 
                                                     position = "right",
                                                     size = "large")
                            ))),
                            # Ranking variable
                            shinyWidgets::prettyRadioButtons(
                              inputId = "ranking_GSEA_rnaseq_norm",
                              label = NULL, 
                              choices = 
                                c(" logFC" = "logFC",
                                  "-log p-value" = "pvalue",
                                  "-log p-value x sign logFC" = "signed_pvalue"),
                              selected = "logFC",
                              status = "danger",
                              fill = TRUE),
                            
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
                              prompter::add_prompt(message = "To link genes to gene sets, 
                              we need to know which gene identifiers are used.", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          
                          # Select organism
                          selectInput(inputId = "organism_ORA_rnaseq_norm",
                                      label = tags$span(
                                        "Organism", 
                                        tags$span(
                                          icon(
                                            name = "question-circle",
                                          ) 
                                        ) |>
                                          prompter::add_prompt(message = "Select the organism. 
                                               This information is needed to match the gene IDs.", 
                                                               position = "right",
                                                               size = "large")
                                      ),
                                      choices = c("Bos taurus",
                                                  "Caenorhabditis elegans",
                                                  "Homo sapiens",
                                                  "Mus musculus", 
                                                  "Rattus norvegicus"),
                                      selected = "Homo sapiens"),
                          
                          uiOutput("UI_geneID_ORA_rnaseq_norm"),
                      
                          br(),
                         
                           # Calculate!
                          h4(strong(tags$span(
                            "5. Perform analysis",
                            tags$span(
                              icon(
                                name = "question-circle",
                              ) 
                            ) |>
                              prompter::add_prompt(message = "Click to start the gene set analysis!", 
                                                   position = "right",
                                                   size = "large")
                          ))),
                          shinyWidgets::actionBttn(inputId = "calculate_ORA_rnaseq_norm",
                                                   label = "Calculate",
                                                   style = "simple",
                                                   color = "warning",
                                                   icon = icon("sync")),
                          br(),
                          hr(),
                          uiOutput("UI_ORAreport_rnaseq_norm"),
                          uiOutput("UI_GSEAreport_rnaseq_norm")
                          
                        ),
                        mainPanel(
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm=='ORA'",
                            uiOutput("UI_output_ORA_rnaseq_norm")
                          ),
                          conditionalPanel(
                            condition = "input.ORA_or_GSEA_rnaseq_norm=='GSEA'",
                            uiOutput("UI_output_GSEA_rnaseq_norm")
                          )
                          
                          
                        )
               ), # Tab panel
               ###################################################################
               
               #  Documentation
               
               ###################################################################
               tabPanel("Documentation", value = "documentation", 
                        icon = icon("info"),
                        
                        tags$iframe(src="docs.html",
                                    width="100%",
                                    style="height: 85vh;",
                                    scrolling="yes",
                                    frameborder="0")
               ) # Tab panel
               
               
               ###################################################################
               
               #  END
               
               ###################################################################
    ) # EO navbarPage
  ) # EO fluidPage
) # EO tagList

