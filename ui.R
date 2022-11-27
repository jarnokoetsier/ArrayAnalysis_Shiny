
ui <- tagList(
  
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(6) {
                           float: right;
                           }
                           .my_style_1{ 
                             background-image: url(Background1.jpg);
                           }
                           
                           .my_style_1 { margin-top: -20px; }
                           
                           .my_style_1 { width: 100%; }
                           
                           .container-fluid { padding-left: 0; padding-right: 0; }
                           
                           .my_style_1 { position: absolute; left: 0; }
                           
                           "))),
  
  
  fluidPage(

  useSweetAlert(),
  
  navbarPage(title = "ArrayAnalysis", id = "navbar",
             
             ###################################################################
             #  Data selection                                                   
             ###################################################################
             tabPanel("Data accession", 
                      value = "panel1", 
                      icon = icon("fas fa-home"), class = "my_style_1",
                      
                      
                      #********************************************************#
                      #   Data selection
                      #********************************************************#
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      br(),
                      
                      fluidRow(
                        
                        column(4, offset = 4, 
                               align = "center", 
                               style = "background-color:#FFFFFF;",
                               
                               br(),
                               
                               img(src = "logo.png", width = "100%"),
                               
                               h1(strong(span(style = "color:#000000", "Welcome to ArrayAnalysis!"))),
                               
                               h5(span(style = "color:#000000", "Get started by entering your database 
                                  accession number or uploading your 
                                  own files.")),
                               
                               br(),

                               #Choose database
                               radioGroupButtons(
                                 inputId = "database",
                                 label = NULL,
                                 choices = c("GEO", 
                                             "ArrayExpress", 
                                             "Upload"),
                                 status = "danger"
                               ),
                               
                               #Enter database accession  
                               uiOutput("getdatabaseout"),
                               
                               #Upload CELs
                               uiOutput("uploadcelsout"),
                               
                               #Upload text file
                               uiOutput("uploadtxtout"),
                               
                               #Normalized or raw data
                               prettyRadioButtons(
                                 inputId = "rawornormalizeddata",
                                 label = NULL, 
                                 choices = c("Raw", "Normalized"),
                                 inline = TRUE, 
                                 status = "danger",
                                 fill = TRUE),
                               
                               #Use example data
                               actionBttn(inputId = "example", 
                                          label = "Example",
                                          style = "jelly",
                                          color = "royal",
                                          icon("fas fa-mouse-pointer")),
                               
                               #Start the analysis
                               actionBttn(inputId = "startAnalysis",
                                          label = "Start",
                                          style = "jelly",
                                          color = "primary",
                                          icon = icon("arrow-right")),

                               #Get information
                               actionBttn(inputId = "infopanel1", 
                                          label = NULL,
                                          style = "simple",
                                          color = "success",
                                          icon = icon("info")),
                               
                               br(),
                               br(),
                               br()
                               
                               )
                        ),
                        
                      
                      #********************************************************#
                      #   Continue with saved data
                      #********************************************************#
                     
                       fluidRow(
                        column(4, offset = 4, align = "center",
                               br(),
                               
                               #Continue with saved data
                               actionBttn(inputId = "continue", 
                                          label = "Continue with saved data",
                                          style = "simple",
                                          color = "warning",
                                          icon = icon("fas fa-sign-in-alt"))
                        )
                        
                      ),
                      
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
                      
                      
                      #********************************************************#
                        
             ),
             
             
             
             
             
             ###################################################################
             #  Grouping
             ###################################################################
             
             tabPanel("Grouping", value = "panel2", 
                      icon = icon("fas fa-layer-group"),
                      
                      #********************************************************#
                      #   Option panel
                      #********************************************************#
                      
                      sidebarPanel(
                        
                        #Title
                        h2(strong("Sample grouping")),
                        
                        #Subtitle (explanation)
                        h5("In the grouping step, you can classify the samples 
                           into multiple experimental groups."),
                        
                        hr(),
                        
                        #Paired data?
                        awesomeCheckbox(
                          inputId = "paired",
                          label = "Paired design", 
                          value = FALSE,
                          status = "danger"
                        ),
                        
                        #Select how to make groups (dataset, description file,
                        #or manual grouping)
                        uiOutput("makegroupsui"),
                        
                        #Upload description file
                        uiOutput("uidescriptionfile"),
                        
                        #Explanation of description file
                        uiOutput("uitextdescription"),

                        #Select grouping variable (from dataset)
                        uiOutput("groups"),
                        
                        #Select pairing variable (from dataset)
                        uiOutput("pairs"),
                        
                        #Select grouping variable (from description file)
                        uiOutput("groups1"),
                        
                        #Select pairing variable (from description file)
                        uiOutput("pairs1"),
                        
                        #Multi-input for manual grouping
                        uiOutput("owngroupsout"),
                        
                        hr(),
                        
                        #Go to previous tab
                        actionBttn(inputId = "meta.back", 
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        #Go to next tab
                        actionBttn(inputId = "meta.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
                        
                        
                      ),
                      
                      #********************************************************#
                      #   Output panel
                      #********************************************************#
                      
                      mainPanel(
                        
                        tabsetPanel(
                          
                          #Meta table
                          tabPanel("Meta table", 
                                   icon = icon("fas fa-file"),
                                   
                                   #Table
                                   dataTableOutput(outputId = "grouping") %>% 
                                     withSpinner(color="#0dc5c1"),
                                   
                                   #Download table
                                   downloadButton("downloadmeta", 
                                                  "Download meta table")
                                   
                                   ),
                          
                          #Meta heatmap
                          tabPanel("Meta heatmap", 
                                   icon = icon("fas fa-mouse-pointer"),
                                   
                                   #Heatmap (interactive)
                                   plotlyOutput("sampleheatmap", 
                                                width = "1200px", 
                                                height="700px") %>% 
                                     withSpinner(color="#0dc5c1")
                                   
                                   ),
                          
                          #Meta bar chart
                          tabPanel("Meta bar chart", 
                                   icon = icon("fas fa-file"),
                                   
                                   br(),
                                   
                                   fluidRow(
                                     column(3, 
                                            uiOutput("uimetaBar1"),
                                            
                                            uiOutput("uimetaBar3")),
                                     
                                     column(3,
                                            uiOutput("uimetaBar2"),
                                            
                                            uiOutput("uimetaBar4"))
                                     ),
                                   
                                   hr(),
                                   
                                   plotOutput("metaBar", 
                                                width = "1000px", 
                                                height="500px") %>% 
                                     withSpinner(color="#0dc5c1")
                                   
                          )
                          
                          
                        )
                        
                        
                      )
                      
                      
                      
             ),
             
             ###################################################################
             #  Pre-processing
             ###################################################################
             
             tabPanel("Pre-processing", value = "panel3", icon = icon("sync"),
                      
                      sidebarPanel(
                        
                        #******************************************************#
                        #   Information
                        #******************************************************#
                        
                        h2(strong("Pre-processing")),
                        
                        h5("In this pre-processing step, you can remove samples 
                           (e.g., outliers), perform normalization, 
                           and choose your desired probeset annotation."),
                        
                        h5(strong("NOTE:"), "After you have pre-processed the 
                        data, you can can always change the grouping variable 
                        in the", em("Grouping"), "tab
                        without needing to perform the pre-processing again."),
                        
                        hr(),
                        
                        #******************************************************#
                        #   Options
                        #******************************************************#
                        
                        #Remove outliers
                        h4(strong("1. Remove samples")),
                        
                        awesomeCheckbox(inputId = "outlier",
                                        label = "Keep all samples",
                                        value = TRUE,
                                        status = "danger"),
                        
                        uiOutput("outliersout"),
                        
                        br(),
                        
                        #Normalization
                        h4(strong("2. Normalization")),
                        
                        radioGroupButtons(inputId = "normMeth", 
                                    label = NULL, 
                                    choices = c("RMA","GCRMA","PLIER","None"),
                                    status = "danger"),
                        
                        prettyRadioButtons(inputId = "perGroup", 
                                    label = NULL, 
                                    choices = c("Use all arrays",
                                                "Per experimental group"),
                                    inline = TRUE, 
                                    status = "danger",
                                    fill = TRUE),
                        
                        br(),
                        
                        #Annotation
                        h4(strong("3. Annotation")),
                        
                        selectInput(inputId = "annotations",
                                    label = NULL,
                                    choices = c("No annotations",
                                                "Custom annotations",
                                                "Upload annotation file"),
                                    selected = "Custom annotations"),
                        
                        conditionalPanel(
                          condition = "input.annotations=='Custom annotations'",
                          
                          selectInput(inputId = "species",
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
                                      
                                      selected = "Homo sapiens"),
                          
                          selectInput(inputId = "CDFtype",
                                      label = "Annotation format",
                                      choices = c("ENTREZG","REFSEQ","ENSG",
                                                  "ENSE","ENST","VEGAG","VEGAE",
                                                  "VEGAT","TAIRG","TAIRT","UG",
                                                  "MIRBASEF","MIRBASEG"))
                        ),
                        conditionalPanel(
                          condition = "input.annotations=='Upload annotation file'",
                          
                          fileInput(inputId = "annot_file",
                                    label = "Upload annotation file",
                                    multiple = FALSE)
                        ),
                        
                        br(),
                        
                        #Pre-processing
                        h4(strong("4. Pre-processing")),
                        
                        actionBttn(inputId = "preprocessing",
                                   label = "Calculate",
                                   style = "simple",
                                   color = "warning",
                                   icon = icon("sync")),
                        
                        br(),
                        
                        br(),
                        
                        br(),
                        
                        #Save data
                        h4(strong("5. Save for later use")),
                        
                        actionBttn(inputId = "save",
                                   label = "Save",
                                   style = "simple",
                                   color = "success",
                                   icon = icon("fas fa-save")),
                        
                        br(),
                        
                        br(),
                        
                        hr(),
                        
                        #Go back
                        actionBttn(inputId = "ann.back",
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        #Go forward
                        actionBttn(inputId = "ann.proceed",
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
   
                      ),
                      
                      #********************************************************#
                      #   Output panel
                      #********************************************************#
                      mainPanel(
                        
                        tabsetPanel(
                          
                          #Boxplot normalized data
                          tabPanel("Boxplot (normalized)", 
                                   icon = icon("fas fa-file"),
                                   
                                   plotOutput("boxplotNorm", 
                                              height = "100%",
                                              width = "100%")  %>% 
                                     withSpinner(color="#0dc5c1", 
                                                 proxy.height = "400px")
                                   
                                   ),
                          
                          #Boxplot raw data
                          tabPanel("Boxplot (raw)", icon = icon("fas fa-file"),
                                   
                                   textOutput("alreadynormalized"),
                                   
                                   plotOutput("boxplotRaw",
                                              height = "100%",
                                              width = "100%") %>% 
                                     withSpinner(color="#0dc5c1", 
                                                 proxy.height = "400px")
                                   
                          ),
                          
                          #Density plot normalized data
                          tabPanel("Density plot", 
                                   icon = icon("fas fa-mouse-pointer"),
                                   
                                   br(),
                                   
                                   plotlyOutput("normhist") %>% 
                                     withSpinner(color="#0dc5c1")
                                   ),
                          
                          
                          #Correlation plot/heatmap    
                          tabPanel("Correlation plot", 
                                   icon = icon("fas fa-mouse-pointer"),
                                   
                                   br(),
                                   
                                   fluidRow(
                                     
                                     #Distance method
                                     column(4,
                                            pickerInput(
                                              inputId = "clusteroption1",
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
                                              inputId = "clusteroption2",
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
                                     
                                     selectInput(inputId = 'heatmaptheme',
                                                 label = NULL,
                                                 choices = c("Default", 
                                                             "Yellow-red", 
                                                             "Dark")),
                                    
                                      circle = TRUE, status = "info",
                                     icon = icon("fas fa-cog"), width = "300px",
                                     
                                     tooltip = tooltipOptions(
                                       title = "Click to change theme!")
                                     
                                   ),
                                   
                                   plotlyOutput("correlNormInt", 
                                                width = "1200px", 
                                                height="1000px") %>% 
                                     withSpinner(color="#0dc5c1", 
                                                 proxy.height = "400px")

                                   ),
                          
                          #Expression matrix
                          tabPanel("Expression matrix", 
                                   icon = icon("fas fa-mouse-pointer"),
                                   
                                   br(),
                                   
                                   dataTableOutput(outputId = "ExpressionMatrix"),
                                   
                                   downloadButton("downloadexpr", 
                                                  "Download expression matrix"),
                                   
                                   plotlyOutput("ExprBoxplot")%>% 
                                     withSpinner(color="#0dc5c1")
                          )
                          
                        ) #Tab set panel
                      
                        
                      ) #Main panel
                      
             ), #Tab panel
             
             
             
             ###################################################################
             #PCA
             ###################################################################
             
             tabPanel("PCA", value = "panel4", icon = icon("fas fa-cube"),
                      
                      sidebarPanel(
                        
                        #******************************************************#
                        #   Title and description
                        #******************************************************#
                        
                        h2(strong("PCA")),
                        
                        h5("Principal component analysis (PCA) shows the 
                        similarity between the samples and is particularly 
                           useful for the detection of outliers."),
                        
                        hr(),
                        
                        #******************************************************#
                        #   Option panel
                        #******************************************************#
                        
                        #3D plot
                        materialSwitch(
                          inputId = "plot3d",
                          label = "3D",
                          value = FALSE, 
                          status = "danger"),
                        
                        #X-axis
                        selectInput(inputId = "xpca", 
                                    label = "x-axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5",
                                                "PC6", "PC7", "PC8"),
                                    selected = "PC1"),
                        
                        #Y-axis
                        selectInput(inputId = "ypca", 
                                    label = "y-axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5", 
                                                "PC6", "PC7", "PC8"),
                                    selected = "PC2"),
                        
                        #Z-axis
                        uiOutput("uizpca"),
                        
                        hr(),
                        
                        #Go back
                        actionBttn(inputId = "pca.back",
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        #Go forward
                        actionBttn(inputId = "pca.ok",
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
                      ),
                      
                      #********************************************************#
                      #   Output panel
                      #********************************************************#
                      mainPanel(
                        
                        #2D PCA plot
                        uiOutput("uipca"),
                        
                        #3D PCA plot
                        uiOutput("uipca3d"),
                        
                        #Explained variances histogram
                        plotlyOutput("variances") %>% 
                          withSpinner(color="#0dc5c1")
                      )
                      
                      
             ),
             
             ###################################################################
             #Statistical analysis
             ###################################################################
             
             tabPanel("Statistical analysis", value = "panel5", 
                      icon = icon("fas fa-chart-bar"),
                      
                      navlistPanel(
                        
                        "Statistical analysis",
                        
                        #******************************************************#
                        #   Top table
                        #******************************************************#
                        tabPanel("Top table", value = "panel6",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "Top table"))),
                                 
                                 #Description
                                 h4("The top table is generated using the limma package and
                                 includes the statistics (e.g., P-value, logFC, etc.) 
                                 of the top-ranked genes."), 
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   
                                   #Select comparison
                                   column(4,
                                          uiOutput("uitoptablecomp")),
                                   
                                   
                                   #Select direction of comparison
                                   column(4,
                                          uiOutput("uitoptabledir"))
                                 ),
                                 
                                 hr(),
                                 
                                 #Top table
                                 DT::dataTableOutput("x1") %>% 
                                   withSpinner(color="#0dc5c1"),
                                 
                                 #Download top table
                                 downloadButton("downloadfinal", 
                                                "Download top table"),
                                 
                                 #Boxplot of selected gene
                                 plotlyOutput("topboxplot")
                                 
                                 ),
                        
                        
                        #******************************************************#
                        #   P-Value histogram
                        #******************************************************#
                        tabPanel("P-value histogram", value = "panel7",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "P-value histogram"))),
                                 
                                 #Description
                                 h4("The P-value histogram shows the distribution of the P-values.
                                    More information about the interpretation of a P-value histogram
                                    can be found", 
                                    a("here.", href="http://varianceexplained.org/statistics/interpreting-pvalue-histogram/")),
                                 
                                 hr(),
                                 
                                 #Select comparison
                                 uiOutput("uipcomp"),
                                 
                                 br(),
                                 
                                 #P-value histogram
                                 plotlyOutput("phist") %>% 
                                   withSpinner(color="#0dc5c1")
                                 
                                 ),
                        
                        
                        #******************************************************#
                        #   logFC histogram
                        #******************************************************#
                        tabPanel("logFC histogram", value = "panel8",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "logFC histogram"))),
                                 
                                 #Description
                                 h4("The logFC histogram shows the distribution of the
                                    log2 fold changes in gene expression levels."),
                                 
                                 hr(),
                                 
                                 #Select comparison
                                 uiOutput("uilogfccomp"),
                                 
                                 br(),
                                 
                                 #logFC histogram
                                 plotlyOutput("logfchist") %>% 
                                   withSpinner(color="#0dc5c1")
                                 
                                 ),
                        
                        #******************************************************#
                        #   Pie chart
                        #******************************************************#
                        tabPanel("Pie chart", value = "panel_piechart",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "Pie chart"))),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   
                                   column(4,
                                          
                                          sliderInput(
                                            inputId = "pthreshold_pie",
                                            label = "P threshold",
                                            value = 0.05,
                                            min = 0,
                                            max = 0.5,
                                            step = 0.01
                                          ),
                                          
                                          sliderInput(
                                            inputId = "logFCthreshold_pie",
                                            label = "logFC threshold",
                                            value = 0,
                                            min = 0,
                                            max = 5,
                                            step = 0.1
                                          ),

                                   ),
                                   
                                   column(4,
                                          
                                          #Raw or adjusted p-value?
                                          prettyRadioButtons(
                                            inputId = "pchoice_pie",
                                            label = "P-value", 
                                            choices = 
                                              c("Raw P-value" = "raw", 
                                                "Adjusted P-value" = "adj")),
                                          
                                          #Select comparison
                                          uiOutput("uipiecomp")
                                   )
                                   
                                 ),
                                 
                                 br(),
                                 
                                 #Pie chart
                                 plotlyOutput("piechart") %>% 
                                   withSpinner(color="#0dc5c1")
                                 
                        ),
                        
                        
                        
                        #******************************************************#
                        #   Volcano plot
                        #******************************************************#
                        tabPanel("Volcano plot", value = "panel9",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "Volcano plot"))),
                                 
                                 #Description
                                 h4("A volcano plot shows the statistical signficance against
                                    the logFC and can thus be used to quickly identify
                                    genes of interest."), 
                                 
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   column(3,
                                          
                                          #Raw or adjusted p-value?
                                          prettyRadioButtons(
                                            inputId = "raworadj",
                                            label = "P-value", 
                                            choices = 
                                              c("Raw P-value" = "raw", 
                                                "Adjusted P-value" = "adj")),
                                          
                                          #Plot
                                          actionBttn(inputId = "volcano.ok", 
                                                     label = "Plot",
                                                     style = "jelly",
                                                     color = "primary")
                                          
                                   ),
                                   
                                   column(4,
                                          
                                          #P-value threshold
                                          numericInput(
                                            inputId = "pthreshold",
                                            label = "P threshold",
                                            value = 0.05),
                                          
                                          #logFC threshold
                                          numericInput(
                                            inputId = "logfcthreshold",
                                            label = "logFC threshold",
                                            value = 1)
                                   ),
                                   
                                   column(4,
                                          
                                          #Select comparison
                                          uiOutput("uivolcanocomp"),
                                          
                                          #Select direction of comparison
                                          uiOutput("uivolcanodir")
                                
                                   )
                                   
                                 ),
                                 
                                 hr(),
                                 
                                 #Volcano plot
                                 plotlyOutput("volcano") %>% withSpinner(color="#0dc5c1"),
                                 
                                 #Top table of significant genes
                                 dataTableOutput("voltable") %>% withSpinner(color="#0dc5c1"),
                                 
                                 #Download the top table of significant genes
                                 downloadButton("downloadvol", "Download table")
                                 ),
                        
                        
                        #******************************************************#
                        #   Venn diagram
                        #******************************************************#
                        tabPanel("Venn diagram", value = "panel10",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", 
                                                "Venn Diagram"))),
                                 
                                 #Description
                                 h4("A venn diagram shows the overlap of significant genes
                                    for different comparisons."),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   
                                   column(3,
                                      
                                          #Raw or adjusted p-value
                                          prettyRadioButtons(
                                            inputId = "raworadjvenn",
                                            label = "P-value", 
                                            choices = 
                                              c("Raw P-value" = "raw", 
                                                "Adjusted P-value" = "adj"))
                                          
                                   ),
                                   
                                   
                                   column(4,
                                          
                                          #P-value threshold
                                          numericInput(
                                            inputId = "pthresholdvenn",
                                            label = "P threshold",
                                            value = 0.05)
                                          
                                   ),
                                   
                                   column(4,
                                          
                                          #logFC theshold
                                          numericInput(
                                            inputId = "logfcthresholdvenn",
                                            label = "logFC threshold",
                                            value = 0)
                                          
                                   )
                                   
                  
                                   
                                 ),
                                 
                                 hr(),
                                 
                                 #Not enough comparisons selected
                                 uiOutput("uicontrast"), 
                                 
                                 #Legend
                                 htmlOutput("legendvenn"),
                                 
                                 br(),
                                 
                                 #Venn diagram
                                 plotOutput("venndiagram")
                                 ),
                        
                        #******************************************************#
                        #   Heatmap
                        #******************************************************#
                        tabPanel("Heatmap", value = "panel12",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", "Heatmap"))),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   column(4,
                                          
                                          #Distance calulcation method
                                          pickerInput(
                                            inputId = "clusteroption3",
                                            label = "Distance calculation method",
                                            choices = c("Pearson","Spearman","Euclidean"),
                                            options = list(
                                              style = "btn-primary")),
                                          
                                          #Clustering method
                                          pickerInput(
                                            inputId = "clusteroption4",
                                            label = "Clustering method",
                                            choices = c("ward.D2","single",
                                                        "complete","average",
                                                        "mcquitty","median","centroid"),
                                            options = list(
                                              style = "btn-info"))
                                   ),
                                   
                                   column(4,
                                          #Select comparison
                                          uiOutput("uiheatmapcomp"),
                                          
                                          #Number of features
                                          sliderInput(
                                            inputId = "topfeatures",
                                            label = "Number of features",
                                            value = 25,
                                            min = 5,
                                            max = 50,
                                            step = 1
                                          ))
                                   
                                 ),
                                 
                                 hr(),
                                 
                                 #Theme
                                 dropdownButton(
                                   tags$h3("Theme"),
                                   
                                   selectInput(inputId = 'heatmaptheme1',
                                               label = NULL,
                                               choices = c("Default", 
                                                           "Yellow-red", 
                                                           "Dark")),
                                   
                                   circle = TRUE, status = "info",
                                   icon = icon("fas fa-cog"), width = "300px",
                                   
                                   tooltip = tooltipOptions(
                                     title = "Click to change theme!")
                                   
                                 ),
                                 
                                 #Heatmap
                                 plotlyOutput("topfeatureheatmap",
                                              height = "1000px",
                                              width = "1200px") %>% withSpinner(color="#0dc5c1")
                                 
                        ),
                        
                        #******************************************************#
                        #   GO analysis
                        #******************************************************#
                        tabPanel("GO analysis", value = "panel11",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C",
                                                "GO analysis"))),
                
                                 #Method: ORA or GSEA
                                 radioGroupButtons(
                                   inputId = "oraGO",
                                   label = NULL,
                                   choices = 
                                     c("Over-representation analysis (ORA)" = "ora",
                                       "Gene set enrichment analysis (GSEA)" = "gsea"),
                                   checkIcon = list(
                                     yes = icon("ok", 
                                                lib = "glyphicon"),
                                     no = icon("remove",
                                               lib = "glyphicon"))
                                 ),
                                 
                                 br(),
                                 
                                 
                                 tabsetPanel(
                                   
                                   tabPanel("Gene ontology (GO)",
                                            
                                            br(),
                                              
                                              fluidRow(
                                                
                                                column(2,
                                                       
                                                       conditionalPanel(
                                                         condition = "input.oraGO=='ora'",
                                                         
                                                         #Adjusted or raw p-value (ORA)
                                                         prettyRadioButtons(
                                                           inputId = "raworadjgoa1",
                                                           label = "P-value", 
                                                           choices = c("Raw P-value" = "raw", 
                                                                       "Adjusted P-value" = "adj"))
                                                       ),
                                                       
                                                       conditionalPanel(
                                                         condition = "input.oraGO=='gsea'",
                                                         
                                                         #Ranking statistic (GSEA)
                                                         prettyRadioButtons(
                                                           inputId = "rankingGO",
                                                           label = "Ranking statistic", 
                                                           choices = c("logFC", 
                                                                       "-log(P-value)"))
                                                       ),
                                                    
                                                     br(),
                                                     
                                                     #Calculate
                                                     actionBttn(inputId = "goa.ok1", 
                                                                label = "Calculate",
                                                                style = "jelly",
                                                                color = "primary")
                                                     
                                              ),
                                              conditionalPanel(
                                                condition = "input.oraGO=='ora'",
                                              
                                              column(3,
                                                       
                                                     #P-value threshold (ORA)
                                                       numericInput(
                                                         inputId = "pthresholdgoa1",
                                                         label = "P threshold",
                                                         value = 0.05),
                                                       
                                                     #logFC threshold (ORA)
                                                       numericInput(
                                                         inputId = "logFCthresholdgoa1",
                                                         label = "logFC threshold",
                                                         value = 0)
                                                     )
                                                     
                                              ),
                                              
                                              column(3,
                                                     
                                                     #Select comparison
                                                     uiOutput("uigoacomp1"),
                                                     
                                                     #Onotology
                                                     pickerInput(
                                                       inputId = "ontology1",
                                                       label = "Ontology", 
                                                       choices = c("Biological process" = "BP", 
                                                                   "Molecular function" = "MF",
                                                                   "Cellular component" = "CC"),
                                                       selected = "Biological process",
                                                       options = list(
                                                         style = "btn-info"))
      
                                              ),
                                              
                                              column(3,
                                                     
                                                     #Select organism
                                                     uiOutput("uigoaOrg"),
                                                     
                                                     #select chiptype
                                                     uiOutput("uigoaChip")
                                                     
                                              )
                                              
                                            ),
                                            
                                            hr(),
                                            
                                            #Error
                                            textOutput("errorgoa1"),
                                            
                                            #Dot plot
                                            plotOutput("GOplot1") %>% 
                                              withSpinner(color="#0dc5c1"),
                                            
                                            #Output table
                                            dataTableOutput("goa1"),
                                            
                                            #GSEA plot
                                            plotOutput("GOplotgsea1")
                                   )
                                 )
                                 
                               
                                 
                        ),
                        
                        #******************************************************#
                        #   Pathway analysis
                        #******************************************************#
                        tabPanel("Pathway analysis", value = "panel11a",
                                 
                                 #Title
                                 h1(strong(span(style = "color:#3A3B3C", "Pathway analysis"))),
                                 
                                 #Method: ORA or GSEA
                                 radioGroupButtons(
                                   inputId = "oraPathway",
                                   label = NULL,
                                   choices = 
                                     c("Over-representation analysis (ORA)" = "ora",
                                       "Gene set enrichment analysis (GSEA)" = "gsea"),
                                   checkIcon = list(
                                     yes = icon("ok", 
                                                lib = "glyphicon"),
                                     no = icon("remove",
                                               lib = "glyphicon"))
                                 ),
                                 
                                 br(),
                                 
                                 tabsetPanel(
                                   
                                   #...........................................#
                                   #   Wikipathways
                                   #...........................................#
                                   tabPanel("WikiPathways",
                                            
                                            br(),
                                            
                                            fluidRow(
                                              
                                              column(2,
                                                     
                                                     conditionalPanel(
                                                       condition = "input.oraPathway=='ora'",
                                                       
                                                       #Raw or adjusted P-value (ORA)
                                                       prettyRadioButtons(
                                                         inputId = "raworadjwp",
                                                         label = "P-value", 
                                                         choices = c("Raw P-value" = "raw", 
                                                                     "Adjusted P-value" = "adj"))
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.oraPathway=='gsea'",
                                                       
                                                       #Ranking statistic (GSEA)
                                                       prettyRadioButtons(
                                                         inputId = "rankingWP",
                                                         label = "Ranking statistic", 
                                                         choices = c("logFC", 
                                                                     "-log(P-value)"))
                                                     ),
                                                     
                                                     br(),
                                                     
                                                     #Calculate
                                                     actionBttn(inputId = "wp.ok", 
                                                                label = "Calculate",
                                                                style = "jelly",
                                                                color = "primary")
                                                     
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.oraPathway=='ora'",
                                                
                                                column(3,
                                                       
                                                       #P-value threshold (ORA)
                                                       numericInput(
                                                         inputId = "pthresholdwp",
                                                         label = "P threshold",
                                                         value = 0.05),
                                                       
                                                       #logFC threshold (GSEA)
                                                       numericInput(
                                                         inputId = "logFCthresholdwp",
                                                         label = "logFC threshold",
                                                         value = 0)
                                                )
                                                
                                              ),
                                              
                                              column(3,
                                                     
                                                     #Select comparison
                                                     uiOutput("uiwpcomp")
                                                     
                                              ),
                                              
                                              column(3,
                                                     
                                                     #Select organism
                                                     uiOutput("uiwpOrg"),
                                                     
                                                     #select chiptype
                                                     uiOutput("uiwpChip")
                                                     
                                              ),
                                              
                                            ),
                                            
                                            hr(),
                                            
                                            #Error
                                            textOutput("errorwp"),
                                            
                                            #Barplot
                                            plotlyOutput("WPplot") %>% 
                                              withSpinner(color="#0dc5c1"),
                                            
                                            #Output table
                                            dataTableOutput("wp" ),
                                            
                                            #GSEA plot
                                            plotOutput("WPplotgsea")
                                            ),
                                   
                                   
                                   #...........................................#
                                   #   KEGG
                                   #...........................................#
                                   
                                   tabPanel("KEGG",
                                            
                                            br(),
                                            
                                            fluidRow(
                                              
                                              column(2,
                                                     
                                                     conditionalPanel(
                                                       condition = "input.oraPathway=='ora'",
                                                       
                                                       #Raw or adjusted p-value (ORA)
                                                       prettyRadioButtons(
                                                         inputId = "raworadjkegg",
                                                         label = "P-value", 
                                                         choices = c("Raw P-value" = "raw", 
                                                                     "Adjusted P-value" = "adj"))
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.oraPathway=='gsea'",
                                                       
                                                       #Ranking statistic (GSEA)
                                                       prettyRadioButtons(
                                                         inputId = "rankingKEGG",
                                                         label = "Ranking statistic", 
                                                         choices = c("logFC", 
                                                                     "-log(P-value)"))
                                                     ),
                                                     
                                                     br(),
                                                     
                                                     #Calculate
                                                     actionBttn(inputId = "kegg.ok", 
                                                                label = "Calculate",
                                                                style = "jelly",
                                                                color = "primary")
                                                     
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.oraPathway=='ora'",
                                                
                                                column(3,
                                                       
                                                       #P-value threshold (ORA)
                                                       numericInput(
                                                         inputId = "pthresholdkegg",
                                                         label = "P threshold",
                                                         value = 0.05),
                                                       
                                                       #logFC threshold (ORA)
                                                       numericInput(
                                                         inputId = "logFCthresholdkegg",
                                                         label = "logFC threshold",
                                                         value = 0)
                                                       
                                                )
                                                
                                              ),
                                              
                                              column(3,
                                                     
                                                     #Select comparison
                                                     uiOutput("uikeggcomp")
                                                     
                                              ),
                                              
                                              
                                              column(3,
                                                     
                                                     #Select organism
                                                     uiOutput("uikeggOrg"),
                                                     
                                                     #select chiptype
                                                     uiOutput("uikeggChip")
                                                     
                                              )
                                              
                                            ),
                                            
                                            hr(),
                                            
                                            #Error
                                            textOutput("errorkegg"),
                                            
                                            #Barplot
                                            plotlyOutput("KEGGplot") %>% 
                                              withSpinner(color="#0dc5c1"),
                                            
                                            #Ouput table
                                            dataTableOutput("kegg" ),
                                            
                                            #GSEA plot
                                            plotOutput("KEGGplotgsea")
                                            )
                                 )
                                 
                        ),
                        
                        
                        
                        #******************************************************#
                        #   Data exploration
                        #******************************************************#
                        tabPanel("Data exploration", value = "panel12a",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "Data exploration"))),
                                 
                                 h4("Explore your data with new interactive visualizations!"),
                                 
                                 br(),
                                 
                                 tabsetPanel(
                                   
                                   tabPanel("Parallel Coordinates Plot", icon = icon("fas fa-mouse-pointer"),
                                            
                                            br(),
                                            
                                            fluidRow(
                                              
                                              column(4, uiOutput("uiparacomp"))
                                              
                                            ),

                                            plotlyOutput("allpara") %>% 
                                              withSpinner(color="#0dc5c1"),
                                            
                                            plotlyOutput("topfeaturepara")
                                   )     
                              )
                        )
                      )
                        
                      
                      
             ),
             
             ###################################################################
             #  Documentation
             ###################################################################
             
             tabPanel("Documentation", value = "panel13", icon = icon("far fa-question-circle")
                      
                      
                      
                      
             )
             
             
             
             
  )
)
)

