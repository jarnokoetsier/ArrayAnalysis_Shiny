
###############################################################################################################################

#USER INTERFACE

###############################################################################################################################

ui <- tagList(
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(6) {
                           float: right;
                           }
                           "))),
  
  
  fluidPage(

  useSweetAlert(),
  
  navbarPage("ArrayAnalysis", id = "navbar",
             ################################################################################################################################
             #Enter database accession number
             ################################################################################################################################
             
             tabPanel("Data accession", value = "panel1", icon = icon("fas fa-home"),
                      
                      fluidRow(
                        column(12, align = "center",
                               img(src = "logo.png", height = "250"))
                      ),
                      
                      hr(),
                      
                      br(),
                      
                      fluidRow(
                        column(12, align = "center",
                               
                               h1(strong("Welcome to ArrayAnalysis!")),
                               
                               h5("Get started by entering your database accession number or uploading your CEL files.")
                        )
                      ),
                      
                      br(),
                      
                      fluidRow(
                        column(4, offset = 4, align = "center", style = "background-color:#E5E4E2;",
                               
                               br(),

                                 radioGroupButtons(
                                   inputId = "database",
                                   label = NULL,
                                   choices = c("GEO", "ArrayExpress", "Upload CELs"),
                                   status = "danger"
                                 ),
                                 
                                 uiOutput("getdatabaseout"),
                                 
                                 uiOutput("uploadcelsout"),
                                 
                                 actionBttn(inputId = "example", 
                                            label = "Example",
                                            style = "jelly",
                                            color = "royal",
                                            icon = icon("play-circle")),
                                 
                                 actionBttn(inputId = "downloaddata", 
                                            label = "Next",
                                            style = "jelly",
                                            color = "primary",
                                            icon = icon("arrow-right")),
                               
                               actionBttn(inputId = "infopanel1", 
                                          label = NULL,
                                          style = "simple",
                                          color = "success",
                                          icon = icon("info")),
                               
                               br(),
                               
                                 
                                br() 
                               ),
                        
                        fluidRow(
                          column(4, offset = 4, align = "center",
                                 br(),
                                 
                                 actionBttn(inputId = "continue", 
                                            label = "Continue with saved data",
                                            style = "simple",
                                            color = "warning",
                                            icon = icon("fas fa-sign-in-alt"))
                                 )
                         
                        )
                       
                        
                      ),
                      
                      
                    
                      
             ),
             
             
             
             
             
             ################################################################################################################################
             #Grouping
             ################################################################################################################################
             
             tabPanel("Grouping", value = "panel2", icon = icon("fad fa-layer-group"),
                      
                      
                      
                      sidebarPanel(
                        
                        h2(strong("Sample grouping")),
                        
                        h5("In the grouping step, you can classify the samples into multiple experimental groups."),
                        
                        hr(),
                        
                        uiOutput("makegroupsui"),
                        
                        uiOutput("uidescriptionfile"),
                        
                        uiOutput("uitextdescription"),
                        
                        uiOutput("uidownloadtemplate"),

                        uiOutput("groups"),
                        
                        uiOutput("owngroupsout"),
                        
                        hr(),
                        
                        actionBttn(inputId = "meta.back", 
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        actionBttn(inputId = "meta.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
                        
                        
                      ),
                      
                      mainPanel(
                        
                        dataTableOutput(outputId = "grouping") %>% withSpinner(color="#0dc5c1"),
                        
                        downloadButton("downloadmeta", "Download meta table")
                      )
                      
                      
                      
             ),
             
             ################################################################################################################################
             #Pre-processing
             ################################################################################################################################
             
             tabPanel("Pre-processing", value = "panel3", icon = icon("refresh"),
                      
                      sidebarPanel(
                        
                        h2(strong("Pre-processing")),
                        
                        h5("In this pre-processing step, you can remove samples (e.g. outliers), perform normalization, and choose your desired probeset annotation."),
                        
                        h5(strong("NOTE:"), "After you have pre-processed the data, you can can always change the grouping variable in the", em("Grouping"), "tab
                        without needing to perform the pre-processing again."),
                        
                        hr(),
                        
                        h4(strong("1. Remove samples")),
                        
                        awesomeCheckbox(inputId = "outlier",
                                        label = "Keep all samples",
                                        value = TRUE,
                                        status = "danger"),
                        
                        uiOutput("outliersout"),
                        
                        br(),
                        
                        h4(strong("2. Normalization")),
                        
                        radioGroupButtons(inputId = "normMeth", 
                                    label = NULL, 
                                    choices = c("RMA","GCRMA","PLIER","None"),
                                    status = "danger"),
                        
                        prettyRadioButtons(inputId = "perGroup", 
                                    label = NULL, 
                                    choices = c("Use all arrays","Per experimental group"),
                                    inline = TRUE, 
                                    status = "danger",
                                    fill = TRUE),
                        
                        br(),
                        
                        h4(strong("3. Annotation")),
                        
                        selectInput(inputId = "annotations",
                                    label = NULL,
                                    choices = c("No annotations","Custom annotations","Upload annotation file")),
                        
                        conditionalPanel(
                          condition = "input.annotations=='Custom annotations'",
                          
                          selectInput(inputId = "species",
                                      label = "Species",
                                      choices = c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus","Caenorhabditis elegans","Canis familiaris", "Danio rerio","Drosophila melanogaster","Gallus gallus","Homo sapiens","Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus","Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa"), 
                                      selected = "Homo sapiens"),
                          
                          selectInput(inputId = "CDFtype",
                                      label = "Annotation format",
                                      choices = c("ENTREZG","REFSEQ","ENSG","ENSE","ENST","VEGAG","VEGAE","VEGAT","TAIRG","TAIRT","UG","MIRBASEF","MIRBASEG"))
                        ),
                        conditionalPanel(
                          condition = "input.annotations=='Upload annotation file'",
                          
                          fileInput(inputId = "annot_file",
                                    label = "Upload annotation file",
                                    multiple = FALSE)
                        ),
                        
                        br(),
                        
                        h4(strong("4. Pre-processing")),
                        
                        actionBttn(inputId = "preprocessing",
                                   label = "Calculate",
                                   style = "simple",
                                   color = "warning",
                                   icon = icon("refresh")),
                        
                        br(),
                        
                        br(),
                        
                        br(),
                        
                        h4(strong("5. Save for later use")),
                        
                        actionBttn(inputId = "save",
                                   label = "Save",
                                   style = "simple",
                                   color = "success",
                                   icon = icon("fas fa-save")),
                        
                        br(),
                        
                        br(),
                        
                        hr(),
                        
                        actionBttn(inputId = "ann.back",
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        actionBttn(inputId = "ann.proceed",
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
   
                      ),
                      
                      mainPanel(
                        
                        tabsetPanel(
                          
                          tabPanel("Boxplot (normalized)", icon = icon("fas fa-file"),
                                   
                                   
                                   plotOutput("boxplotNorm", 
                                              height = "100%",
                                              width = "100%")  %>% withSpinner(color="#0dc5c1", proxy.height = "400px")
                                   
                                   ),
                          
                          tabPanel("Boxplot (raw)", icon = icon("fas fa-file"),
                                   
                                   plotOutput("boxplotRaw",
                                              height = "100%",
                                              width = "100%") %>% withSpinner(color="#0dc5c1", proxy.height = "400px")
                                   
                          ),
                          
                          tabPanel("Density plot", icon = icon("fas fa-mouse-pointer"),
                                   
                                   br(),
                                   
                                   plotlyOutput("normhist") %>% withSpinner(color="#0dc5c1")
                                   
                                   ),
                          
                              
                          tabPanel("Correlation plot", icon = icon("fas fa-file"),
                                   
                                   br(),
                                   
                                   fluidRow(
                                     column(4,
                                            pickerInput(
                                              inputId = "clusteroption1",
                                              label = "Distance calculation method",
                                              choices = c("Pearson","Spearman","Euclidean"),
                                              options = list(
                                                style = "btn-primary"))
                                     ),
                                     
                                     column(4,
                                            pickerInput(
                                              inputId = "clusteroption2",
                                              label = "Clustering method",
                                              choices = c("ward.D2","single","complete","average","mcquitty","median","centroid"),
                                              options = list(
                                                style = "btn-info"))
                                     )),
                                   
                                   hr(),
                                     
                                   plotOutput("correlNorm",
                                              height = "100%",
                                              width = "100%") %>% withSpinner(color="#0dc5c1", proxy.height = "400px")
                                   
                                   )                                   
                        )
                      
                        
                      )
                      
             ),
             
             
             
             ################################################################################################################################
             #PCA
             ################################################################################################################################
             
             tabPanel("PCA", value = "panel4", icon = icon("fas fa-cube"),
                      
                      sidebarPanel(
                        
                        h2(strong("PCA")),
                        
                        h5("Principal component analysis (PCA) shows the similarity between the samples 
                           and is particularly useful for the detection of outliers."),
                        
                        hr(),
                
                        materialSwitch(
                          inputId = "plot3d",
                          label = "3D",
                          value = FALSE, 
                          status = "danger"),
                        
                        selectInput(inputId = "xpca", 
                                    label = "x-axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", "PC7", "PC8"),
                                    selected = "PC1"),
                        
                        selectInput(inputId = "ypca", 
                                    label = "y-axis",
                                    choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", "PC7", "PC8"),
                                    selected = "PC2"),
                        
                        uiOutput("uizpca"),
                        
                        hr(),
                        
                        actionBttn(inputId = "pca.back",
                                   label = "Back",
                                   style = "jelly",
                                   color = "danger",
                                   icon = icon("arrow-left")),
                        
                        actionBttn(inputId = "pca.ok",
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary",
                                   icon = icon("arrow-right"))
                      ),
                      
                      mainPanel(
                        uiOutput("uipca"),
                        
                        uiOutput("uipca3d"),
                        
                        plotlyOutput("variances") %>% withSpinner(color="#0dc5c1")
                      )
                      
                      
             ),
             
             ################################################################################################################################
             #Statistical analysis
             ################################################################################################################################
             
             tabPanel("Statistical analysis", value = "panel5", icon = icon("fas fa-chart-bar"),
                      
                      navlistPanel(
                        
                        "Statistical analysis",
                        
                        #Top table
                        tabPanel("Top table", value = "panel6",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "Top table"))),
                                 
                                 hr(),
                                 
                                 fluidRow(

                                   column(4,
                                          uiOutput("uitoptablecomp")),
                                   
                                   column(4,
                                          uiOutput("uitoptabledir"))
                                 ),
                                 
                                 hr(),
                                 
                                 dataTableOutput("finaltable") %>% withSpinner(color="#0dc5c1"),
                                 
                                 downloadButton("downloadfinal", "Download top table")
                                 
                                 ),
                        
                        
                        #P-value analysis
                        tabPanel("P-value analysis", value = "panel7",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "P-value analysis"))),
                                 
                                 hr(),
                                 
                                 uiOutput("uipcomp"),
                                 
                                 hr(),
                                 
                                 plotlyOutput("phist") %>% withSpinner(color="#0dc5c1")
                                 
                                 ),
                        
                        
                        #logFC analysis
                        tabPanel("Fold change analysis", value = "panel8",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "Fold change analysis"))),
                                 
                                 hr(),
                                 
                                 uiOutput("uilogfccomp"),
                                 
                                 hr(),
                                 
                                 plotlyOutput("logfchist") %>% withSpinner(color="#0dc5c1")
                                 
                                 ),
                        
                        
                        #Volcano plot
                        tabPanel("Volcano plot", value = "panel9",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "Volcano plot"))),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   column(3,
                                          
                                          prettyRadioButtons(
                                            inputId = "raworadj",
                                            label = "P-value", 
                                            choices = c("Raw P-value" = "raw", "Adjusted P-value" = "adj")),
                                          
                                          actionBttn(inputId = "volcano.ok", 
                                                     label = "Plot",
                                                     style = "jelly",
                                                     color = "primary"),
                                          
                                   ),
                                   
                                   column(4,
                                          
                                          numericInput(
                                            inputId = "pthreshold",
                                            label = "P threshold",
                                            value = 0.05),
                                          
                                          numericInput(
                                            inputId = "logfcthreshold",
                                            label = "logFC threshold",
                                            value = 1)
                                   ),
                                   
                                   column(4,
                                          
                                          uiOutput("uivolcanocomp"),
                                          
                                          uiOutput("uivolcanodir")
                                
                                   )
                                   
                                 ),
                                 
                                 hr(),
                                 
                                 plotlyOutput("volcano") %>% withSpinner(color="#0dc5c1"),
                                 
                                 dataTableOutput("voltable") %>% withSpinner(color="#0dc5c1"),
                                 
                                 downloadButton("downloadvol", "Download table")
                                 ),
                        
                        
                        tabPanel("Venn diagram", value = "panel10",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "Venn Diagram"))),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   
                                   column(3,
                                      
                                          prettyRadioButtons(
                                            inputId = "raworadjvenn",
                                            label = "P-value", 
                                            choices = c("Raw P-value" = "raw", "Adjusted P-value" = "adj"))
                                          
                                   ),
                                   
                                   
                                   column(4,
                                          
                                          numericInput(
                                            inputId = "pthresholdvenn",
                                            label = "P threshold",
                                            value = 0.05)
                                          
                                   ),
                                   
                                   column(4,
                                          
                                          numericInput(
                                            inputId = "logfcthresholdvenn",
                                            label = "logFC threshold",
                                            value = 0)
                                          
                                   )
                                   
                  
                                   
                                 ),
                                 
                                 hr(),
                                 
                                 uiOutput("uicontrast"), 
                         
                                 htmlOutput("legendvenn"),
                                 
                                 br(),
                                 
                                 plotOutput("venndiagram")
                                 ),
                        
                        #GO analysis
                        tabPanel("GO analysis", value = "panel11",
                                 
                                 h1(strong(span(style = "color:#3A3B3C", "GO enrichment analysis"))),
                                 
                                 hr(),
                                 
                                 fluidRow(
                                   
                                   column(3,
                                          
                                          prettyRadioButtons(
                                            inputId = "raworadjgoa",
                                            label = "P-value", 
                                            choices = c("Raw P-value" = "raw", "Adjusted P-value" = "adj")),
                                          
                                          br(),
                                          
                                          actionBttn(inputId = "goa.ok", 
                                                     label = "Calculate",
                                                     style = "jelly",
                                                     color = "primary")
                                          
                                   ),
                                   
                                   
                                   column(4,
                                          
                                          numericInput(
                                            inputId = "pthresholdgoa",
                                            label = "P threshold",
                                            value = 0.05),
                                          
                                          numericInput(
                                            inputId = "logFCthresholdgoa",
                                            label = "logFC threshold",
                                            value = 0)
                                          
                                   ),
                                   
                                   column(4,
                                          
                                          pickerInput(
                                            inputId = "ontology",
                                            label = "Ontology", 
                                            choices = c("Biological process" = "BP", 
                                                        "Molecular function" = "MF",
                                                        "Cellular component" = "CC"),
                                            options = list(
                                              style = "btn-danger")),
                                          
                                          uiOutput("uigoacomp")
                                          
                                   )
                
                                 ),
                                 
                                 hr(),
                                 
                                 dataTableOutput("goa" )%>% withSpinner(color="#0dc5c1"),
                                 
                                 plotOutput("GOplot",
                                            height = "100%",
                                            width = "100%")
                                 
                        )
                        
                      )
                      
                      
             ),
             
             ################################################################################################################################
             #Help
             ################################################################################################################################
             
             tabPanel("Documentation", value = "panel12", icon = icon("far fa-question-circle"),
                      
                      
                      
                      
             )
             
             
             
             
  )
)
)

