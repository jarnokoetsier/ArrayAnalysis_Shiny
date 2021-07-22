
###############################################################################################################################

#USER INTERFACE

###############################################################################################################################

ui <- fluidPage(
  
  useSweetAlert(),
  
  navbarPage("ArrayAnalysis", id = "navbar",
             ################################################################################################################################
             #Enter database accession number
             ################################################################################################################################
             
             tabPanel("Data accession", value = "panel1",
                      
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
                                            color = "royal"),
                                 
                                 actionBttn(inputId = "downloadGEO", 
                                            label = "Next",
                                            style = "jelly",
                                            color = "primary"),
                                 
                                 actionBttn(inputId = "infopanel1", 
                                            label = NULL,
                                            style = "simple",
                                            color = "success",
                                            icon = icon("info")),
                               
                               br(),
                               
                                 
                                br() 
                               )
                      ),
                      
                      
                    
                      
             ),
             
             
             
             
             
             ################################################################################################################################
             #get meta data
             ################################################################################################################################
             
             tabPanel("Meta data", value = "panel2",
                      
                      
                      
                      sidebarPanel(
                        
                        radioGroupButtons(
                          inputId = "makegroups",
                          label = "Make sample groups from:", 
                          choices = c("Dataset", "Description file", "Manual grouping"),
                          selected = "Dataset",
                          status = "danger"),
                        
                        uiOutput("uidescriptionfile"),
                        
                        uiOutput("uitextdescription"),
                        
                        uiOutput("groups"),
                        
                        uiOutput("owngroupsout"),
                        
                        actionBttn(inputId = "meta.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary"),
                        
                        actionBttn(inputId = "infopanel2", 
                                   label = NULL,
                                   style = "simple",
                                   color = "success",
                                   icon = icon("info"))
                        
                        
                      ),
                      
                      mainPanel(
                        
                        dataTableOutput(outputId = "grouping") %>% withSpinner(color="#0dc5c1"),
                        
                        downloadButton("downloadmeta", "Download meta table")
                      )
                      
                      
                      
             ),
             
             ################################################################################################################################
             #Pre-processing
             ################################################################################################################################
             
             tabPanel("Pre-processing", value = "panel3",
                      
                      sidebarPanel(
                        
                        awesomeCheckbox(inputId = "outlier",
                                        label = "Keep all samples",
                                        value = TRUE,
                                        status = "danger"),
                        
                        uiOutput("outliersout"),
                        
                        radioGroupButtons(
                          inputId = "rma",
                          label = "Pre-processing",
                          choices = c("RMA", "Advanced"),
                          selected = "RMA",
                          status = "primary"
                        ),
                        
                        uiOutput("uibg"),
                        
                        uiOutput("uinorm"),
                        
                        uiOutput("uibgcorrectmethod"),
                        
                        uiOutput("uinormmethod"),
                        
                        uiOutput("uipmcorrectmethod"),
                        
                        uiOutput("uisummarymethod"),
                        
                        actionBttn(inputId = "ann.ok",
                                   label = "Calculate",
                                   style = "fill",
                                   color = "danger"),
                        
                        br(),
                        
                        uiOutput("proceedann"),
                        
                        
                      ),
                      
                      mainPanel(
                        
                        tabsetPanel(
                          tabPanel("Density plot", 
                                   plotOutput("normhist") %>% withSpinner(color="#0dc5c1")),
                          tabPanel("Boxplot", 
                                   plotOutput("normboxplot") %>% withSpinner(color="#0dc5c1"))
                        )
                      
                        
                      )
                      
             ),
             
             
             
             ################################################################################################################################
             #PCA
             ################################################################################################################################
             
             tabPanel("PCA", value = "panel4",
                      
                      sidebarPanel(
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
                        
                        
                        actionBttn(inputId = "pca.ok", 
                                   label = "Next",
                                   style = "jelly",
                                   color = "primary")
                      ),
                      
                      mainPanel(
                        uiOutput("uipca"),
                        
                        uiOutput("uipca3d"),
                        
                        plotlyOutput("variances") %>% withSpinner(color="#0dc5c1")
                      )
                      
                      
             ),
             
             ################################################################################################################################
             #differential isoform expression
             ################################################################################################################################
             
             tabPanel("Differential expression analysis", value = "panel5",
                      
                      navlistPanel(
                        
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
                        tabPanel("Fold change analysis", value = "panel8"),
                        
                        
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
                                 
                                 uiOutput("uicontrast"),
                                 
                                 hr(),
                                 
                                 plotOutput("venndiagram"))
                      )
                      
                      
             )
             
             
             
             
  )
)


