

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

server <- function(input, output, session){
  options(shiny.maxRequestSize=125*1024^2)
  
  ##############################################################################
  #Data selection
  ##############################################################################
  
  
  #Hide tabs
  hideTab("navbar", target = "panel2")
  hideTab("navbar", target = "panel3")
  hideTab("navbar", target = "panel3b")
  hideTab("navbar", target = "panel4")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel6")
  
  
  #Welcome message
  #sendSweetAlert(
  #  session = session,
  #  title = "ArrayAnalysis",
  #  text = "Welcome to the ArrayAnalysis app!",
  #  btn_labels = "Click here to start"
  #)
  
  
  #Info panel
  observeEvent(input$infopanel1, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "ArrayAnalysis has been developed for the user-friendly analysis 
      of microarray data. Click on the documentation tab for more information.",
      type = "info"
    )
  })
  
  
  #****************************************************************************#
  #Saved data
  #****************************************************************************#
  
  #Enter ID to retrieve saved data
  observeEvent(input$continue, {
    inputSweetAlert(
      session = session,
      inputId = "idcontinue",
      title = "Enter ID",
      text = "After you saved the data, you received an ID. 
      This ID is required to retrieve the saved data.",
      input = "text",
      inputPlaceholder = "e.g. aa123456789",
      btn_labels = c("Continue", "Cancel")
    )
    
  })
  
  #Success/Error message
  observe({
    req(input$idcontinue)
    
    if (file.exists(paste0("meta1", input$idcontinue))){
      showTab("navbar", target = "panel3")
      
      updateNavbarPage(session, "navbar",
                       selected = "panel3")
      
    }
    
    if (!file.exists(paste0("meta1", input$idcontinue))){
      sendSweetAlert(
        session = session,
        title = "Error!!",
        text = "Invalid ID. Please try again.",
        type = "error")
      
    }
    
  })
  
  
  #get saved meta data
  meta1.saved <- reactive({
    meta1 <- NULL
    
    if(length(input$idcontinue) > 0){
      if (file.exists(paste0("meta1", input$idcontinue))){
        meta1 <- readRDS(paste0("meta1", input$idcontinue))
      }
    }
    return(meta1)
  })
  
  #Get saved normalized expression data
  normData.saved <- reactive({
    normData <- NULL
    
    if(length(input$idcontinue) > 0){
      if (file.exists(paste0("normData", input$idcontinue))){
        normData <- readRDS(paste0("normData", input$idcontinue))
      }
      
    }
    return(normData)
  })
  
  #Get saved info data
  info.saved <- reactive({
    info <- NULL
    
    if(length(input$idcontinue) > 0){
      if (file.exists(paste0("info", input$idcontinue))){
        info <- readRDS(paste0("info", input$idcontinue))
      }
      
    }
    return(info)
  })
  

  #****************************************************************************#
  #   Data from database
  #****************************************************************************#
  
  #Select database accession
  output$getdatabaseout <- renderUI({
    req(input$database)
    if (input$database != "Upload"){
      textInput(inputId = "dbAccession", 
                label = NULL, 
                value = "",
                placeholder = "Accession number")
      
    }
  })
  
  
  #Example dataset (GSE36980)
  observeEvent(input$example, {
    
    if (input$database == "GEO"){
      updateTextInput(session, "dbAccession",
                      label = NULL,
                      value = "GSE36980")
      
    }
    
    if (input$database == "ArrayExpress"){
      updateTextInput(session, "dbAccession",
                      label = NULL,
                      value = "E-MTAB-9988")
    }
    
    
  })
  
  
  #Get database
  database <- eventReactive(input$startAnalysis, {
    input$database
    
  })
  
  #Get accession
  accession <- eventReactive(input$startAnalysis, {
    input$dbAccession
  })
  
  
  #Get raw or normalized data
  rawornormalized <- reactive({
    if (is.null(normData.saved())){
      rawornormalized <- input$rawornormalizeddata
    }
    
    if (!is.null(normData.saved())){
      rawornormalized <- "Normalized"
    }
    return(rawornormalized)
  })
  

  #****************************************************************************#
  #   CEL files (raw data)
  #****************************************************************************#
  
  
  #Upload CELs (raw data)
  output$uploadcelsout <- renderUI({
    req(input$database)
    if (input$database == "Upload"){
      if (input$rawornormalizeddata == "Raw"){
        fileInput(inputId = "uploadcelsin",
                  label = NULL,
                  accept = ".zip",
                  placeholder = "Select .zip data file")
      }
    }
  })
  
  #Get data path.
  zzip <- eventReactive(input$startAnalysis, {
    req(input$uploadcelsin)
    input$uploadcelsin$datapath
  })
  
  
  #****************************************************************************#
  #   txt file (normalized data)
  #****************************************************************************#
  
  #Upload txt (normalized data)
  output$uploadtxtout <- renderUI({
    req(input$database)
    if (input$database == "Upload"){
      if (input$rawornormalizeddata == "Normalized"){
        fileInput(inputId = "uploadtxtin",
                  label = NULL,
                  accept = ".txt",
                  placeholder = "Select .txt data file")
      }
    }
  })
  
  gxData <- eventReactive(input$startAnalysis, {
    req(input$uploadtxtin$datapath)
    gxData <- read.delim(input$uploadtxtin$datapath, 
                         as.is = TRUE,
                         row.names = 1)
    gxData <- as.matrix(gxData)
    return(gxData)
  })
  
  
  #****************************************************************************#
  #   Get Gene ExpressionSet
  #****************************************************************************#
  
  #Required for gene expression data of normalized data only
  #Required for meta data of normalized and raw data
  
  #Download data for retrieval of meta data
  gset <- eventReactive(input$startAnalysis, {
    
    #Loading message
    if (database() != "Upload"){
      showModal(modalDialog(
        title = h4(paste0(accession(), " dataset from the ", database(), 
                          " database is being loaded...."), align = "center"), 
        footer = NULL,
        h5("Please be patient. This might take a while.", align = "center")
      ))
    }
    
    #GEO  
    if (database() == "GEO"){
      gset <- NULL
      if (grepl("GSE", accession())){
        tryCatch({
          gset <- getGEO(accession(), GSEMatrix =TRUE, getGPL = FALSE)
        },
        error=function(cond) {
          return(NULL)
        })
      }
    }
    
    #ArrayExpress
    if (database() == "ArrayExpress"){
      gset <- NULL
      if (grepl("E-MTAB-", accession()) | grepl("E-GEOD-", accession())){
        tryCatch({
          gset <- ArrayExpress(accession())
        },
        error=function(cond) {
          return(NULL)
        })
        
        if (!is.null(gset)){
          gset <- rma(gset)
        }
      }
    }
    
    #Upload
    if (database() == "Upload"){
      if (rawornormalized() == "Raw"){
        gset <- "empty"
      }
      
      if (rawornormalized() == "Normalized"){
        
        gset <- new(Class = "ExpressionSet", exprs = gxData())
      }
      
    }
    
    return(gset)
  })
  
  
  #****************************************************************************#
  #   Get CEL files (only for raw data)
  #****************************************************************************#

  celfiles <- eventReactive(input$startAnalysis, {
    
    celfiles <- NULL
    database <- database()
    
    if ((rawornormalized() == "Raw") & (!is.null(gset()))){
      if (database == "GEO"){
        
        accession <- accession()
        exdir <- paste0("data_", accession)
        
        if (!file.exists(exdir)){
          
          getGEOSuppFiles(accession)
          tarfile <- paste0(accession, "/", accession, "_RAW.tar")
          untar(tarfile, exdir = exdir)
          
        }
        
        
        celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", 
                              full.names = TRUE)
        
      }
      
      if (database == "ArrayExpress"){
        
        accession <- accession()
        exdir <- paste0("data_", accession)
        
        if (!file.exists(exdir)) {
          getAE(accession, type = "raw", extract = FALSE)
          unzip(paste0(accession, ".raw.1.zip"), exdir = exdir)
        }
        
        
        celfiles <- list.files(paste0(exdir, "/"), pattern = "CEL")
        
      }
      
      if (database == "Upload") {
        zzip <- zzip()
        celfiles <- unzip(zzip)
        
      }
    }
    
    if ((rawornormalized() == "Normalized") & (!is.null(gset()))){
      celfiles <- "empty"
    }
    
    return(celfiles)
    
  })
  

  #****************************************************************************#
  #Data successfully downloaded
  #****************************************************************************#

  #Success message
  observeEvent(if (length(celfiles())>0){input$startAnalysis}, {
    
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "Success!!",
      text = "Dataset successfully selected! Now you can group the samples.",
      type = "success")
    
  })
  
  
  #go to next tab
  observeEvent(if (length(celfiles())>0){input$startAnalysis}, {
    showTab("navbar", target = "panel2")
  })
  
  
  observeEvent(if (length(celfiles())>0){input$startAnalysis}, {
    updateNavbarPage(session, "navbar",
                     selected = "panel2")
    
  })
  
  
  
  #****************************************************************************#
  #Incorrect data selected
  #****************************************************************************#
  
  #Error message
  observeEvent(if (length(celfiles())==0){input$startAnalysis}, {
      
    removeModal()
    
    if ((rawornormalized() == "Raw")){
      sendSweetAlert(
        session = session,
        title = "Error!",
        text = "Raw data does not exist for this dataset. 
        Consider using normalized data.",
        type = "error")
      
    } else{
      sendSweetAlert(
        session = session,
        title = "Error!",
        text = "No valid dataset selected. Please try again.",
        type = "error")
    }
  
  })
  
  


  
  
  
  ##############################################################################
  #Grouping of the data
  ##############################################################################
  
  #****************************************************************************#
  #Get samples
  #****************************************************************************#
  
  #This information will be used to allow manual grouping as well as a final
  #quality check.
  
  samplelist <- reactive({
    
    if(database() == "GEO"){
      samplelist <- gset()[[paste0(accession(),"_series_matrix.txt.gz")]]@phenoData@data[["geo_accession"]]
    }
    
    if(database() == "ArrayExpress"){
      samplelist <- rownames(gset()@phenoData@data)
    }
    
    if (database() == "Upload"){ 
      if (rawornormalized() == "Raw"){
        samplelist <- gsub(".*/","",celfiles())
      }
      
      if (rawornormalized() == "Normalized"){
        samplelist <- colnames(gxData())
      }
    }
    
    return(samplelist)
  })
  
  
  #****************************************************************************#
  #Options
  #****************************************************************************#

  #Options for grouping the samples
  
  #Use can select either retrieve the metadata from the dataset, description 
  #file or manual input.
  
  output$makegroupsui <- renderUI({
    
    if (length(input$paired) > 0){
      
      if ((database() != "Upload") & (!input$paired)){
        choices = c("Dataset", "Description file", "Manual grouping")
        selected = "Dataset"
      }
      
      if ((database() == "Upload") & (!input$paired)){
        choices = c("Description file", "Manual grouping")
        selected = "Description file"
      } 
      
      if ((database() != "Upload") & (input$paired)){
        choices = c("Dataset", "Description file")
        selected = "Dataset"
      }
      
      if ((database() == "Upload") & (input$paired)){
        choices = "Description file"
        selected = "Description file"
      }   
      
      radioGroupButtons(
        inputId = "makegroups",
        label = "Make sample groups from:", 
        choices = choices,
        selected = selected,
        status = "danger")
      
    }

  })
 
  
  #Upload description file
  output$uidescriptionfile <- renderUI({
    
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Description file") {
        
        fileInput(inputId = "descriptionfile",
                  label = NULL,
                  accept = ".txt",
                  placeholder = "Select description file")
        
      }
    }
  })
  
  #Explanation of description file
  output$textdescription <- renderText({
    "<p> The description file needs to be an txt file containing at least two 
    columns:<br /> 
    <p> 1) The first column should contain the <b> sample IDs</b>. Note that 
    these IDs should match to the CEL file names.<br />
    <p> 2) The other columns should include the <b> grouping </b> of the 
    samples. <br />
    <p> <b> IMPORTANT: </b> All columns should have an header. <br />
    <br />"
  })
  
  #Add explanation to user interface
  output$uitextdescription <- renderUI({
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Description file") {
        
        htmlOutput("textdescription")
        
      }
    }
    
  })
  
  
  #If dataset: provide possible grouping variable(s)
  output$groups <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        checkboxGroupInput(inputId = "groupselect", 
                           label = "Select grouping variable(s):", 
                           choices = 
                             colnames(
                               get_grouping(gset(), database = database())),
                           selected = auto_group(
                             get_grouping(gset(), database = database())))
      }
    }
   
  })
  
  #If description file: provide possible grouping variable(s)
  output$groups1 <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Description file"){
        req(input$descriptionfile)
        checkboxGroupInput(inputId = "groupselect1", 
                           label = "Select grouping variable(s):", 
                           choices = colnames(
                             read.table(input$descriptionfile$datapath, 
                                        header = TRUE))[-1],
                           selected = colnames(
                             read.table(input$descriptionfile$datapath, 
                                        header = TRUE))[2])
      }
    }
    
  })

  
  #If dataset: provide possible pairing variable(s)
  output$pairs <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        if (input$paired) {
          checkboxGroupInput(inputId = "pairselect", 
                             label = "Select pairing variable(s):", 
                             choices = colnames(
                               get_grouping(gset(), database = database())),
                             selected = auto_group(
                               get_grouping(gset(), database = database())))
        }
      }
    }
    
  })
  
  
  #If description file: provide possible pairing variable(s)
  output$pairs1 <- renderUI({
    
    if (length(input$makegroups) > 0){
    
        if (input$makegroups == "Description file"){
          if (input$paired) {
            req(input$descriptionfile)
            checkboxGroupInput(inputId = "pairselect1", 
                               label = "Select pairing variable(s):", 
                               choices = colnames(
                                 read.table(input$descriptionfile$datapath, 
                                            header = TRUE))[-1],
                               selected = colnames(
                                 read.table(input$descriptionfile$datapath, 
                                            header = TRUE))[2])
          }
          
        }
      
    }
    
  })
  
  
  #Multi-input for manual grouping
  output$owngroupsout <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Manual grouping") {
        
        multiInput(
          inputId = "owngroupsin", label = NULL,
          choices = samplelist(),
          width = "500px", 
          options = list(
            enable_search = TRUE,
            non_selected_header = "Group 1",
            selected_header = "Group 2"))
        
      }
    }
  })
  

  
  
  #****************************************************************************#
  # Metadata dataframe
  #****************************************************************************#
  
  #Make "metadata" dataframe
  meta <- reactive({
    
    meta <- NULL
    
    #Metadata from dataset
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        if (length(input$groupselect) > 0){
          
          #Independent samples
          if (!input$paired){
            groups <- get_grouping(gset(), database = database())
            meta <- get_meta(gset = gset(), 
                             grouping_column = groups[,input$groupselect], 
                             database = database())
          }
          
          #Paired samples
          if (input$paired){
            groups <- get_grouping(gset(), database = database())
            meta <- get_meta(gset = gset(), 
                             grouping_column = groups[,input$groupselect], 
                             pairing_column = groups[, input$pairselect], 
                             database = database())
          }

        }
      }
      
      #Metadata from manual grouping
      if (input$makegroups == "Manual grouping"){
        
        group2 <- cbind(input$owngroupsin, rep("group 2", 
                                               length(input$owngroupsin)))
        
        group1 <- setdiff(samplelist(), input$owngroupsin)
        group1 <- cbind(group1, rep("group 1", length(group1)))
        
        meta <- as.data.frame(rbind(group1, group2))
        rownames(meta) <- NULL
        colnames(meta) <- c("Sample.ID", "Grouping")
      }
      
      #Metadata from description file
      if (input$makegroups == "Description file"){
        if (length(input$descriptionfile$datapath) > 0){
          req(input$groupselect1)
          
          #Read description file
          dscfile <- read.table(input$descriptionfile$datapath, header = TRUE)
          
          if (length(input$groupselect1) > 0){
            
            #Samples from description file
            samples <- dscfile[,1]
            
            #Grouping from description file
            grouping <- dscfile[,input$groupselect1]
            if (length(input$groupselect1) > 1){
              grouping <- unite(grouping, col = grouping, sep = ".")
            }
            
            #Independent samples
            if (!input$paired){
              meta <- cbind.data.frame(samples,grouping)
              colnames(meta) <- c("Sample.ID", "Grouping")
            }
            
            #Dependent samples
            if (input$paired){
              if (length(input$pairselect1) > 0){
                pairing <- dscfile[, input$pairselect1]
                if (length(input$pairselect1) > 1){
                  pairing <- unite(pairing, col = pairing, sep = ".")
                }
                meta <- cbind.data.frame(samples,grouping, pairing)
                colnames(meta) <- c("Sample.ID", "Grouping", "Pairing")
              }
            }
            
            rownames(meta) <- NULL
          }
 
        }
      }
    }
    
    return(meta)
  })
  
  
  
  #****************************************************************************#
  #Output
  #****************************************************************************#
  
  #Table of metadata
  output$grouping <- renderDataTable({
    
    meta()
    
  })
  
  
  #Metadata heatmap
  output$sampleheatmap <- renderPlotly(NULL)
  
  output$sampleheatmap <- renderPlotly({
    if (input$makegroups != "Manual grouping"){
      if (input$makegroups == "Dataset") {
        groups <- get_grouping(gset(), database = database())
      }
      
      if (input$makegroups == "Description file") {
        groups <- read.table(input$descriptionfile$datapath, 
                             header = TRUE, row.names = 1)
      }
      
      groups1 <- groups
      for (i in 1:ncol(groups)){
        groups1[,i] <- factor(groups[,i])
      }
      
      groups2 <- sapply(groups1, unclass)
      
      
      for (i in 1:ncol(groups2)){
        groups2[,i] <- groups2[,i]/(max(groups2[,i]))
      }
      
      rownames(groups2) <- rownames(groups1)
      
      plot_ly(z = groups2, x = colnames(groups2), y = rownames(groups2), 
              type = "heatmap", showscale = F, 
              text = as.matrix(groups), hoverinfo = c("x", "y", "text"))
    }
  })
  
  
  
  #Bar charts of metadata
  
  #If dataset: 
  
  # provide possible grouping variable(s)
  output$uimetaBar1 <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        pickerInput(
          inputId = "xselect",
          label = "Bars",
          choices = colnames(
            get_grouping(gset(), database = database())),
          selected = auto_group(
            get_grouping(gset(), database = database())),
          options = list(
            style = "btn-primary"))
      }
    }

  })
  
  # provide possible colour variable(s)
  output$uimetaBar2 <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        pickerInput(
          inputId = "fillselect",
          label = "Colour",
          choices = colnames(
            get_grouping(gset(), database = database())),
          selected = auto_group(
            get_grouping(gset(), database = database())),
          options = list(
            style = "btn-info"))
      }
    }
    
  })
  
  #If description file: 
  
  # provide possible grouping variable(s)
  output$uimetaBar3 <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Description file"){
        req(input$descriptionfile)
        
        pickerInput(
          inputId = "xselect",
          label = "Bars",
          choices = colnames(
            read.table(input$descriptionfile$datapath, 
                       header = TRUE))[-1],
          selected = colnames(
            read.table(input$descriptionfile$datapath, 
                       header = TRUE))[2],
          options = list(
            style = "btn-primary"))

      }
    }
    
  })
  
  # provide possible colour variable(s)
  output$uimetaBar4 <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Description file"){
        req(input$descriptionfile)
        
        pickerInput(
          inputId = "fillselect",
          label = "Colour",
          choices = colnames(
            read.table(input$descriptionfile$datapath, 
                       header = TRUE))[-1],
          selected = colnames(
            read.table(input$descriptionfile$datapath, 
                       header = TRUE))[2],
          options = list(
            style = "btn-info"))
        
      }
    }
    
  })

    
  output$metaBar <- renderPlot(NULL)
  
  output$metaBar <- renderPlot({
    req(input$xselect)
    req(input$fillselect)
    
    if (input$makegroups != "Manual grouping"){
      if (input$makegroups == "Dataset") {
        groups <- get_grouping(gset(), database = database())
      }
      
      if (input$makegroups == "Description file") {
        groups <- read.table(input$descriptionfile$datapath, 
                             header = TRUE, row.names = 1)
      }
      
      barPlot <- ggplot(data = groups, aes(x = groups[,input$xselect],
                                           fill = groups[,input$fillselect])) +
        geom_bar(color = "black") +
        ylab("Count") +
        xlab(NULL) +
        labs(fill = " ") +
        theme_classic() +
        theme(legend.position="bottom") + 
        scale_fill_brewer(palette="Dark2")
      
      return(barPlot)
      
      
    }
    
    
  })
  
  
  #****************************************************************************#
  #Error/success message
  #****************************************************************************#
  
  observeEvent(input$meta.ok, {
    samplelist <- samplelist()
    samplelist <- as.data.frame(samplelist)
    colnames(samplelist) <- "names"
    
    #Match sample names to CEL name
      test <- samplelist %>% 
        fuzzyjoin::fuzzy_inner_join(meta(), by = c("names" = "Sample.ID"), 
                                    match_fun = str_detect)
      
      if (nrow(test) < length(samplelist)){
        test <- meta() %>% 
          fuzzyjoin::fuzzy_inner_join(samplelist, by = c("Sample.ID" = "names"), 
                                      match_fun = str_detect)
        
        rownames(test) <- test$names
        test <- test[samplelist$names,]
        rownames(test) <- NULL
      }
    
      #Display error if the matching is not successful
      if (nrow(test) < length(samplelist)){
        sendSweetAlert(
          session = session,
          title = "Error!!",
          text = "Something went wrong. Please try again.",
          type = "error"
        )
      } else{
        sendSweetAlert(
          session = session,
          title = "Success!!",
          text = "Samples are successfully grouped! 
          Now you can start pre-processing the data.",
          type = "success"
        )
        
        
        updateNavbarPage(session, "navbar",
                         selected = "panel3")
        
        showTab("navbar", target = "panel3")
        
        
      }

  })
  
  #download table of meta data
  output$downloadmeta <- downloadHandler(
    filename = "MetaData.txt",
    content = function(file){
      write.table(meta(), file)
    }
  )
  
  #Go to previous tab
  observeEvent(input$meta.back, {
    updateNavbarPage(session, "navbar",
                     selected = "panel1")
    
    
  })
  
  
  ##############################################################################
  #   Pre-processing
  ##############################################################################

  
  #****************************************************************************#
  #Options and messages
  #****************************************************************************#
 
  #No raw data if continued with saved data
   output$alreadynormalized <- renderText({
    
    if (rawornormalized() == "Normalized"){
      
      "You continued with normalized data. So, no raw data available."
    }
  })
  
  
  #Select outliers for removal
  output$outliersout <- renderUI({
    
    if (length(input$outlier) > 0){
      if(input$outlier == FALSE){
        samples <- meta()[,1]
        pickerInput(inputId = "outliersin", 
                    label = "Select samples to be removed", 
                    choices = as.vector(meta()[,1]),
                    multiple = TRUE)
        
      }
    }
  })

  
  #Identify outliers
  outliers <- eventReactive(input$preprocessing, {
    if (input$outlier == TRUE){
      outliers <- NULL
    }
    
    if (input$outlier == FALSE){
      if (length(input$outliersin) < 1){
        outliers <- NULL
      }
      if (length(input$outliersin) > 0){
        outliers <- input$outliersin
      }
    }
    
    return(outliers)
  })
  
  
  #****************************************************************************#
  #Raw data
  #****************************************************************************#
  
  #Get raw data:
  #1. remove outliers
  #2. read CEL files
  rawData <- eventReactive(input$preprocessing, {
    rawData <- NULL
    
    if (rawornormalized() == "Raw"){
      
      #Show loading message
      showModal(modalDialog(title = h4("Dataset is being pre-processed....", 
                                       align = "center"), 
                            footer = NULL,
                            h5("Please be patient. This might take a while.", 
                               align = "center"),
                            br()))
      
      
      #Remove outliers
      outliers = outliers()
      celfiles = celfiles()
      
      if (!is.null(outliers)){
        
        outliers_select <- outliers[1]
        
        if (length(outliers) > 1){
          for (i in 1:(length(outliers)-1)){
            outliers_select <- paste(outliers_select, outliers[1+i], sep = "|")
          }
          
        }
        
        celfiles = celfiles[-grep(outliers_select, celfiles)]
      }
      
      #Read CEL files
      rawData = ReadAffy(filenames=celfiles)
      
    }

    return(rawData)
    
  })
  
  #****************************************************************************#
  #Metadata
  #****************************************************************************#
  
  #meta1: CEL names + sample names + grouping (+ pairing)
  meta1 <- reactive({
    
    #If not continued with saved data
    if (is.null(meta1.saved())){
      
      if (rawornormalized() == "Raw"){
        rawData <- rawData()
      }
      
      if (rawornormalized() == "Normalized"){
        rawData <- normData()
      }
      
      meta <- meta()
      
      #add CEL names to meta
      names <- as.data.frame(sampleNames(rawData))
      
      colnames(names) <- "names"
      
      meta1 <- names %>% 
        fuzzyjoin::fuzzy_inner_join(meta, by = c("names" = "Sample.ID"), 
                                    match_fun = str_detect)
    }
    
    #if continued with saved data
    if (!is.null(meta1.saved())){
      meta1 <- meta1.saved()
    }

    return(meta1)
    
  })
  
  
  #Experiment factor: required for colouring of plots
  experimentFactor <- reactive({
    meta1 <- meta1()
    experimentFactor <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
    return(experimentFactor)
  })
  
  
  
  
  #****************************************************************************#
  #Normalization options
  #****************************************************************************#
  
  #aType
  aType <- reactive({
    Data <- rawData()
    aType <- "PMMM"  
    
    # Test whether the dataset is of aType "PM-only"
    mismatches <- mm(Data[,1])  
    
    if(is.null(mismatches)) {
      # mm does not exist
      aType <- "PMonly"
    } else {
      if(sum(is.na(mismatches))>=length(mismatches)/2){ 
        # mm is always NA or there are more NA values in the mm probes than 
        # defined values (assuming these would just be controls)
        aType <- "PMonly"
      } else {
        matches <- pm(Data[,1])
        notNA <- !is.na(mismatches) & !is.na(matches)
        if(sum(mismatches[notNA]!=matches[notNA])==0){
          # MM contains a copy of PM, which indicates a PMonly array
          aType <- "PMonly"
        }
      }
    }
    return(aType)
  })
  
  #Probeset annotation
  PSannotation <- reactive({
    if (is.null(info.saved())){
      req(input$annotations)
      PSannotation <- input$annotations
    }
    
    if (!is.null(info.saved())){
      PSannotation <- info.saved()[[1]]
    }
    return(PSannotation)
  })
  
  #CDF type
  CDFtype <- reactive({
    if (is.null(info.saved())){
      req(input$CDFtype)
      CDFtype <- input$CDFtype
    }
    
    if (!is.null(info.saved())){
      CDFtype <- info.saved()[[2]]
    }
    return(CDFtype)
  })
  
  #Species
  species <- reactive({
    if (is.null(info.saved())){
      req(input$species)
      species <- input$species
    }
    
    if (!is.null(info.saved())){
      species <- info.saved()[[3]]
    }
    return(species)
  })
  
  
  #****************************************************************************#
  #Normalization
  #****************************************************************************#
  
  #Normalize data
  normData <- eventReactive(((input$preprocessing) | (if(rawornormalized() == "Normalized"){input$meta.ok})), {
    
    
    #Start from raw data
    #..........................................................................#
    if (rawornormalized() == "Raw"){
      
      normData <- tryCatch({
        normMeth <- input$normMeth
        CDFtype <- input$CDFtype
        species <- input$species
        
        experimentFactor <- experimentFactor()
        
        ifelse(input$annotations=="Custom annotations",
               customCDF<-TRUE,customCDF<-FALSE)
        ifelse(input$annotations=="Upload annotation file",
               uploadCDF<-TRUE,uploadCDF<-FALSE)
        ifelse(input$perGroup == "Use all arrays",
               perGroup<-FALSE,perGroup<-TRUE)
        
        normMeth <- toupper(normMeth)
        #if customCDF option is chosen, apply to copy of Data, in order not to change the original data object
        Data.copy <- rawData()
        if(customCDF){
          print ("Changing CDF before pre-processing")
          Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
        }
        
        if(uploadCDF){
          print ("Changing CDF before pre-processing")
          Data.copy <- uploadcdfenv(Data.copy, input$annot_file$datapath)
        }
        
        print ("Pre-processing is running")
        
        nGroups <- 1
        if(perGroup) {
          nGroups <- max(1,length(levels(experimentFactor)))
          if(nGroups==1) warning("normalization per group requested, 
                               but no groups indicated in data set")
        }
        
        #if per group normalization required, or a method selected that does not return an ExpressionSet object,
        #make a model of class ExpressionSet to paste real values in, use the relatively fast RMA method
        #note that binding of ExpressionSet objects is NOT possible
        if((nGroups>1)) { # || (normMeth=="MAS5")) {
          normData <- affy::rma(Data.copy)
          exprs(normData)[] <- NA
        }
        
        for(group in 1:nGroups) {
          if(nGroups==1) {
            Data.tmp <- Data.copy
          } else {
            Data.tmp <- Data.copy[,experimentFactor==(levels(experimentFactor)[group])]
          }
          switch(normMeth, 
                 "MAS5" = {
                   #doesn't work
                   normData.tmp <- mas5(Data.tmp) 
                 },
                 "GCRMA" = {
                   if(customCDF) {
                     #probe library needed, first try whether this has been intalled, otherwise do so
                     probeLibrary <- tolower(paste(Data@annotation,species,
                                                   CDFtype,"probe",sep=""))
                     loaded <- suppressWarnings(
                       try(eval(parse("",-1,paste("library(",probeLibrary,")", 
                                                  sep=""))),TRUE))
                     if(class(loaded)=="try-error") {
                       install.packages(probeLibrary, 
                                        repos="http://brainarray.mbni.med.umich.edu/bioc")
                     }
                   }
                   if(aType() == "PMMM") ntype = "fullmodel"
                   if(aType() == "PMonly") ntype = "affinities" # good results if most of the genes are not expressed
                   normData.tmp <- gcrma(Data.tmp, type=ntype, fast = FALSE)
                 },
                 "RMA" = {
                   normData.tmp <- affy::rma(Data.tmp)
                 },
                 "PLIER" = {
                   if(aType() == "PMMM") ntype = "together"
                   if(aType() == "PMonly") ntype = "pmonly"
                   normData.tmp <- justPlier(Data.tmp, normalize=TRUE, 
                                             norm.type = ntype)
                 }
          )
          if(nGroups==1) {
            normData <- normData.tmp
            if(normMeth=="MAS5") exprs(normData)<-log2(exprs(normData))
          } else {
            try(
              if(normMeth=="MAS5"){
                exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- log2(exprs(normData.tmp))
              }else{
                exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- exprs(normData.tmp)
              },TRUE)
          }
          rm(normData.tmp, Data.tmp)
        }
        
        rm(Data.copy) 
        
        normData
      },
      error=function(cond) {
        return(NULL)
      })
      
    }
    
    
    #Start from normalized data
    #..........................................................................#
    if ((rawornormalized() == "Normalized")){
      if (is.null(normData.saved())){
        if (database() == "GEO"){
          normData <- gset()[[paste0(accession(),"_series_matrix.txt.gz")]]
        }
        if (database() == "ArrayExpress"){
          normData <- gset()
        }
        if (database() == "Upload"){
          normData <- new(Class = "ExpressionSet", exprs = gxData())
        }
      }
      
      if (!is.null(normData.saved())){
        normData <- normData.saved()
      }
    }

    #..........................................................................#
    return(normData)
  })
  
  
  #Error message
  observeEvent(if (is.null(normData())){input$preprocessing},{
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "Error!",
      text = "Something went wrong. Try using a different annotation.",
      type = "error")
    
  })
  
  #Get expression matrix
  
  data.expr <- reactive({
    
    req(normData())
    data.expr <- exprs(normData())
    
    if (PSannotation() == "Custom annotations"){
     rownames(data.expr) <- str_remove(rownames(data.expr), "_at")
    }
    removeModal()
    
    return(data.expr)
    
  })
  
  
  #****************************************************************************#
  #Plots
  #****************************************************************************#
  
  #Normalized boxplot: empty
  output$boxplotNorm <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    png(outfile, width = 400, height = 300)
    NULL
    dev.off()
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  
  #Raw boxplot: empty
  output$boxplotRaw <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    # Generate the PNG
    png(outfile, width = 400, height = 300)
    NULL
    dev.off()
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  
  #Normalized heatmap: empty
  output$correlNormInt <- renderPlotly(NULL)
  
  #Normalized density plot: empty
  output$normhist <- renderPlotly(NULL)
  
  #Normalized density plot: empty
  output$ExprBoxplot <- renderPlotly(NULL)
  
  #Generate Plots
  observeEvent(if (length(data.expr()) > 0) TRUE, {
    
    #..........................................................................#
    # Normalized boxplot
    #..........................................................................#
    
    output$boxplotNorm<-renderImage({
      
      # Width of all the plots
      WIDTH <- 1000
      
      # Height of all the plots
      HEIGHT <- 1414
      
      # Point sizes for the plots
      POINTSIZE <- 24
      
      # Maximum number of arrays that can be computed
      MAXARRAY <- 41 
      
      
      myPalette <- colorsByFactor(experimentFactor())
      
      plotColors <- myPalette$plotColors
      
      # Legend colros
      legendColors <- myPalette$legendColors
      
      # Plot symbols
      plotSymbols <- 18-as.numeric(experimentFactor())
      
      # Legend symbols
      legendSymbols <- sort(plotSymbols, decreasing=TRUE)
      
      Type <- "Norm"
      tmain <- "Boxplot of normalized intensities"
      tmtext2 <- "Normalized log intensity\n\n\n"
      
      DataBoxplot<- tempfile(fileext='.png')
      png(file = DataBoxplot,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
      par(oma=c(17,0,0,0), cex.axis=1) 
      if (class(normData())[[1]] != "GeneFeatureSet"){
        suppressWarnings(boxplot(normData(), col=plotColors ,main=tmain, 
                                 axes=FALSE, pch = 20, cex=0.7))
      }
      if (class(normData())[[1]] == "GeneFeatureSet"){
        suppressWarnings(boxplot(normData(), target = "core", col=plotColors,
                                 main=tmain, axes=FALSE, pch = 20, cex=0.7))
      }
      if(length(levels(experimentFactor()))>1){ 
        legend("topright", levels(experimentFactor()),
               col=legendColors,fill=legendColors, cex = 0.7, bg = "white", 
               bty = "o")
      }
      if(length(sampleNames(normData()))<MAXARRAY){
        cexval <- 0.65
      }else{
        cexval <- 0.45
      }  
      axis(1,at=1:length(sampleNames(normData())),las=2,
           labels=sampleNames(normData()), cex.axis=cexval)        
      axis(2, cex.axis=0.7)  
      mtext(tmtext2, side=2, cex=0.8)  	   
      mtext("Distributions should be comparable between arrays\n", side=3, 
            font=1, cex=0.7)
      dev.off()
      
      list(src = DataBoxplot,width = WIDTH,height = HEIGHT,
           alt = "This is alternate text")
    },deleteFile = TRUE)
    
    
    
    
    #..........................................................................#
    # Raw boxplot
    #..........................................................................#
    
    
    output$boxplotRaw<-renderImage({
      if (!is.null(rawData())){
        rawData <- rawData()
        
        # Width of all the plots
        WIDTH <- 1000
        
        # Height of all the plots
        HEIGHT <- 1414
        
        # Point sizes for the plots
        POINTSIZE <- 24
        
        # Maximum number of arrays that can be computed
        MAXARRAY <- 41 
        
        
        myPalette <- colorsByFactor(experimentFactor())
        
        plotColors <- myPalette$plotColors
        
        # Legend colros
        legendColors <- myPalette$legendColors
        
        # Plot symbols
        plotSymbols <- 18-as.numeric(experimentFactor())
        
        # Legend symbols
        legendSymbols <- sort(plotSymbols, decreasing=TRUE)
        
        Type <- "Raw"
        tmain <- "Boxplot of raw intensities"
        tmtext2 <- "Raw log intensity\n\n\n"
        
        DataBoxplot<- tempfile(fileext='.png')
        png(file = DataBoxplot,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
        par(oma=c(17,0,0,0), cex.axis=1) 
        suppressWarnings(boxplot(rawData, col=plotColors ,main=tmain, 
                                 axes=FALSE, pch = 20, cex=0.7))
        if(length(levels(experimentFactor()))>1){ 
          legend("topright", levels(experimentFactor()),
                 col=legendColors,fill=legendColors, cex = 0.7, bg = "white", 
                 bty = "o")
        }
        if(length(sampleNames(rawData))<MAXARRAY){
          cexval <- 0.65
        }else{
          cexval <- 0.45
        }  
        axis(1,at=1:length(sampleNames(rawData)),las=2,
             labels=sampleNames(rawData), cex.axis=cexval)        
        axis(2, cex.axis=0.7)  
        mtext(tmtext2, side=2, cex=0.8)  	   
        mtext("Distributions should be comparable between arrays\n", side=3, 
              font=1, cex=0.7)
        dev.off()
        
        list(src = DataBoxplot,width = WIDTH,height = HEIGHT,
             alt = "This is alternate text")
      } else{
        outfile <- tempfile(fileext = '.png')
        
        # Generate the PNG
        png(outfile, width = 400, height = 300)
        NULL
        dev.off()
        
        # Return a list containing the filename
        list(src = outfile,
             contentType = 'image/png',
             width = 400,
             height = 300,
             alt = "This is alternate text")
      }
     
    },deleteFile = TRUE)
    
    
  
    
    #..........................................................................#
    # Interactive heatmap
    #..........................................................................#
    
    
    output$correlNormInt <- renderPlotly({
      myPalette <- colorsByFactor(experimentFactor())
      
      plotColors <- myPalette$plotColors
      
      # Legend colros
      legendColors <- myPalette$legendColors
      
      # Plot symbols
      plotSymbols <- 18-as.numeric(experimentFactor())
      
      # Legend symbols
      legendSymbols <- sort(plotSymbols, decreasing=TRUE)
      clusterOption1 <- input$clusteroption1
      clusterOption2 <- input$clusteroption2 
      normMeth <- input$normMeth
      
      
      Type <- "Norm"
      
      #note: for computing array correlation, euclidean would not make sense
      #only use euclidean distance to compute the similarity of the correlation 
      #vectors for the arrays
      COpt1 <- "pearson"
      if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
      crp <- cor(data.expr(), use="complete.obs", method=COpt1)
      
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
      
      my.hclust <- function(d) hclust(d, method=clusterOption2)

      
      names(legendColors) <- levels(experimentFactor())
      sidecolors <- data.frame(experimentFactor())
      colnames(sidecolors) <- "Grouping"
      
      
      
      if (input$heatmaptheme == "Default"){
       gradient = viridis(n = 256, alpha = 1, begin = 0, end = 1, 
                          option = "viridis")
      }
      
      if (input$heatmaptheme == "Yellow-red"){
        gradient = heat.colors(100)
      }
      
      if (input$heatmaptheme == "Dark"){
        gradient = RColorBrewer::brewer.pal(8, "Dark2")
      }
      
      heatmaply(crp, plot_method = "plotly", distfun = my.dist, 
                hclustfun = my.hclust, symm = TRUE, seriate = "mean", 
                titleX = FALSE, titleY = FALSE, key.title = NULL, 
                show_dendrogram = c(TRUE, FALSE), col_side_colors = sidecolors, 
                col_side_palette = legendColors, column_text_angle = 90,
                colors = gradient)
      
   })

    
    
    #..........................................................................#
    # Normalized Desity plot
    #..........................................................................#
    
    output$normhist <- renderPlotly({
      
      data.expr <- data.expr()
      meta1 <- meta1()
      
      #make dataframe with samples + intensities
      logData <- cbind.data.frame(as.vector(data.expr), 
                                  rep(colnames(data.expr), 
                                      each = nrow(data.expr)))
      
      colnames(logData) <- c("logInt", "sampleName")
      
      logData <- logData %>% inner_join(meta1, by = c("sampleName" = "names"))
      
      logData$Grouping <- factor(logData$Grouping, 
                                 levels = unique(meta1$Grouping))
      logData$Sample.ID <- factor(logData$Sample.ID, 
                                  levels = unique(meta1$Sample.ID))

      #Make colours
      
      myPalette <- colorsByFactor(experimentFactor())


      #Make density plot
      plot <- 
      ggplot(logData, aes(x = logInt, colour = Sample.ID, shape = Grouping)) +
        geom_density(size = 1) +
        scale_colour_manual(values = myPalette$plotColors) +
        labs(title = "Density plot of normalized intensities",
             subtitle = "Distributions should be comparable between arrays") +
        xlab("Normalized log intensity") +
        ylab("Density") +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.title = element_blank()) +
        guides(shape=guide_legend(title="Grouping", 
                                  override.aes = list(
                                    colour = myPalette$legendColors)), 
               colour = "none")
      
      
      return(ggplotly(plot))
      
    })
    
    
    #..........................................................................#
    # Expression matrix
    #..........................................................................#
    
    output$ExpressionMatrix <- DT::renderDataTable(
      data.expr(), selection = list(mode = "single", selected = 1)
    )
    
    output$ExprBoxplot <- renderPlotly({
      
      myPalette <- colorsByFactor(experimentFactor())
      legendColors <- myPalette$legendColors
      
      selected_row <- as.data.frame(data.expr())[input$ExpressionMatrix_rows_selected,]
      
      plotExpr <- selected_row %>%
        gather(key = "Sample", value = "logExpr") %>%
        inner_join(meta1(), by = c("Sample" = "names"))
      
      plotExpr %>%
        plot_ly(x = ~Grouping,y = ~as.numeric(logExpr), 
                color = ~Grouping, colors = legendColors, type = "box") %>%
        layout(xaxis = list(title = " "), 
               yaxis = list(title = 'logExpr'), 
               legend = list(title=list(text='<b> Species of Iris </b>')),
               showlegend = FALSE)
        
    })
    
    
  })
  
  

  
  
  #****************************************************************************#
  #After normalization
  #****************************************************************************#
  
  
  observeEvent(if (length(data.expr()) > 0) TRUE, {
    
  #Download expression matrix
  output$downloadexpr <- downloadHandler(
    filename = "ExpressionMatrix.txt",
    content = function(file){
    write.table(data.expr(), file, sep = "\t")
  })
    
  #Remove loading page
  removeModal()
    
  #Success message
  if (rawornormalized() == "Raw"){
        sendSweetAlert(
          session = session,
          title = "Success!!",
          text = "Data successfully pre-processed! 
          Please wait for the QC plots to be rendered.",
          type = "success"
        )
    }
  })
  
  #Error message
  observeEvent(input$preprocessing,{
    if (rawornormalized() == "Normalized"){
      sendSweetAlert(
        session = session,
        title = "Error!!",
        text = "Pre-processing is not possible because you already started with 
        normalized data!",
        type = "error"
      )
    }
  })
  
  #Go to previous tab
  observeEvent(input$ann.back, {
    if (length(meta()) > 0) {
      updateNavbarPage(session, "navbar",
                       selected = "panel2")
    }

  })
  
  #Go to next tab
  observeEvent((input$ann.proceed), {
    if (length(data.expr()) > 0){
      updateNavbarPage(session, "navbar",
                       selected = "panel4")
    }
    
  })
  
  #Show next tab
  observeEvent(input$ann.proceed, {
    if (length(data.expr()) > 0){
      showTab("navbar", target = "panel4")
    }
  })
  
  
  #****************************************************************************#
  # Save normalized data
  #****************************************************************************#
  
  observeEvent(input$save, {
    
    if (length(data.expr()) > 0){
      id <- str_remove_all(Sys.time(), "-")
      id <- str_remove_all(id, " ")
      id <- str_remove_all(id, ":")
      id <- paste0("aa", id)
      
      id.normData <- paste0("normData", id)
      id.meta1 <- paste0("meta1", id)
      id.info <- paste0("info", id)
      
      saveRDS(normData(), id.normData)
      saveRDS(meta1(), id.meta1)
      
      info <- list(PSannotation(), CDFtype(), species())
      saveRDS(info, id.info)
      
      sendSweetAlert(
        session = session,
        title = "Data successfully saved!!",
        text = paste0("Please use the following ID to access your data: ", id),
        type = "success"
      )
    }
    
    if (length(normData()) == 0){
      sendSweetAlert(
        session = session,
        title = "Error!!",
        text = "Data needs to be pre-processed first!!",
        type = "Error"
      )
    }
    
  })

  
  
  ##############################################################################
  # PCA
  ##############################################################################
  
  #Calculate PCs
  data.PC <- eventReactive(input$ann.proceed, {
    
    data.expr <- data.expr()
    
    #Remove rows with constant values
    rowVar <- apply(data.expr, 1, var)
    data.expr <- data.expr[(rowVar != 0),]
    
    #Perform PCA
    prcomp(t(data.expr),scale.=TRUE)
  })
  
  #Selection z-axis for 3D plot
  output$uizpca <- renderUI({
    if (input$plot3d == TRUE) {
      selectInput(inputId = "zpca", 
                  label = "z-axis",
                  choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", 
                              "PC7", "PC8"),
                  selected = "PC3")
    }
  })
  
  #PC: x-axis
  pc.x <- reactive({
    
    switch(input$xpca,
           "PC1" = 1,
           "PC2" = 2,
           "PC3" = 3,
           "PC4" = 4,
           "PC5" = 5,
           "PC6" = 6,
           "PC7" = 7,
           "PC8" = 8)
    
  })
  
  #PC: y-axis
  pc.y <- reactive({
    switch(input$ypca,
           "PC1" = 1,
           "PC2" = 2,
           "PC3" = 3,
           "PC4" = 4, 
           "PC5" = 5,
           "PC6" = 6,
           "PC7" = 7,
           "PC8" = 8)
    
  })
  
  #PC: z-axis
  pc.z <- reactive({
    
    pc.z <- NULL
    
    if(input$plot3d == TRUE){
      if (length(input$zpca) > 0) {
        pc.z <- switch(input$zpca,
                       "PC1" = 1,
                       "PC2" = 2,
                       "PC3" = 3,
                       "PC4" = 4, 
                       "PC5" = 5,
                       "PC6" = 6,
                       "PC7" = 7,
                       "PC8" = 8)
        
      }
    }
    
    return(pc.z)

  })
  
  #Render 2D PCA plot
  output$pca <- renderPlotly({
    if (length(input$plot3d) > 0){
      if (input$plot3d == FALSE){
        pca.plot(data.PC(), meta1(), pc.x(), pc.y())
      }
    }
  })
  
  
  output$uipca = renderUI({
    if (length(input$plot3d) > 0) {
      if (input$plot3d == FALSE){
        plotlyOutput("pca") %>% withSpinner(color="#0dc5c1")
      }
    }
  })
  
  #Render 3D PCA plot
  output$pca3d <- renderPlotly({
    if (length(input$plot3d) > 0){
      if (input$plot3d == TRUE){
        pca3d(data.PC(), meta1(), x = pc.x(), y = pc.y(), z = pc.z())
      }
    }
  })
  
  output$uipca3d = renderUI({
    if (length(input$plot3d) > 0) {
      if (input$plot3d == TRUE){
        plotlyOutput("pca3d") %>% withSpinner(color="#0dc5c1")
      }
    }
  })
  
  
  #Render explained variances plot
  output$variances <- renderPlotly({
    PCvariances(data.PC(), x = pc.x(), y = pc.y(), z = pc.z())
  })
  
  observeEvent(input$pca.ok, {
    sendSweetAlert(
      session = session,
      title = "Success!!",
      text = "Now you can perform various statistical analyses!",
      type = "success")
    
  })
  
  #Go to previous tab
  observeEvent(input$pca.back, {
    updateNavbarPage(session, "navbar",
                     selected = "panel3")
    
    
  })
  
  #Go to next tab
  observeEvent(input$pca.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel5")
    
  })
  
  #Show next tab
  observeEvent(input$pca.ok, {
    showTab("navbar", target = "panel5")
  })
  
  ##############################################################################
  # Statistical analysis
  ##############################################################################
  
  
  #get statistics
  
  top.table <- reactive({
    
    top.table <- diff_expr(data.expr(), meta1(), 
                           comparisons = get_contrasts(meta1()))
    
    return(top.table)
  })
  
  
  ##############################################################################
  # Top table
  ##############################################################################
  
  
  #Select comparison
  output$uitoptablecomp <- renderUI({
    pickerInput(
      inputId = "toptablecomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #Get the two possible comparison directions
  bothdirectionstop <- reactive({
    req(input$toptablecomp)
    
    contrast <- input$toptablecomp
    
    switch <- strsplit(contrast, split= " - ")
    
    contrast2 <- paste0(switch[[1]][2], " - ", switch[[1]][1])
    
    bothdirections <- c(contrast, contrast2)
    
    return(bothdirections)
  })
  
  #Ask for direction of comparison
  output$uitoptabledir <- renderUI({
    pickerInput(
      inputId = "toptabledir",
      label = "Direction of comparison",
      choices = bothdirectionstop(),
      options = list(
        style = "btn-info")
    )
  })

  #Comparison of interest
  toptablecomp <- reactive(
    input$toptablecomp
  )
  
  #Direction of interest
  toptabledir <- reactive(
    input$toptabledir
  )
  
  #Top table
  toptable <- reactive(
    if (length(toptabledir()) > 0){
      toptable <- top.table()[[toptablecomp()]]
      
      if (toptablecomp() != toptabledir()){
        toptable$logFC <- -1 * (toptable$logFC)
        toptable$t <- -1 * (toptable$t)
      }
      
      return(toptable)
    }
  )
  
  #print table
  output$x1 <- DT::renderDataTable({
   
    req(toptable())
    
    toptable <- toptable()
    
    if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
      if(CDFtype() == "ENTREZG"){
        
        toptable$Probeset.ID = paste0(
          '<a ',
          'href=',
          paste(
            "https://www.ncbi.nlm.nih.gov/gene/",
            toptable$Probeset.ID,
            sep = ''
          ),
          '>',
          toptable$Probeset.ID,
          '</a>'
        )
        
      }
    }
    
    return(toptable)
    
  }, selection = list(mode = "single", selected = 1), rownames = FALSE, escape = FALSE)
  
  #Download table
  output$downloadfinal <- downloadHandler(
    filename = "TopTable.txt",
    content = function(file){
      write.table(toptable(), file, sep = "\t", row.names = FALSE, 
                  quote = FALSE)
    }
  )
  
  #Boxplot
  output$topboxplot <- renderPlotly({
    if (length(toptabledir()) > 0){
      selected_row = toptable()[input$x1_rows_selected,]
      
      comparison <- strsplit(toptablecomp(), " - ")[[1]]
      
      meta2 <- meta1()
      meta2$Grouping <- make.names(meta2$Grouping)
      
      
      samples <- meta2 %>%
        filter(Grouping == comparison[1] | Grouping == comparison[2])
      
      data.expr1 <- data.expr()[,samples$names]
      data.expr1 <- data.expr1[selected_row$Probeset.ID,]
      
      if (length(as.vector(data.expr1)) == length(samples$Grouping)){
        if (toptablecomp() == toptabledir()){
          dat <- data.frame(logExpr = as.vector(data.expr1),
                            group = factor(samples$Grouping))
        }
        
        if (toptablecomp() != toptabledir()){
          dat <- data.frame(logExpr = as.vector(data.expr1),
                            group = factor(samples$Grouping, 
                                           levels = 
                                             c(unique(samples$Grouping)[2], 
                                               unique(samples$Grouping)[1])))
        }
        
        dat %>%
          plot_ly() %>% 
          add_trace(x = ~as.numeric(group),y = ~logExpr, color = ~group, 
                    colors = c("red", "blue"), type = "box", 
                    hoverinfo = 'name+y') %>%
          add_markers(x = ~jitter(as.numeric(group)), y = ~logExpr, 
                      color = ~group,
                      marker = list(size = 6),
                      hoverinfo = "text",
                      text = ~paste0("Group: ",group,
                                     "<br>logExpr: ",logExpr),
                      showlegend = FALSE) %>% 
          layout(legend = list(orientation = "h",
                               x =0.5, xanchor = "center",
                               y = 1, yanchor = "bottom"
          ),
          xaxis = list(title = "Group",
                       showticklabels = FALSE))
      }
      
    }
  })
  
  ##############################################################################
  # P-value analysis
  ##############################################################################
  
  
  #Select comparison
  output$uipcomp <- renderUI({
    pickerInput(
      inputId = "pcomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #P-Value histogram
  output$phist <- renderPlotly({
    req(input$pcomp)
    phist <- 
      ggplot(data = top.table()[[input$pcomp]], aes(x = P.Value)) +
      geom_histogram(bins = 100, colour = "#696969", fill = "#d3d3d3") +
      labs(title = "P-value histogram") +
      xlab("p-values") +
      ylab("Count") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    return(ggplotly(phist))
    
  })
  
 
  ##############################################################################
  # FC analysis
  ##############################################################################
  
  
  #Select comparison
  output$uilogfccomp <- renderUI({
    pickerInput(
      inputId = "logfccomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #logFC histogram
  output$logfchist <- renderPlotly({
    req(input$logfccomp)
    logfchist <- 
      ggplot(data = top.table()[[input$logfccomp]], aes(x = logFC)) +
      geom_histogram(bins = 100, colour = "#696969", fill = "#d3d3d3") +
      labs(title = "logFC histogram") +
      xlab("logFC") +
      ylab("Count") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    return(ggplotly(logfchist))
    
  })
  
  
  ##############################################################################
  # volcano plot
  ##############################################################################
  
  #Select comparison
  output$uivolcanocomp <- renderUI({
    pickerInput(
      inputId = "volcanocomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  #Get directions of comparison
  bothdirectionsvol <- reactive({
    req(input$volcanocomp)
    
    contrast <- input$volcanocomp
    
    switch <- strsplit(contrast, split= " - ")
    
    contrast2 <- paste0(switch[[1]][2], " - ", switch[[1]][1])
    
    bothdirections <- c(contrast, contrast2)
    
    return(bothdirections)
  })
  

  #Select direction of comparison
  output$uivolcanodir <- renderUI({
    pickerInput(
      inputId = "volcanodir",
      label = "Direction of comparison",
      choices = bothdirectionsvol(),
      options = list(
        style = "btn-info")
    )
  })
  
  
  #Selected comparison
  volcanocomp <- reactive({
    input$volcanocomp
  })
  
  #Raw or adjusted P-value
  p.choice <- reactive({
    input$raworadj
  })
  
  #p-value threshold
  p.threshold <- reactive({
    input$pthreshold
  })
  
  #logFC threshold
  logFC.threshold <- reactive({
    input$logfcthreshold
  })
  
  #Direction of comparison
  volcanodir <- reactive({
    input$volcanodir
  })
  
  
  #Make volcano plot
  volcano <- eventReactive(input$volcano.ok, {
    
    if (length(volcanodir()) > 0){
      top.table <- top.table()[[volcanocomp()]]
      
      if (volcanocomp() != volcanodir()) {
        
        top.table$logFC <- -1 * (top.table$logFC)
        top.table$t <- -1 * (top.table$t)
      }
      
      plotvolcano(top.table, 
                  p = p.choice(), 
                  p.threshold = p.threshold(),
                  logFC.threshold = logFC.threshold())
      
      
    }
    
  }, ignoreNULL = FALSE)
  
  
  output$volcano <- renderPlotly({
    volcano()
  })
  
  
  
  #Top table of significant genes
  voltable <- eventReactive(input$volcano.ok, {
    
    if (length(volcanodir()) > 0){
      voldata <- top.table()[[volcanocomp()]]
      
      if (p.choice() == "raw"){
        voltable <- voldata %>%
          filter(P.Value <= p.threshold()) %>%
          filter(abs(logFC) >= logFC.threshold())
        
        
      }
      
      if (p.choice() == "adj"){
        voltable <- voldata %>%
          filter(adj.P.Val <= p.threshold()) %>%
          filter(abs(logFC) >= logFC.threshold())
        
      }
      
      if (volcanocomp() != volcanodir()) {
        
        voltable$logFC <- -1 * (voltable$logFC)
      }
      
      return(voltable)
      
    }
    
  }, ignoreNULL = FALSE)
  
    
  #Render top table with significant genes only
  output$voltable <- renderDataTable({
    voltable()
  })
  

  #Download table with significant genes
  output$downloadvol <- downloadHandler(
    filename = "VolcanoTopTable.txt",
    content = function(file){
      write.table(voltable(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  ##############################################################################
  # venn diagram
  ##############################################################################
  
  
  #Select comparisons (between 2 and 4)
  output$uicontrast <- renderUI({
    pickerInput(
      inputId = "contrastselection",
      label = "Select between 2 and 4 comparisons", 
      choices = get_contrasts(meta1()),
      multiple = TRUE,
      options = list(
        style = "btn-primary")
    )
  })
  
  
  #Plot venn diagram
  output$venndiagram <- renderPlot({
    req(input$contrastselection)
    
    if ((length(input$contrastselection) < 5) & (length(input$contrastselection) > 1)){
      contrasts <- input$contrastselection
      
      toptableselection <- list()
      for (j in contrasts) {
        toptableselection[[j]] <- top.table()[[j]]
      }
      
      genelist <- NULL
      for (i in names(toptableselection)) {
        if (input$raworadjvenn == "raw") {
          genes <- toptableselection[[i]] %>% 
            filter(P.Value <= input$pthresholdvenn) %>%
            filter(logFC > input$logfcthresholdvenn | logFC < (-1 * input$logfcthresholdvenn))
        }
        
        if (input$raworadjvenn == "adj") {
          genes <- toptableselection[[i]] %>% 
            filter(adj.P.Val <= input$pthresholdvenn) %>%
            filter(logFC > input$logfcthresholdvenn | logFC < (-1 * input$logfcthresholdvenn))
        }
        
        
        genelist[[i]] <- genes$Probeset.ID
      }
      
      names(genelist) <- c("A", "B", "C", "D")[1:length(contrasts)]
      
      ggvenn(genelist, set_name_size = 4, text_size = 3, stroke_size = 0.7)
      
    }
    
  })
  
  
  #Plot legend of venn diagram
  output$legendvenn <- renderText({
    if (length(input$contrastselection) == 2) {
      text <- paste0("<p> <b> A: </b> ", 
                     input$contrastselection[1], 
                     "<br /> <p> <b> B: </b> ", 
                     input$contrastselection[2])
    }
    
    if (length(input$contrastselection) == 3) {
      text <- paste0("<p> <b> A: </b> ", 
                     input$contrastselection[1], 
                     "<br /> <p> <b> B: </b> ", 
                     input$contrastselection[2], 
                     "<br /> <p> <b> C: </b> ", 
                     input$contrastselection[3])
    }
    
    if (length(input$contrastselection) == 4) {
      text <- paste0("<p> <b> A: </b> ", 
                     input$contrastselection[1], 
                     "<br /> <p> <b> B: </b> ", 
                     input$contrastselection[2], 
                     "<br /> <p> <b> C: </b> ", 
                     input$contrastselection[3], 
                     "<br /> <p> <b> D: </b> ", 
                     input$contrastselection[4])
    }
   
    if ((length(input$contrastselection) > 4)){
      text <- "More than four comparisons selected. Please select between two and four comparisons."
    }
    
    if (length(input$contrastselection) < 2){
      text <- "Less than two comparisons selected. Please select between two and four comparisons."
    }
    
    return(text)
  })
  
  ##############################################################################
  #Heatmap
  ##############################################################################
  
  #Select comparison
  output$uiheatmapcomp <- renderUI({
    pickerInput(
      inputId = "heatmapcomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  #Plot heatmap
  output$topfeatureheatmap <- renderPlotly({
    
    req(input$heatmapcomp)
    
    comparison <- input$heatmapcomp
    comparison <- strsplit(comparison, " - ")[[1]]
    
    meta2 <- meta1()
    meta2$Grouping <- make.names(meta2$Grouping)
    
    samples <- meta2 %>%
      filter(Grouping == comparison[1] | Grouping == comparison[2])
    
    
    topfeatures <- head(top.table()[[input$heatmapcomp]], input$topfeatures)
    data.expr1 <- data.expr()[topfeatures$Probeset.ID,]
    data.expr1 <- data.expr1[,samples$names]
    
    
    clusterOption1 <- input$clusteroption3
    clusterOption2 <- input$clusteroption4 
    
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
    
    my.hclust <- function(d) hclust(d, method=clusterOption2)
    
    experimentFactor <- factor(samples$Grouping)
    myPalette <- colorsByFactor(experimentFactor)
    legendColors <- myPalette$legendColors
    names(legendColors) <- levels(experimentFactor)
    sidecolors <- data.frame(experimentFactor)
    colnames(sidecolors) <- "Grouping"
    
    
    if (input$heatmaptheme1 == "Default"){
      gradient = viridis(n = 256, alpha = 1, begin = 0, end = 1, 
                         option = "viridis")
    }
    
    if (input$heatmaptheme1 == "Yellow-red"){
      gradient = heat.colors(100)
    }
    
    if (input$heatmaptheme1 == "Dark"){
      gradient = RColorBrewer::brewer.pal(8, "Dark2")
    }
    
    
    heatmaply(data.expr1, plot_method = "plotly", distfun = my.dist, hclustfun = my.hclust, seriate = "mean", 
              titleX = TRUE, titleY = TRUE, xlab = "Samples", ylab = "Features", main = paste0("Heatmap of the top ", input$topfeatures, " features"), 
              key.title = "Expression", column_text_angle = 90, col_side_colors = sidecolors, col_side_palette = legendColors, colors = gradient)
  })
  
  ##############################################################################
  # GO analysis
  ##############################################################################
  
  
  #****************************************************************************#
  #   topGO
  #****************************************************************************#
  
  #Select comparison
  output$uigoacomp <- renderUI({
    pickerInput(
      inputId = "goacomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  #Make topGO data object
  sampleGOdata <- eventReactive(input$goa.ok, {
    
    #Show loading message
    showModal(modalDialog(title = h4("GO enrichment analysis is being performed....", 
                                     align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", 
                             align = "center"),
                          br()))
    
    #Select top table
    top.table <- top.table()[[input$goacomp]]
    
    #Get significant genes
    if(input$raworadjgoa == "raw"){
      geneList <- ifelse((top.table$P.Value < input$pthresholdgoa) 
                         & (abs(top.table$logFC) > input$logFCthresholdgoa), 1, 0)
    }
    
    if(input$raworadjgoa == "adj"){
      geneList <- ifelse((top.table$adj.P.Val < input$pthresholdgoa)  
                         & (abs(top.table$logFC) > input$logFCthresholdgoa), 1, 0)
    }
    
    names(geneList) <- top.table$Probeset.ID
    
    
    #Make topGOdata object
    sampleGOdata <- NULL
    
    #No custom annotation
    if ((PSannotation() == "No annotations") & (rawornormalized() == "Raw")){

      if (annotation(normData()) == "hugene10stv1"){
        affyLib <- "hugene10stprobeset.db"
        
      } else {
        affyLib <- paste(annotation, "db", sep = ".")
        
      }
      
      #install annotation package from bioconductor
      if (!requireNamespace(affyLib, quietly = TRUE))
        BiocManager::install(affyLib, ask = FALSE)
      

      #If installation successful
      if(require(as.character(affyLib), character.only = TRUE)){
        
        #Make topGOdata object
        sampleGOdata <- new("topGOdata",
                            description = "Simple session", ontology = input$ontology,
                            allGenes = geneList, geneSelectionFun = function(x)(x == 1),
                            nodeSize = 10,
                            annot = annFUN.db, affyLib = affyLib)
      }
    } 
    
    
    #Custom annotation
    if ((PSannotation() == "Custom annotations") & (rawornormalized() == "Raw")){
      if((CDFtype() == "ENTREZG") & (species() == "Homo sapiens")){
        sampleGOdata <- new("topGOdata",
                            description = "Simple session", ontology = input$ontology,
                            allGenes = geneList, geneSelectionFun = function(x)(x == 1),
                            nodeSize = 10,
                            annot = annFUN.org,
                            mapping = "org.Hs.eg.db")
      } 
    } 
    
    return(sampleGOdata)
  })
  
  
  #GO enrichment
  resultFisher <- eventReactive(input$goa.ok, {
    if (!is.null(sampleGOdata())){
      resultFisher <- runTest(sampleGOdata(), algorithm = "classic", statistic = "fisher")
      return(resultFisher)
    }
  })
  
  
  #Outputs
  output$goa <- renderDataTable(NULL)
  output$errorgoa <- renderText(NULL)
  
  observeEvent(input$goa.ok, {
    
    #error message
    if (is.null(sampleGOdata())){
      removeModal()
      output$errorgoa <- renderText({
        "GO enrichment analysis is currently not possible using the selected 
        probeset annotation."
      })
    }
    
    if (!is.null(sampleGOdata())){
      removeModal()
      
      #Table output
      output$goa <- renderDataTable({
        
        allRes <- GenTable(sampleGOdata(), classicFisher = resultFisher(), 
                           orderBy = "classicFisher", ranksOf = "classicFisher", 
                           topNodes = 30)
        
        return(allRes)
        
      }, options = list(lengthMenu = c(5, 10, 20, 30), pageLength = 5))
      
      
      #Plot output
      output$GOplot <-renderImage({
        # Width of all the plots
        WIDTH <- 1000
        
        # Height of all the plots
        HEIGHT <- 1414
        
        # Point sizes for the plots
        POINTSIZE <- 24
        
        GOplot <- tempfile(fileext='.png')
        png(file = GOplot,width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)  
        par(cex = 0.3) 
        showSigOfNodes(sampleGOdata(), score(resultFisher()), firstSigNodes = 3, useInfo = 'all')
        title(main = list("Subgraph of the top 3 GO terms", cex = 3, font = 2), line = -1)
        dev.off()
        
        list(src = GOplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
        
        
      },deleteFile = TRUE)

    }

  })
  
  
  
  
  #****************************************************************************#
  #   clusterProfiler
  #****************************************************************************#
  
  #Select comparison
  output$uigoacomp1 <- renderUI({
    pickerInput(
      inputId = "goacomp1",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #Select organism
  output$uigoaOrg <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "goaOrg",
        label = "Organism",
        choices = "hsapiens_gene_ensembl",
        options = list(
          style = "btn-danger")
      )
    }
  })
  
  #Select chiptype
  output$uigoaChip <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "goaChip",
        label = "Chip type",
        choices = c("affy_hc_g110", "affy_hg_focus", "affy_hg_u133a", "affy_hg_u133a_2", 
                    "affy_hg_u133b", "affy_hg_u133_plus_2", "affy_hg_u95a", "affy_hg_u95av2",
                    "affy_hg_u95b", "affy_hg_u95c", "affy_hg_u95d", "affy_hg_u95d", "affy_hg_u95e",
                    "affy_hta_2_0", "affy_huex_1_0_st_v2", "affy_hugenefl", "affy_hugene_1_0_st_v1",
                    "affy_hugene_2_0_st_v1", "affy_hugene_2_1_st_v1", "affy_primeview", "affy_u133_x3p"),
        selected = get_chiptype(annotation(normData())),
        options = list(
          style = "btn-warning")
      )
    }
  })
  
  output$goa1 <- renderDataTable(NULL)
  output$GOplot1 <- renderPlot(NULL)
  output$errorgoa1 <- renderText(NULL)
  

  
  #Get gene list
  observeEvent(input$goa.ok1, {
    
    #Show loading message
    showModal(modalDialog(title = h4("GO analysis is being performed....", 
                                     align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", 
                             align = "center"),
                          br()))
    
    #Get gene list
    genelist <- eventReactive(input$goa.ok1,{
      
      genelist <- NULL
    
      
      #Custom annoation
      if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
        top.table <- top.table()[[input$goacomp1]]
        
        
        #ORA
        if(input$oraGO == "ora"){
          bgenelist <- top.table$Probeset.ID
          
          if(input$raworadjgoa1 == "raw"){
            siggenelist <- top.table %>%
              filter(P.Value < input$pthresholdgoa1) %>%
              filter(abs(logFC) > input$logFCthresholdgoa1)
          }
          
          if(input$raworadjgoa1 == "adj"){
            siggenelist <- top.table %>%
              filter(adj.P.Val < input$pthresholdgoa1) %>%
              filter(abs(logFC) > input$logFCthresholdgoa1)
          }
          
          siggenelist <- siggenelist$Probeset.ID
          
          
          genelist <- list(bgenelist, siggenelist)
          
        }

        
        #GSEA
        if(input$oraGO == "gsea"){
          if (input$rankingGO == "logFC") {
            genelist <- top.table$logFC
          }
          if (input$rankingGO == "-log(P-value)") {
            genelist <- -log10(top.table$P.Value)
          }
          
          # name the vector
          names(genelist) <- top.table$Probeset.ID
          
          # omit any NA values 
          genelist <-na.omit(genelist)
          
          # sort the list in decreasing order (required for clusterProfiler)
          genelist <- sort(genelist, decreasing = TRUE)
        }
        
      }
      
      #No annotations
      if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
        
        for (i in 1:10) {
          top.table <- tryCatch({
            top.table1 <- top.table()[[input$goacomp1]]
            
            ensembl <- useMart("ensembl")
            ensembl <- useDataset(input$goaOrg,mart=ensembl)
            
            annotation1 <- getBM(attributes = c(input$goaChip, 'ensembl_gene_id'), 
                                 filters = input$goaChip, 
                                 values = top.table1$Probeset.ID, 
                                 mart = ensembl)
            
            colnames(annotation1) <- c("Probeset.ID", "Gene.ID")
            annotation1$Probeset.ID <- as.character(annotation1$Probeset.ID)
            
            top.table <- top.table1 %>%
              inner_join(annotation1, by = c("Probeset.ID" = "Probeset.ID"))
            
            top.table
            
          },
          error=function(cond) {
            return(NULL)
          })
          if (!is.null(top.table)){
            break
          } 
        }
        
        
        if (!is.null(top.table)){
          
          #ORA
          if(input$oraGO == "ora"){
            bgenelist <- top.table$Gene.ID
            
            if(input$raworadjgoa1 == "raw"){
              siggenelist <- top.table %>%
                filter(P.Value < input$pthresholdgoa1) %>%
                filter(abs(logFC) > input$logFCthresholdgoa1)
            }
            
            if(input$raworadjgoa1 == "adj"){
              siggenelist <- top.table %>%
                filter(adj.P.Val < input$pthresholdgoa1) %>%
                filter(abs(logFC) > input$logFCthresholdgoa1)
            }
            
            siggenelist <- siggenelist$Gene.ID
            
            
            genelist <- list(bgenelist, siggenelist)
            
          }
          
          
          #GSEA
          if(input$oraGO == "gsea"){
            if (input$rankingGO == "logFC") {
              genelist <- top.table$logFC
            }
            if (input$rankingGO == "-log(P-value)") {
              genelist <- -log10(top.table$P.Value)
            }
            
            # name the vector
            names(genelist) <- top.table$Gene.ID
            
            # omit any NA values 
            genelist <-na.omit(genelist)
            
            # sort the list in decreasing order (required for clusterProfiler)
            genelist <- sort(genelist, decreasing = TRUE)
          }
          
        }
        
        
      }
      
      return(genelist)
    })
    
    
    #Get GO results
    GOresults <- eventReactive(input$goa.ok1,{
      GOresults <- NULL

      if (!is.null((genelist()))){
        if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
          if ((CDFtype() == "ENTREZG") & (species() == "Homo sapiens")){
            
            if(input$oraGO == "ora"){
              GOresults <- enrichGO(
                gene = genelist()[[2]], 
                universe = genelist()[[1]], 
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = input$ontology1,
                pAdjustMethod = "fdr",
                pvalueCutoff = 1,
                qvalueCutoff = 1)
            }
            
            if(input$oraGO == "gsea"){
              GOresults <- gseGO(
                geneList = genelist(),
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = input$ontology1,
                pAdjustMethod = "fdr",
                pvalueCutoff = 1)
            }
            
          }
        }
        
        if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
          if(input$oraGO == "ora"){
            GOresults <- enrichGO(
              gene = genelist()[[2]], 
              universe = genelist()[[1]], 
              keyType = "ENSEMBL",
              OrgDb = org.Hs.eg.db,
              ont = input$ontology1,
              pAdjustMethod = "fdr",
              pvalueCutoff = 1,
              qvalueCutoff = 1)
          }
          
          if(input$oraGO == "gsea"){
            GOresults <- gseGO(
              geneList = genelist(),
              keyType = "ENSEMBL",
              OrgDb = org.Hs.eg.db,
              ont = input$ontology1,
              pAdjustMethod = "fdr",
              pvalueCutoff = 1)
          }

        }
        
      }
      return(GOresults)
    })
    
    
    if (is.null(GOresults())){
      removeModal()
      
      #Error
      output$errorgoa1 <- renderText({
        "GO analysis is unfortunately not possible using the current probeset annotation."
      })
    }
    
    if (!is.null(GOresults())){
      removeModal()
      
      #Table
      output$goa1 <- renderDataTable({
        
        allRes <- GOresults()@result
        rownames(allRes) <- NULL
        
        return(allRes)
        
      }, options = list(lengthMenu = c(5, 10, 20, 30), pageLength = 5),
      selection = list(mode = "single", selected = 1))
      
      
      #Dotplot
      output$GOplot1 <- renderPlot({
        
        dotplot(GOresults(), showCategory = 20, orderBy = "x")

      })
      
      #GSEA plot
      output$GOplotgsea1 <- renderPlot({
        req(input$oraGO)
        if(input$oraGO == "gsea"){
          gseaplot2(GOresults(), geneSetID = input$goa1_rows_selected, 
                    title = GOresults()$Description[input$goa1_rows_selected])
        }
      })
    }
    
  })
    

  
 
  ##############################################################################
  #Pathway enrichment analysis
  ##############################################################################
  
  OraGSEA_WP <- eventReactive(input$wp.ok, {
    input$oraPathway
  })
    
  OraGSEA_KEGG <- eventReactive(input$kegg.ok, {
    input$oraPathway
  })
  
  #****************************************************************************#
  #   KEGG
  #****************************************************************************#
  
  #Select comparison
  output$uikeggcomp <- renderUI({
    pickerInput(
      inputId = "keggcomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #Select organism
  output$uikeggOrg <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "keggOrg",
        label = "Organism",
        choices = "hsapiens_gene_ensembl",
        options = list(
          style = "btn-danger")
      )
    }
  })
  
  #Select chiptype
  output$uikeggChip <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "keggChip",
        label = "Chip type",
        choices = c("affy_hc_g110", "affy_hg_focus", "affy_hg_u133a", "affy_hg_u133a_2", 
                    "affy_hg_u133b", "affy_hg_u133_plus_2", "affy_hg_u95a", "affy_hg_u95av2",
                    "affy_hg_u95b", "affy_hg_u95c", "affy_hg_u95d", "affy_hg_u95d", "affy_hg_u95e",
                    "affy_hta_2_0", "affy_huex_1_0_st_v2", "affy_hugenefl", "affy_hugene_1_0_st_v1",
                    "affy_hugene_2_0_st_v1", "affy_hugene_2_1_st_v1", "affy_primeview", "affy_u133_x3p"),
        selected = get_chiptype(annotation(normData())),
        options = list(
          style = "btn-warning")
      )
    }
  })
  
  output$kegg <- renderDataTable(NULL)
  output$KEGGplot <- renderPlotly(NULL)
  output$errorkegg <- renderText(NULL)
  
  
  
  observeEvent(input$kegg.ok, {
    
    #Show loading message
    showModal(modalDialog(title = h4("Pathway analysis is being performed....", 
                                     align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", 
                             align = "center"),
                          br()))
    
    kegggenelist <- eventReactive(input$kegg.ok,{
      
      genelist <- NULL
      
      
      #Custom annotation
      if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
        if(CDFtype() == "ENTREZG"){
          top.table <- top.table()[[input$keggcomp]]
        
          
          #ORA
          if(input$oraPathway == "ora"){
            bgenelist <- top.table$Probeset.ID
            
            if(input$raworadjkegg == "raw"){
              siggenelist <- top.table %>%
                filter(P.Value < input$pthresholdkegg) %>%
                filter(abs(logFC) > input$logFCthresholdkegg)
            }
            
            if(input$raworadjkegg == "adj"){
              siggenelist <- top.table %>%
                filter(adj.P.Val < input$pthresholdkegg) %>%
                filter(abs(logFC) > input$logFCthresholdkegg)
            }
            
            siggenelist <- siggenelist$Probeset.ID
            
            
            genelist <- list(bgenelist, siggenelist)
          }

          
          #GSEA
          if(input$oraPathway == "gsea"){
            if (input$rankingKEGG == "logFC") {
              genelist <- top.table$logFC
            }
            if (input$rankingKEGG == "-log(P-value)") {
              genelist <- -log10(top.table$P.Value)
            }
            
            
            # name the vector
            names(genelist) <- top.table$Probeset.ID
            
            # omit any NA values 
            genelist <-na.omit(genelist)
            
            # sort the list in decreasing order (required for clusterProfiler)
            genelist <- sort(genelist, decreasing = TRUE)
          }
          
        }
      }
      
      if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
        
        
          top.table <- tryCatch({
            top.table1 <- top.table()[[input$keggcomp]]
            
            ensembl <- useMart("ensembl")
            ensembl <- useDataset(input$keggOrg,mart=ensembl)
            
            annotation1 <- getBM(attributes = c(input$keggChip, 'entrezgene_id'), 
                                 filters = input$keggChip, 
                                 values = top.table1$Probeset.ID, 
                                 mart = ensembl)
            
            colnames(annotation1) <- c("Probeset.ID", "Gene.ID")
            annotation1$Probeset.ID <- as.character(annotation1$Probeset.ID)
            
            top.table <- top.table1 %>%
              inner_join(annotation1, by = c("Probeset.ID" = "Probeset.ID"))
            
            top.table
            
          },
          error=function(cond) {
            return(NULL)
          })
        
        
        
        if (!is.null(top.table)){
          
          #ORA
          if(input$oraPathway == "ora"){
            bgenelist <- top.table$Gene.ID
            
            if(input$raworadjkegg == "raw"){
              siggenelist <- top.table %>%
                filter(P.Value < input$pthresholdkegg) %>%
                filter(abs(logFC) > input$logFCthresholdkegg)
            }
            
            if(input$raworadjkegg == "adj"){
              siggenelist <- top.table %>%
                filter(adj.P.Val < input$pthresholdkegg) %>%
                filter(abs(logFC) > input$logFCthresholdkegg)
            }
            
            siggenelist <- siggenelist$Gene.ID
            
            
            genelist <- list(bgenelist, siggenelist)
          }
          
          
          #GSEA
          if(input$oraPathway == "gsea"){
            if (input$rankingKEGG == "logFC") {
              genelist <- top.table$logFC
            }
            if (input$rankingKEGG == "-log(P-value)") {
              genelist <- -log10(top.table$P.Value)
            }
            
            
            # name the vector
            names(genelist) <- top.table$Gene.ID
            
            # omit any NA values 
            genelist <-na.omit(genelist)
            
            # sort the list in decreasing order (required for clusterProfiler)
            genelist <- sort(genelist, decreasing = TRUE)
          }
        }
      }
      
      
      return(genelist)
    })
    

    
    
      #Results of enrichment analysis
      KEGGresults <- eventReactive(input$kegg.ok, {
        
        KEGGresults <- NULL
        
        if (!is.null((kegggenelist()))){
          if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
            if (species() == "Homo sapiens"){
              
              if(input$oraPathway == "ora"){
                KEGGresults <- enrichKEGG(gene = kegggenelist()[[2]],
                                          universe = kegggenelist()[[1]],
                                          organism = 'hsa',
                                          pvalueCutoff = 1,
                                          qvalueCutoff = 1)
              }
              
              if(input$oraPathway == "gsea"){
                KEGGresults <- gseKEGG(geneList = kegggenelist(),
                                       organism = 'hsa',
                                       pvalueCutoff = 1)
              }
            }
          }
          
          if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
            if (input$keggOrg == "hsapiens_gene_ensembl"){
              
              if(input$oraPathway == "ora"){
                KEGGresults <- enrichKEGG(gene = kegggenelist()[[2]],
                                          universe = kegggenelist()[[1]],
                                          organism = 'hsa',
                                          pvalueCutoff = 1,
                                          qvalueCutoff = 1)
              }
              
              if(input$oraPathway == "gsea"){
                KEGGresults <- gseKEGG(geneList = kegggenelist(),
                                       organism = 'hsa',
                                       pvalueCutoff = 1)
              }
            }
          }

        }

        return(KEGGresults)
        
      })
      
      
      #Error message
      if (is.null(KEGGresults())){
        removeModal()
        output$errorkegg <- renderText({
          "Pathway analysis is unfortunately not possible using the current probeset annotation."
        })
      }
      
      
      if (!is.null(KEGGresults())){
        removeModal()
        
        #Output table
        output$kegg <- renderDataTable({
          
          allRes <- KEGGresults()@result
          
          allRes$ID = paste0(
            '<a ',
            'href=',
            paste(
              "https://www.genome.jp/pathway/",
              allRes$ID,
              sep = ''
            ),
            '>',
            allRes$ID,
            '</a>'
          )
          
          return(allRes)
          
        }, options = list(lengthMenu = c(5, 10, 20, 30), pageLength = 5),
        selection = list(mode = "single", selected = 1), rownames = FALSE, 
        escape = FALSE)
        
        
        #Bar plot
        output$KEGGplot <- renderPlotly({
          req(KEGGresults())
          req(OraGSEA_KEGG())
          
          plotWP <- as.data.frame(KEGGresults())

          if(OraGSEA_KEGG() == "ora"){

            plotWP$GeneRatio <- with(read.table(text = plotWP$GeneRatio, sep = "/"), V1 / V2)
            
            #Arrange based on p-value
            plotWP<- plotWP %>%
              dplyr::arrange(pvalue)
            
            #Get class
            test <- NULL
            for (i in 1:20){
              test[i] <- keggGet(plotWP[i,1])[[1]]$CLASS
            }
            
            plotWP1 <- plotWP[1:20,]
            plotWP1$Class <- test

            
            #Plot
            fig <- ggplot(plotWP1, aes(x = reorder(Description,GeneRatio), y = GeneRatio, 
                                       fill = Class)) + 
              geom_bar(stat="identity") +
              theme_bw(base_size = 14) +
              coord_flip() +
              ylab("Gene ratio") +
              xlab(NULL) +
              theme(legend.position = "none")
            
            return(ggplotly(fig, tooltip = "fill"))
          }
          
          if(OraGSEA_KEGG() == "gsea"){
            #Arrange based on p-value
            plotWP <- plotWP %>%
              dplyr::arrange(pvalue)
            
            #Get class
            test <- NULL
            for (i in 1:20){
              test[i] <- keggGet(plotWP[i,1])[[1]]$CLASS
            }
            
            plotWP1 <- plotWP[1:20,]
            plotWP1$Class <- test
            
            #Plot
            fig <- ggplot(plotWP1, aes(x = reorder(Description,NES), y = NES,
                                       fill = Class)) +
              geom_bar(stat="identity") +
              theme_bw(base_size = 14) +
              coord_flip() +
              ylab("Normalized enrichment score (NES)") +
              xlab(NULL) +
              theme(legend.position = "none")
            
            return(ggplotly(fig, tooltip = "fill"))
          }
        })
        
        #GSEA plot
        output$KEGGplotgsea <- renderPlot({
          req(input$oraPathway)
          if(input$oraPathway == "gsea"){
            gseaplot2(KEGGresults(), geneSetID = input$kegg_rows_selected, 
                      title = KEGGresults()$Description[input$kegg_rows_selected])
          }
        })
      }
  })
  
  #****************************************************************************#
  #   Wikipathways
  #****************************************************************************#
  
  #Select comparison
  output$uiwpcomp <- renderUI({
    pickerInput(
      inputId = "wpcomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  #Select organism
  output$uiwpOrg <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "wpOrg",
        label = "Organism",
        choices = "hsapiens_gene_ensembl",
        options = list(
          style = "btn-danger")
      )
    }
  })
  
  #Select chiptype
  output$uiwpChip <- renderUI({
    if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
      pickerInput(
        inputId = "wpChip",
        label = "Chip type",
        choices = c("affy_hc_g110", "affy_hg_focus", "affy_hg_u133a", "affy_hg_u133a_2", 
                    "affy_hg_u133b", "affy_hg_u133_plus_2", "affy_hg_u95a", "affy_hg_u95av2",
                    "affy_hg_u95b", "affy_hg_u95c", "affy_hg_u95d", "affy_hg_u95d", "affy_hg_u95e",
                    "affy_hta_2_0", "affy_huex_1_0_st_v2", "affy_hugenefl", "affy_hugene_1_0_st_v1",
                    "affy_hugene_2_0_st_v1", "affy_hugene_2_1_st_v1", "affy_primeview", "affy_u133_x3p"),
        selected = get_chiptype(annotation(normData())),
        options = list(
          style = "btn-warning")
      )
    }
  })
  
  output$wp <- renderDataTable(NULL)
  output$WPplot <- renderPlotly(NULL)
  output$errorwp <- renderText(NULL)
  
  
  
  observeEvent(input$wp.ok, {
    
    #Show loading message
    showModal(modalDialog(title = h4("Pathway analysis is being performed....", 
                                     align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", 
                             align = "center"),
                          br()))
    
    wpgenelist <- eventReactive(input$wp.ok,{
      
      genelist <- NULL
      
      #Custom annotation
      if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
        if(CDFtype() == "ENTREZG"){
          top.table <- top.table()[[input$wpcomp]]
          
          
          #ORA
          if(input$oraPathway == "ora"){
            bgenelist <- top.table$Probeset.ID
            
            if(input$raworadjwp == "raw"){
              siggenelist <- top.table %>%
                filter(P.Value < input$pthresholdwp) %>%
                filter(abs(logFC) > input$logFCthresholdwp)
            }
            
            if(input$raworadjwp == "adj"){
              siggenelist <- top.table %>%
                filter(adj.P.Val < input$pthresholdwp) %>%
                filter(abs(logFC) > input$logFCthresholdwp)
            }
            
            siggenelist <- siggenelist$Probeset.ID
            
            
            genelist <- list(bgenelist, siggenelist)
          }

          
          #GSEA
          if(input$oraPathway == "gsea"){
            if (input$rankingWP == "logFC") {
              genelist <- top.table$logFC
            }
            if (input$rankingWP == "-log(P-value)") {
              genelist <- -log10(top.table$P.Value)
            }
            
            # name the vector
            names(genelist) <- top.table$Probeset.ID
            
            # omit any NA values 
            genelist <-na.omit(genelist)
            
            # sort the list in decreasing order (required for clusterProfiler)
            genelist <- sort(genelist, decreasing = TRUE)
          }
          
        }
      }
      
      if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
        
        
          top.table <- tryCatch({
            top.table1 <- top.table()[[input$wpcomp]]
            
            ensembl <- useMart("ensembl")
            ensembl <- useDataset(input$wpOrg,mart=ensembl)
            
            annotation1 <- getBM(attributes = c(input$wpChip, 'entrezgene_id'), 
                                 filters = input$wpChip, 
                                 values = top.table1$Probeset.ID, 
                                 mart = ensembl)
            
            colnames(annotation1) <- c("Probeset.ID", "Gene.ID")
            annotation1$Probeset.ID <- as.character(annotation1$Probeset.ID)
            
            top.table <- top.table1 %>%
              inner_join(annotation1, by = c("Probeset.ID" = "Probeset.ID"))
            
            top.table
            
          },
          error=function(cond) {
            return(NULL)
          })
     
        
        
        if (!is.null(top.table)){
          
          #ORA
          if(input$oraPathway == "ora"){
            bgenelist <- top.table$Gene.ID
            
            if(input$raworadjwp == "raw"){
              siggenelist <- top.table %>%
                filter(P.Value < input$pthresholdwp) %>%
                filter(abs(logFC) > input$logFCthresholdwp)
            }
            
            if(input$raworadjwp == "adj"){
              siggenelist <- top.table %>%
                filter(adj.P.Val < input$pthresholdwp) %>%
                filter(abs(logFC) > input$logFCthresholdwp)
            }
            
            siggenelist <- siggenelist$Gene.ID
            
            
            genelist <- list(bgenelist, siggenelist)
          }
          
          
          #GSEA
          if(input$oraPathway == "gsea"){
            if (input$rankingWP == "logFC") {
              genelist <- top.table$logFC
            }
            if (input$rankingWP == "-log(P-value)") {
              genelist <- -log10(top.table$P.Value)
            }
            
            # name the vector
            names(genelist) <- top.table$Gene.ID
            
            # omit any NA values 
            genelist <-na.omit(genelist)
            
            # sort the list in decreasing order (required for clusterProfiler)
            genelist <- sort(genelist, decreasing = TRUE)
          }
        }
      }
      
      return(genelist)
    })
    
    
    
    
    #Results of enrichment analysis
    WPresults <-eventReactive(input$wp.ok, {
      
      WPresults <- NULL
      
      if (!is.null((wpgenelist()))){
        if ((PSannotation() == "Custom annotations") & ((rawornormalized() == "Raw") | (!is.null(normData.saved())))){
          if (species() == "Homo sapiens"){
            if(input$oraPathway == "ora"){
              WPresults <- enrichWP(gene = wpgenelist()[[2]],
                                    universe = wpgenelist()[[1]],
                                    organism = "Homo sapiens",
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1)
            }
            
            if(input$oraPathway == "gsea"){
              WPresults <- gseWP(geneList = wpgenelist(),
                                 organism = "Homo sapiens",
                                 pvalueCutoff = 1)
            }
          }
        }
        if ((PSannotation() == "No annotations") | ((rawornormalized() == "Normalized") & (is.null(normData.saved())))){
          if (input$wpOrg == "hsapiens_gene_ensembl"){
            if(input$oraPathway == "ora"){
              WPresults <- enrichWP(gene = wpgenelist()[[2]],
                                    universe = wpgenelist()[[1]],
                                    organism = "Homo sapiens",
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1)
            }
            
            if(input$oraPathway == "gsea"){
              WPresults <- gseWP(geneList = wpgenelist(),
                                 organism = "Homo sapiens",
                                 pvalueCutoff = 1)
            }
          }
        }
        

      }
      
      return(WPresults)
      
    })
    
    
    #Error message
    if (is.null(WPresults())){
      removeModal()
      output$errorwp <- renderText({
        "Pathway analysis is unfortunately not possible using the current probeset annotation."
      })
    }
    
    
    if (!is.null(WPresults())){
      removeModal()
      
      #Output table
      output$wp <- DT::renderDataTable({
        
        allRes <- WPresults()@result
        allRes$ID = paste0(
          '<a ',
          'href=',
          paste(
            "https://www.wikipathways.org/index.php/Pathway:",
            allRes$ID,
            sep = ''
          ),
          '>',
          allRes$ID,
          '</a>'
        )
        return(allRes)
        
      }, options = list(lengthMenu = c(5, 10, 20, 30), pageLength = 5),  
      selection = list(mode = "single", selected = 1), rownames = FALSE, 
      escape = FALSE)
      
      #Bar plot
      output$WPplot <- renderPlotly({
        
        plotWP <- as.data.frame(WPresults())
        
        if(OraGSEA_WP() == "ora"){
          plotWP$GeneRatio <- with(read.table(text = plotWP$GeneRatio, sep = "/"), V1 / V2)
          
          #Arrange based on p-value
          plotWP<- plotWP %>%
            dplyr::arrange(pvalue)
          
          #Get the first 20 elements
          plotWP1 <- plotWP[1:20,]
          
          #Plot
          fig <- ggplot(plotWP1, aes(x = reorder(Description,GeneRatio), y = GeneRatio, 
                                     fill = p.adjust)) + 
            geom_bar(stat="identity") +
            theme_bw(base_size = 14) +
            coord_flip() +
            scale_fill_continuous(low="blue", high="red") +
            ylab("Gene ratio") +
            xlab(NULL)
          
          return(ggplotly(fig, tooltip = "fill"))
        }
        
        if(OraGSEA_WP() == "gsea"){
          #Arrange based on p-value
          plotWP <- plotWP %>%
            dplyr::arrange(pvalue)
          
          #Get the first 20 elements
          plotWP1 <- plotWP[1:20,]
          
          #Plot
          fig <- ggplot(plotWP1, aes(x = reorder(Description,NES), y = NES,
                                     fill = p.adjust)) +
            geom_bar(stat="identity") +
            theme_bw(base_size = 14) +
            coord_flip() +
            scale_fill_continuous(low="blue", high="red") +
            ylab("Normalized enrichment score (NES)") +
            xlab(NULL)
          
          return(ggplotly(fig, tooltip = "fill"))
        }
      })
      
      #GSEA plot
      output$WPplotgsea <- renderPlot({
        req(input$oraPathway)
        if(input$oraPathway == "gsea"){
          gseaplot2(WPresults(), geneSetID = input$wp_rows_selected, 
                    title = WPresults()$Description[input$wp_rows_selected])
        }
      })
    }
  })
  
  
  ##############################################################################
  #Parallel coordinates plot
  ##############################################################################
  
  #Select comparison
  output$uiparacomp <- renderUI({
    pickerInput(
      inputId = "paracomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  

  
  output$allpara <- renderPlotly({
    
    req(input$paracomp)
    
    top.table1 <- top.table()[[input$paracomp]]
    
    fig1 <- top.table1 %>% plot_ly(type = 'parcoords', 
                                   line = list(color = ~t, colorscale = "Jet"),
                                   dimensions = list(
                                     list(label = "-log(P-value)", 
                                          values = -log10(top.table1$P.Value),
                                          constraintrange = c(1.3,max(-log10(top.table1$P.Value)))),
                                     list(label = "logFC", 
                                          values = ~logFC),
                                     list(label = "AveExpr", 
                                          values = ~AveExpr)
                                     
                                   )
    )
    
    return(fig1)
  })
  
  
  
  #Parrallel coordinates plot
  output$topfeaturepara <- renderPlotly({
    
    req(input$paracomp)
    
    comparison <- input$paracomp
    comparison <- strsplit(comparison, " - ")[[1]]
    
    meta2 <- meta1()
    meta2$Grouping <- make.names(meta2$Grouping)
    
    samples <- meta2 %>%
      filter(Grouping == comparison[1] | Grouping == comparison[2])
    
    
    topfeatures <- head(top.table()[[input$paracomp]], 10)
    data.expr1 <- data.expr()[topfeatures$Probeset.ID,]
    data.expr1 <- data.expr1[,samples$names]
    
    experimentFactor <- factor(samples$Grouping)
    sidecolors <- data.frame(experimentFactor)
    colnames(sidecolors) <- "Grouping"
    
    
    plot_ly(type = 'parcoords', 
            line = list(color = as.numeric(sidecolors$Grouping), colorscale = "Jet"),
            dimensions = list(
              list(range = c(min(data.expr1[1,]),max(data.expr1[1,])),
                   label = rownames(data.expr1)[1], 
                   values = as.vector(data.expr1[1,])),
              list(range = c(min(data.expr1[2,]),max(data.expr1[2,])),
                   label = rownames(data.expr1)[2], 
                   values = as.vector(data.expr1[2,])),
              list(range = c(min(data.expr1[3,]),max(data.expr1[3,])),
                   label = rownames(data.expr1)[3], 
                   values = as.vector(data.expr1[3,])),
              list(range = c(min(data.expr1[4,]),max(data.expr1[4,])),
                   label = rownames(data.expr1)[4], 
                   values = as.vector(data.expr1[4,])),
              list(range = c(min(data.expr1[5,]),max(data.expr1[5,])),
                   label = rownames(data.expr1)[5], 
                   values = as.vector(data.expr1[5,])),
              list(range = c(min(data.expr1[6,]),max(data.expr1[6,])),
                   label = rownames(data.expr1)[6], 
                   values = as.vector(data.expr1[6,])),
              list(range = c(min(data.expr1[7,]),max(data.expr1[7,])),
                   label = rownames(data.expr1)[7], 
                   values = as.vector(data.expr1[7,])),
              list(range = c(min(data.expr1[8,]),max(data.expr1[8,])),
                   label = rownames(data.expr1)[8], 
                   values = as.vector(data.expr1[8,])),
              list(range = c(min(data.expr1[9,]),max(data.expr1[9,])),
                   label = rownames(data.expr1)[9], 
                   values = as.vector(data.expr1[9,])),
              list(range = c(min(data.expr1[10,]),max(data.expr1[10,])),
                   label = rownames(data.expr1)[10], 
                   values = as.vector(data.expr1[10,]))
              
            )
    )
    
  })
  
  ##############################################################################
  #Pie chart
  ##############################################################################
  
  #Select comparison
  output$uipiecomp <- renderUI({
    pickerInput(
      inputId = "comp_pie",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  

  output$piechart <- renderPlotly({
    
    #Requirements
    req(input$comp_pie)
    req(input$pchoice_pie)
    req(input$pthreshold_pie)
    req(input$logFCthreshold_pie)
    
    #Settings
    comp <- input$comp_pie
    p.choice <- input$pchoice_pie
    p.threshold <- input$pthreshold_pie
    logFC.threshold <- input$logFCthreshold_pie
    
    #Top.table
    top.table <- top.table()[[comp]]
    
    #Filtered table
    if (p.choice == "raw"){
      top.table$Class[(top.table$logFC < logFC.threshold & top.table$logFC > (-1 * logFC.threshold)) | top.table$P.Value > p.threshold] = "Insignificant"
      top.table$Class[top.table$logFC < (-1 * logFC.threshold) & top.table$P.Value <= p.threshold] = "Downregulated"
      top.table$Class[top.table$logFC >= logFC.threshold & top.table$P.Value <= p.threshold] = "Upregulated"
    }
    
    if (p.choice == "adj"){
      top.table$Class[(top.table$logFC < logFC.threshold & top.table$logFC > (-1 * logFC.threshold)) | top.table$adj.P.Val > p.threshold] = "Insignificant"
      top.table$Class[top.table$logFC < (-1 * logFC.threshold) & top.table$adj.P.Val <= p.threshold] = "Downregulated"
      top.table$Class[top.table$logFC >= logFC.threshold & top.table$adj.P.Val <= p.threshold] = "Upregulated"
      
    }
    
    
    df <- data.frame("Category" = names(table(top.table$Class)), table(top.table$Class))
    
    fig <- plot_ly(df, labels = ~Category, 
                   values = ~Freq, type = 'pie', textposition = 'inside',
                   insidetextfont = list(color = '#FFFFFF'),
                   marker = list(colors = c( 'blue', 'darkgrey', 'red'),
                                 line = list(color = '#FFFFFF', width = 1)),
                   showlegend = TRUE)
    fig <- fig %>% layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
    fig
  })
 
  
  
  ##############################################################################
  #Remainder
  ##############################################################################

  
  
  
  output$topfeatureradar <- renderPlotly({
    
    req(input$paracomp)
    
    comparison <- input$paracomp
    comparison <- strsplit(comparison, " - ")[[1]]
    
    meta2 <- meta1()
    meta2$Grouping <- make.names(meta2$Grouping)
    
    samples <- meta2 %>%
      filter(Grouping == comparison[1] | Grouping == comparison[2])
    
    
    topfeatures <- head(top.table()[[input$paracomp]], 5L)
    data.expr1 <- data.expr()[topfeatures$Probeset.ID,]
    data.expr1 <- data.expr1[,samples$names]
    
    experimentFactor <- factor(samples$Grouping)
    sidecolors <- data.frame(experimentFactor)
    colnames(sidecolors) <- "Grouping"
    
    
    Group1 <- samples %>%
      filter(Grouping == comparison[1])
    
    Group1 <- Group1$names
    
    Group2 <- samples %>%
      filter(Grouping == comparison[2])
    
    Group2 <- Group2$names
    
    fig <- plot_ly(
      type = 'scatterpolar',
      fill = 'toself'
    ) 
    fig <- fig %>%
      add_trace(
        r = c(mean(data.expr1[1,Group1]), 
              mean(data.expr1[2,Group1]), 
              mean(data.expr1[3,Group1]), 
              mean(data.expr1[4,Group1]), 
              mean(data.expr1[5,Group1]),
              mean(data.expr1[1,Group1])),
        
        theta = paste0("Probeset: ", c(rownames(data.expr1), rownames(data.expr1)[1])),
        name = comparison[1],
        marker = list(color = "blue"),
        line = list(color = "blue"),
        fillcolor = "blue",
        opacity = 0.5
      ) 
    fig <- fig %>%
      add_trace(
        r = c(mean(data.expr1[1,Group2]), 
              mean(data.expr1[2,Group2]), 
              mean(data.expr1[3,Group2]), 
              mean(data.expr1[4,Group2]), 
              mean(data.expr1[5,Group2]),
              mean(data.expr1[1,Group2])),
        theta = paste0("Probeset: ", c(rownames(data.expr1), rownames(data.expr1)[1])),
        name = comparison[2],
        marker = list(color = "red"),
        line = list(color = "red"),
        fillcolor = "red",
        opacity = 0.5
      ) 
    fig <- fig %>%
      layout(
        polar = list(
          radialaxis = list(
            visible = T,
            range = c(min(data.expr1),max(data.expr1))
          )
        )
      )
    
    return(fig)
    
  })
  

  
}




