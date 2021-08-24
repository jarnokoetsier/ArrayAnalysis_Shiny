

###############################################################################################################################

#SERVER

###############################################################################################################################
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

server <- function(input, output, session){
  options(shiny.maxRequestSize=50*1024^2)
  
  ################################################################################################################################
  #data accession
  ################################################################################################################################
  
  #Hide tabs
  hideTab("navbar", target = "panel2")
  hideTab("navbar", target = "panel3")
  hideTab("navbar", target = "panel4")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel5")
  hideTab("navbar", target = "panel6")
  
  
  #Welcome message
  sendSweetAlert(
    session = session,
    title = "ArrayAnalysis",
    text = "Welcome to the ArrayAnalysis app!",
    btn_labels = "Click here to start"
  )
  
  
  #Info panel
  observeEvent(input$infopanel1, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "ArrayAnalysis has been developed for the user-friendly analysis of microarray data. 
      Click on the documentation tab for more information.",
      type = "info"
    )
  })
  
  
  #Select data access (i.e., GEO/ArrayExpress database or upload own CELs)
  output$getdatabaseout <- renderUI({
    req(input$database)
    if (input$database != "Upload CELs"){
      textInput(inputId = "getGEO", 
                label = NULL, 
                value = "",
                placeholder = "Accession number")
      
    }
  })
  
  output$uploadcelsout <- renderUI({
    req(input$database)
    if (input$database == "Upload CELs"){
      fileInput(inputId = "uploadcelsin",
                label = NULL,
                accept = ".zip",
                placeholder = "Select .zip data file")
    }
  })
  
  
  #Get accession and database
  database <- eventReactive(input$downloaddata, {
    input$database
   
  })
  
  
  accession <- eventReactive(input$downloaddata, {
    input$getGEO
  })
  

  #Example
  observeEvent(input$example, {
    
    if (input$database == "GEO"){
      updateTextInput(session, "getGEO",
                      label = NULL,
                      value = "GSE36980")
      
    }
    
    if (input$database == "ArrayExpress"){
      updateTextInput(session, "getGEO",
                      label = NULL,
                      value = "E-MTAB-9988")
    }
    
    
  })
  
  
  #Download data for retrieval of meta data
  gset <- eventReactive(input$downloaddata, {
    if (database() == "GEO"){
      gset <- getGEO(accession(), GSEMatrix =TRUE, getGPL = FALSE)
      
    }
    if (database() == "ArrayExpress"){
      gset <- ArrayExpress(accession())
      
    }
    
    if (database() == "Upload CELs"){
      gset <- NULL
      
    }
    
    return(gset)
  })
  
  
  #Loading message
  observeEvent(input$downloaddata, {
    
    if (database() != "Upload CELs"){
      showModal(modalDialog(
        title = h4(paste0(accession(), " dataset from the ", database(), " database is being loaded...."), align = "center"), 
        footer = NULL,
        h5("Please be patient. This might take a while.", align = "center")
      ))
        
                           
      if(!is.null(gset())) {
        get_grouping(gset(), database())
        
      }
      removeModal()
      
    }
    
  })
  
  
  #Get data path
  zzip <- eventReactive(input$downloaddata, {
    input$uploadcelsin$datapath
  })
  
  
  #Check cel files
  celfiles <- eventReactive(input$downloaddata, {
    
    database = database()
    
    if (database == "GEO"){
      
      accession = accession()
      exdir = paste0("data_", accession)
      
      if (!file.exists(exdir)){
        
        getGEOSuppFiles(accession)
        tarfile = paste0(accession, "/", accession, "_RAW.tar")
        untar(tarfile, exdir = exdir)
        
      }
      
      
      celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", full.names = TRUE)
      
    }
    
    if (database == "ArrayExpress"){
      
      accession = accession()
      exdir = paste0("data_", accession)
      
      if (!file.exists(exdir)) {
        getAE(accession, type = "raw", extract = FALSE)
        unzip(paste0(accession, ".raw.1.zip"), exdir = exdir)
      }
      
      
      celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", full.names = TRUE)
      
    }
    
    if (database == "Upload CELs") {
      zzip <- zzip()
      celfiles <- unzip(zzip)
     
    }
    
    
    return(celfiles)
    
  })
  


  #Data successfully downloaded
  observeEvent(if (length(celfiles()) > 0) TRUE, {
    sendSweetAlert(
      session = session,
      title = "Success!!",
      text = "Dataset successfully selected! Now you can group the samples.",
      type = "success")
    
  })
  
  
  #Incorrect dataset selected
  observeEvent(if (length(celfiles()) <= 0) TRUE, {
    sendSweetAlert(
      session = session,
      title = "Error!",
      text = "No valid dataset selected. Please try again.",
      type = "error")
    
  })
  
  
  #go to next tab
  observeEvent(if (length(celfiles()) > 0) TRUE, {
    showTab("navbar", target = "panel2")
  })
  
  
  observeEvent(if (length(celfiles()) > 0) TRUE,  {
    updateNavbarPage(session, "navbar",
                     selected = "panel2")
    
  })
  

  
  #Enter ID to retrieve saved data
  observeEvent(input$continue, {
    inputSweetAlert(
      session = session,
      inputId = "idcontinue",
      title = "Enter ID",
      text = "After you saved the data, you received an ID. This ID is required to retrieve the saved data.",
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
  

  
  #get saved data
  meta1.saved <- reactive({
    meta1 <- NULL
    
    if(length(input$idcontinue) > 0){
      if (file.exists(paste0("meta1", input$idcontinue))){
        meta1 <- readRDS(paste0("meta1", input$idcontinue))
      }
    }
    return(meta1)
  })
  
  normData.saved <- reactive({
    normData <- NULL
    
    if(length(input$idcontinue) > 0){
      if (file.exists(paste0("normData", input$idcontinue))){
        normData <- readRDS(paste0("normData", input$idcontinue))
      }
      
    }
    return(normData)
  })
  
  
  
  ################################################################################################################################
  #get meta
  ################################################################################################################################
  

  
  output$makegroupsui <- renderUI({
    
    if (!is.null(gset())){
      choices = c("Dataset", "Description file", "Manual grouping")
      selected = "Dataset"
    }
     
    if (is.null(gset())){
      choices = c("Description file", "Manual grouping")
      selected = "Description file"
    }   
    
    radioGroupButtons(
      inputId = "makegroups",
      label = "Make sample groups from:", 
      choices = choices,
      selected = selected,
      status = "danger")
  
  })
 
  
  #select Grouping variable(s)
  
  output$groups <- renderUI({
    
    if (length(input$makegroups) > 0){
      if (input$makegroups == "Dataset"){
        checkboxGroupInput(inputId = "groupselect", 
                           label = "Select grouping variable(s):", 
                           choices = colnames(get_grouping(gset(), database = database())),
                           selected = auto_group(get_grouping(gset(), database = database())))
      }
    }
   
  })
  
  
  #Make own sample groups
  
  samplelist <- reactive({
    if(database() == "GEO"){
      
      samplelist <- gset()[[paste0(accession(),"_series_matrix.txt.gz")]]@phenoData@data[["geo_accession"]]
      
    }
    
    if(database() == "ArrayExpress"){
      
      samplelist <- rownames(gset@phenoData@data)
      
    }
    
    if (database() == "Upload CELs"){ 
      
      samplelist <- gsub(".*/","",celfiles())
      
    }
    
    return(samplelist)
  })
  
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
            selected_header = "Group 2"
            
          )
        )
      }
    }
  })
  
  
  output$uidescriptionfile <- renderUI({
    
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Description file") {
        
        fileInput(inputId = "descriptionfile",
                  label = NULL,
                  accept = ".xlsx",
                  placeholder = "Select excel description file")
      
      }
    }
  })
  
  output$textdescription <- renderText({
    "<p> The description file needs to be an excel file containing two columns:<br /> 
    <p> 1) The first column should contain the <b> sample IDs</b>. Note that these IDs should match to the CEL file names.<br />
    <p> 2) The second column should include the <b> grouping </b> of the samples. <br />
    <p> <b> IMPORTANT: </b> Both columns should have an header. <br />
    <br />"
  })
  
  
  
  output$uitextdescription <- renderUI({
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Description file") {
        
        htmlOutput("textdescription")
        
      }
    }

  })
  
  output$downloadtemplate <- downloadHandler(
    filename = "template.xlsx",
    content = function(file) {
      file.copy("template.xlsx", file)
    }
  )
  
  output$uidownloadtemplate <- renderUI({
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Description file") {
        
        downloadButton("downloadtemplate", "Download template")
        
      }
    }
    
  })
  
  
  #Make "meta" dataframe
  meta <- reactive({
    
    meta <- NULL
    
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Dataset"){
        
        if (length(input$groupselect) > 0){
          groups <- get_grouping(gset(), database = database())
          meta <- get_meta(gset = gset(), grouping_column = groups[,input$groupselect], database = database())
        }
      }
      
      if (input$makegroups == "Manual grouping"){
        
        group2 <- cbind(input$owngroupsin, rep("group 2", length(input$owngroupsin)))
        
        group1 <- setdiff(samplelist(), input$owngroupsin)
        group1 <- cbind(group1, rep("group 1", length(group1)))
        
        meta <- as.data.frame(rbind(group1, group2))
        rownames(meta) <- NULL
        colnames(meta) <- c("Sample.ID", "Grouping")
      }
      
      if (input$makegroups == "Description file"){
        if (length(input$descriptionfile$datapath) > 0){
          meta <- read_excel(input$descriptionfile$datapath)
          rownames(meta) <- NULL
          colnames(meta) <- c("Sample.ID", "Grouping")
        }
      }
    }
    
    return(meta)
  })
  
  
  
  #make table of meta data
  output$grouping <- renderDataTable({
    
    meta()
    
  })
  
  
  #Error and success message
  observeEvent(input$meta.ok, {
    samplelist <- gsub(".*/","",celfiles())
    samplelist <- as.data.frame(samplelist)
    colnames(samplelist) <- "names"
    
    test <- samplelist %>% 
      fuzzyjoin::fuzzy_inner_join(meta(), by = c("names" = "Sample.ID"), match_fun = str_detect)
    
    
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
        text = "Samples are successfully grouped! Now you can start pre-processing the data.",
        type = "success"
      )
      
      updateNavbarPage(session, "navbar",
                       selected = "panel3")
      
      showTab("navbar", target = "panel3")
    }

  })
  
  #download table of meta data
  output$downloadmeta <- downloadHandler(
    filename = "meta data",
    content = function(file){
      write.table(meta(), file)
    }
  )
  
  #Go to previous tab
  observeEvent(input$meta.back, {
    updateNavbarPage(session, "navbar",
                     selected = "panel1")
    
    
  })
  
  
  ################################################################################################################################
  #Pre-processing
  ################################################################################################################################
  

  
  
  #Remove outliers
  
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

  
  outliers <- eventReactive(input$preprocessing, {
    if (input$outlier == TRUE){
      outliers = NULL
    }
    
    if (input$outlier == FALSE){
      if (length(input$outliersin) < 1){
        outliers = NULL
      }
    }
    
    if (input$outlier == FALSE){
      if (length(input$outliersin) > 0){
        outliers = input$outliersin
      }
    }
    
    return(outliers)
  })
  
  
  #rawData
  rawData <- eventReactive(input$preprocessing, {
    
    showModal(modalDialog(title = h4("Dataset is being pre-processed....", align = "center"), 
                          footer = NULL,
                          h5("Please be patient. This might take a while.", align = "center"),
                          br()))
   
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
    
    rawData = ReadAffy(filenames=celfiles)
    
    return(rawData)
    
  })
  
  
  #meta1
  meta1 <- reactive({
    
    if (is.null(meta1.saved())){
      meta <- meta()
      rawData <- rawData()
      
      
      #add CEL names to meta
      names <- as.data.frame(sampleNames(rawData))
      colnames(names) <- "names"
      
      meta1 <- names %>% 
        fuzzyjoin::fuzzy_inner_join(meta, by = c("names" = "Sample.ID"), match_fun = str_detect)
    }
    
    if (!is.null(meta1.saved())){
      meta1 <- meta1.saved()
    }

    return(meta1)
    
  })
  
  experimentFactor <- reactive({
    meta1 <- meta1()
    experimentFactor <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
    return(experimentFactor)
  })
  
  
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
  
  
  # Normalize data
  normData <- eventReactive((if(!is.null(normData.saved())) TRUE) | input$preprocessing, {
    
    if (is.null(normData.saved())){
      normMeth <- input$normMeth
      CDFtype <- input$CDFtype
      species <- input$species
      
      experimentFactor <- experimentFactor()
      
      ifelse(input$annotations=="Custom annotations",customCDF<-TRUE,customCDF<-FALSE)
      ifelse(input$annotations=="Upload annotation file",uploadCDF<-TRUE,uploadCDF<-FALSE)
      ifelse(input$perGroup == "Use all arrays",perGroup<-FALSE,perGroup<-TRUE)
      
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
        if(nGroups==1) warning("normalization per group requested, but no groups indicated in data set")
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
                   probeLibrary <- tolower(paste(Data@annotation,species,CDFtype,"probe",sep=""))
                   loaded <- suppressWarnings(try(eval(parse("",-1,paste("library(",probeLibrary,")", sep=""))),TRUE))
                   if(class(loaded)=="try-error") {
                     install.packages(probeLibrary, repos="http://brainarray.mbni.med.umich.edu/bioc")
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
                 normData.tmp <- justPlier(Data.tmp, normalize=TRUE, norm.type = ntype)
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
    }
    
    
    if (!is.null(normData.saved())){
      normData <- normData.saved()
    }

    return(normData)
  })
  
  
  data.expr <- eventReactive(if (length(normData()) > 0) TRUE, {
    
    data.expr <- exprs(normData())
    removeModal()
    
    return(data.expr)
    
  })
  
  
  
  
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
  
  
  
  output$correlNorm <- renderImage({
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
  
  
  output$normhist <- renderPlot(NULL)
  
  
  
  
  observeEvent(if (length(data.expr()) > 0) TRUE, {
    
    #########################
    ###Norm boxplot###
    #########################
    
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
      suppressWarnings(boxplot(normData(), col=plotColors ,main=tmain, axes=FALSE, pch = 20, cex=0.7))
      if(length(levels(experimentFactor()))>1){ 
        legend("topright", levels(experimentFactor()),
               col=legendColors,fill=legendColors, cex = 0.7, bg = "white", bty = "o")
      }
      if(length(sampleNames(normData()))<MAXARRAY){
        cexval <- 0.65
      }else{
        cexval <- 0.45
      }  
      axis(1,at=1:length(sampleNames(normData())),las=2,labels=sampleNames(normData()), cex.axis=cexval)        
      axis(2, cex.axis=0.7)  
      mtext(tmtext2, side=2, cex=0.8)  	   
      mtext("Distributions should be comparable between arrays\n", side=3, font=1, 
            cex=0.7)
      dev.off()
      
      list(src = DataBoxplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    },deleteFile = TRUE)
    
    
    #########################
    ###Raw boxplot###
    #########################
    
    output$boxplotRaw<-renderImage({
      
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
      suppressWarnings(boxplot(rawData, col=plotColors ,main=tmain, axes=FALSE, pch = 20, cex=0.7))
      if(length(levels(experimentFactor()))>1){ 
        legend("topright", levels(experimentFactor()),
               col=legendColors,fill=legendColors, cex = 0.7, bg = "white", bty = "o")
      }
      if(length(sampleNames(rawData))<MAXARRAY){
        cexval <- 0.65
      }else{
        cexval <- 0.45
      }  
      axis(1,at=1:length(sampleNames(rawData)),las=2,labels=sampleNames(rawData), cex.axis=cexval)        
      axis(2, cex.axis=0.7)  
      mtext(tmtext2, side=2, cex=0.8)  	   
      mtext("Distributions should be comparable between arrays\n", side=3, font=1, 
            cex=0.7)
      dev.off()
      
      list(src = DataBoxplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    },deleteFile = TRUE)
    
    
    #########################
    ###Heatmap
    #########################
    
    output$correlNorm <- renderImage({
      
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
      clusterOption1 <- input$clusteroption1
      clusterOption2 <- input$clusteroption2 
      normMeth <- input$normMeth
      
      
      Type <- "Norm"
      text1 <- paste("Array correlation plot\nafter",normMeth,"normalization")
      
      if(length(sampleNames(normData()))<2) {
        warning("Only one array in dataset, no correlation plot made")
      } else {
        normdataCorrelation <- tempfile(fileext = ".png")
        png(file = normdataCorrelation,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
        if(length(sampleNames(normData()))<MAXARRAY) {
          par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
          #subval <- 10
        } else {
          par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
          #subval <- 16
        }        
        
        #note: for computing array correlation, euclidean would not make sense
        #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
        COpt1 <- "pearson"
        if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
        crp <- cor(data.expr(), use="complete.obs", method=COpt1)
        
        text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)
        
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
        
        sideColors <- legendColors[as.numeric(experimentFactor())]
        
        heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
                  main=text1, dendrogram="row", ColSideColors=sideColors)
        
        dev.off()
        list(src = normdataCorrelation,
             width = WIDTH,
             height = HEIGHT,
             alt = "This is alternate text")
      }
    }, deleteFile = TRUE)
    
    
    #########################
    ###Normalized boxplot### GGPLOT
    #########################
    output$normboxplot <- renderPlot({
      
      
      data.expr <- data.expr()
      meta1 <- meta1()
      
      #make dataframe with samples + intensities
      logData <- cbind.data.frame(as.vector(data.expr), rep(colnames(data.expr), each = nrow(data.expr)))
      colnames(logData) <- c("logInt", "sampleName")
      
      logData <- logData %>% inner_join(meta1, by = c("sampleName" = "names"))
      
      logData$Grouping <- factor(logData$Grouping, levels = unique(meta1$Grouping))
      logData$Sample.ID <- factor(logData$Sample.ID, levels = unique(meta1$Sample.ID))
      
      #Make colours
      myPalette <- colorsByFactor(experimentFactor())

      #Make Boxplot
      
      ggplot(logData, aes(x = Sample.ID, y = logInt, fill = Sample.ID, shape = Grouping)) +
        geom_boxplot() +
        scale_fill_manual(values = myPalette$plotColors) +
        labs(title = "Boxplot of normalized intensities",
             subtitle = "Distributions should be comparable between arrays") +
        xlab(NULL) +
        ylab("Normalized log intensity") +
        theme_classic() +
        theme(axis.text.x=element_text(angle=90),
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.background = element_rect(fill="#d3d3d3", colour = "black"),
              legend.title = element_text(face = "bold", hjust = 0.5)) +
        guides(shape=guide_legend(title="Grouping", override.aes = list(fill = myPalette$legendColors)), fill = "none")
    })
      

    
    #Normalized histogram
    output$normhist <- renderPlotly({
      
      data.expr <- data.expr()
      meta1 <- meta1()
      
      #make dataframe with samples + intensities
      logData <- cbind.data.frame(as.vector(data.expr), rep(colnames(data.expr), each = nrow(data.expr)))
      colnames(logData) <- c("logInt", "sampleName")
      
      logData <- logData %>% inner_join(meta1, by = c("sampleName" = "names"))
      
      logData$Grouping <- factor(logData$Grouping, levels = unique(meta1$Grouping))
      logData$Sample.ID <- factor(logData$Sample.ID, levels = unique(meta1$Sample.ID))

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
        guides(shape=guide_legend(title="Grouping", override.aes = list(colour = myPalette$legendColors)), colour = "none")
      
      
      return(ggplotly(plot))
      
    })
  })
  
  
  
  
  #After normalization....
  
  observeEvent(if (length(data.expr()) > 0) TRUE, {
    
    
    #Download data expr
    output$downloadexpr <- downloadHandler(
      filename = "data expr",
      content = function(file){
        write.table(data.expr(), file, sep = "\t")
      })
    
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "Success!!",
      text = "Data successfully pre-processed! Please wait for the QC plots to be rendered.",
      type = "success"
    )
  })
  
  observeEvent(input$ann.back, {
    if (length(meta()) > 0) {
      updateNavbarPage(session, "navbar",
                       selected = "panel2")
    }

  })
  
  
  observeEvent((input$ann.proceed), {
    if (length(data.expr()) > 0){
      updateNavbarPage(session, "navbar",
                       selected = "panel4")
    }
    
  })
  
  observeEvent(input$ann.proceed, {
    if (length(data.expr()) > 0){
      showTab("navbar", target = "panel4")
    }
  })
  
  
  observeEvent(input$save, {
    
    if (length(data.expr()) > 0){
      id <- str_remove_all(Sys.time(), "-")
      id <- str_remove_all(id, " ")
      id <- str_remove_all(id, ":")
      id <- paste0("aa", id)
      
      id.normData <- paste0("normData", id)
      id.meta1 <- paste0("meta1", id)
      
      saveRDS(normData(), id.normData)
      saveRDS(meta1(), id.meta1)
      
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

  
  
  ################################################################################################################################
  #PCA
  ################################################################################################################################
  
  data.PC <- eventReactive(input$ann.proceed, {
    prcomp(t(data.expr()),scale.=TRUE)
  })
  
  output$uizpca <- renderUI({
    if (input$plot3d == TRUE) {
      selectInput(inputId = "zpca", 
                  label = "z-axis",
                  choices = c("PC1","PC2","PC3", "PC4", "PC5", "PC6", "PC7", "PC8"),
                  selected = "PC3")
    }
  })
  
  
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
  
  
  #make reactive for input
  
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
  
  observeEvent(input$pca.back, {
    updateNavbarPage(session, "navbar",
                     selected = "panel3")
    
    
  })
  
  observeEvent(input$pca.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel5")
    
  })
  
  observeEvent(input$pca.ok, {
    showTab("navbar", target = "panel5")
  })
  
  ################################################################################################################################
  #Statistical analysis
  ################################################################################################################################
  
  
  #get statistics
  
  top.table <- reactive({
    
    top.table <- diff_expr(data.expr(), meta1(), comparisons = get_contrasts(meta1()))
    
    return(top.table)
  })
  
  
  ################################################################################################################################
  #Top table
  ################################################################################################################################
  
  
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
  
  bothdirectionstop <- reactive({
    req(input$toptablecomp)
    
    contrast <- input$toptablecomp
    
    switch <- strsplit(contrast, split= " - ")
    
    contrast2 <- paste0(switch[[1]][2], " - ", switch[[1]][1])
    
    bothdirections <- c(contrast, contrast2)
    
    return(bothdirections)
  })
  
  
  
  output$uitoptabledir <- renderUI({
    pickerInput(
      inputId = "toptabledir",
      label = "Direction of comparison",
      choices = bothdirectionstop(),
      options = list(
        style = "btn-info")
    )
  })

  #comparison of interest
  toptablecomp <- reactive(
    input$toptablecomp
  )
  
  toptabledir <- reactive(
    input$toptabledir
  )
  
  
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
  output$finaltable <- renderDataTable({
   toptable()
  })
  
 
  #Download table
  output$downloadfinal <- downloadHandler(
    filename = "top table",
    content = function(file){
      write.table(toptable(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  
 
  ################################################################################################################################
  #P-value analysis
  ################################################################################################################################
  
  
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
 
  ################################################################################################################################
  #FC analysis
  ################################################################################################################################
  
  
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
  
  
  ################################################################################################################################
  #volcano plot
  ################################################################################################################################
  
  
  output$uivolcanocomp <- renderUI({
    pickerInput(
      inputId = "volcanocomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  bothdirectionsvol <- reactive({
    req(input$volcanocomp)
    
    contrast <- input$volcanocomp
    
    switch <- strsplit(contrast, split= " - ")
    
    contrast2 <- paste0(switch[[1]][2], " - ", switch[[1]][1])
    
    bothdirections <- c(contrast, contrast2)
    
    return(bothdirections)
  })
  

  
  output$uivolcanodir <- renderUI({
    pickerInput(
      inputId = "volcanodir",
      label = "Direction of comparison",
      choices = bothdirectionsvol(),
      options = list(
        style = "btn-info")
    )
  })
  
  
  
  volcanocomp <- reactive({
    input$volcanocomp
  })
  
  
  p.choice <- reactive({
    input$raworadj
  })
  
  
  p.threshold <- reactive({
    input$pthreshold
  })
  
  logFC.threshold <- reactive({
    input$logfcthreshold
  })
  
  volcanodir <- reactive({
    input$volcanodir
  })
  
  
  

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
  
    
  output$voltable <- renderDataTable({
    voltable()
  })
  

  output$downloadvol <- downloadHandler(
    filename = "volcano plot table",
    content = function(file){
      write.table(voltable(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  ################################################################################################################################
  #venn diagram
  ################################################################################################################################
  
  
  
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
  
  
  output$legendvenn <- renderText({
    if (length(input$contrastselection) == 2) {
      text <- paste0("<p> <b> A: </b> ", input$contrastselection[1], "<br /> <p> <b> B: </b> ", input$contrastselection[2])
    }
    
    if (length(input$contrastselection) == 3) {
      text <- paste0("<p> <b> A: </b> ", input$contrastselection[1], "<br /> <p> <b> B: </b> ", input$contrastselection[2], "<br /> <p> <b> C: </b> ", input$contrastselection[3])
    }
    
    if (length(input$contrastselection) == 4) {
      text <- paste0("<p> <b> A: </b> ", input$contrastselection[1], "<br /> <p> <b> B: </b> ", input$contrastselection[2], "<br /> <p> <b> C: </b> ", input$contrastselection[3], "<br /> <p> <b> D: </b> ", input$contrastselection[4])
    }
   
    if ((length(input$contrastselection) > 4)){
      text <- "More than four comparisons selected. Please select between two and four comparisons."
    }
    
    if (length(input$contrastselection) < 2){
      text <- "Less than two comparisons selected. Please select between two and four comparisons."
    }
    
    return(text)
  })
  
  
  ################################################################################################################################
  #GO analysis
  ################################################################################################################################
  
  
  output$uigoacomp <- renderUI({
    pickerInput(
      inputId = "goacomp",
      label = "Comparison",
      choices = get_contrasts(meta1()),
      options = list(
        style = "btn-primary")
    )
  })
  
  
  affyLib <- reactive({
    req(normData())
    annotation <- annotation(normData())
    
    if (annotation == "hugene10stv1"){
      affyLib <- "hugene10stprobeset.db"
    } else {
      affyLib <- paste(annotation, "db", sep = ".")
    }
    
    #install annotation package from bioconductor
    
  })
  
  sampleGOdata <- eventReactive(input$goa.ok, {
    
    top.table <- top.table()[[input$goacomp]]
    
    if(input$raworadjgoa == "raw"){
      geneList <- ifelse((top.table$P.Value < input$pthresholdgoa) & (abs(top.table$logFC) > input$logFCthresholdgoa), 1, 0)
    }
    
    if(input$raworadjgoa == "adj"){
      geneList <- ifelse((top.table$adj.P.Val < input$pthresholdgoa)  & (abs(top.table$logFC) > input$logFCthresholdgoa), 1, 0)
    }
    
    names(geneList) <- top.table$Probeset.ID
    
    
    sampleGOdata <- new("topGOdata",
                        description = "Simple session", ontology = input$ontology,
                        allGenes = geneList, geneSelectionFun = function(x)(x == 1),
                        nodeSize = 10,
                        annot = annFUN.db, affyLib = affyLib())
    
    return(sampleGOdata)
  })
  
  
  resultFisher <- eventReactive(input$goa.ok, {
    resultFisher <- runTest(sampleGOdata(), algorithm = "classic", statistic = "fisher")
    return(resultFisher)
  })
  
  output$goa <- renderDataTable(NULL)
  
  observeEvent(input$goa.ok, {
    output$goa <- renderDataTable({
      
      allRes <- GenTable(sampleGOdata(), classicFisher = resultFisher(), orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
      return(allRes)
      
    }, options = list(lengthMenu = c(5, 10, 20, 30), pageLength = 5))
    
    
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
      dev.off()
      
      list(src = GOplot,width = WIDTH,height = HEIGHT,alt = "This is alternate text")
    },deleteFile = TRUE)
    
  })
  
  
  
  
}




