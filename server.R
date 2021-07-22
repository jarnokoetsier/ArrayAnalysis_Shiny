

###############################################################################################################################

#SERVER

###############################################################################################################################


server <- function(input, output, session){
  
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
    text = "Welcome to the ArrayAnalysis app! 
    Get started by entering a database accession number or selecting your own CEL files."
  )
  
  
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
  
  
  #Example dataset
  
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
  
  
  #Info panel
  
  observeEvent(input$infopanel1, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "Enter the database accession number. Find your accession number at GEO or ArrayExpress. 
      If you select 'upload CELs', you can also perform the analysis on your own CEL files.",
      type = "info"
    )
  })
  
  
  #Load data from database
  
  database <- eventReactive(input$downloadGEO, {
    input$database
  })
  
  
  accession <- eventReactive(input$downloadGEO, {
    input$getGEO
  })
  
  gset <- eventReactive(input$downloadGEO, {
    if (database() == "GEO"){
      gset <- getGEO(accession(), GSEMatrix =TRUE, getGPL = FALSE)
    }
    if (database() == "ArrayExpress"){
      gset <- ArrayExpress(accession())
    }
    return(gset)
  })
  
  
  
  #Loading message
  
  observeEvent(input$downloadGEO, {
    showModal(modalDialog("This might take a while. Please be patient.", 
                          title = "Dataset is being loaded....", 
                          footer=NULL))
    get_grouping(gset(), database = database())
    removeModal()
  })
  
  
  
  #go to next tab
  
  observeEvent(input$downloadGEO, {
    showTab("navbar", target = "panel2")
  })
  
  
  observeEvent(if ((length(get_grouping(gset(), database = database())) > 0) | (length(get_grouping(gset(), database = database())) > 0)) TRUE,  {
    updateNavbarPage(session, "navbar",
                     selected = "panel2")
    
  })
  
  
  
  
  
  ################################################################################################################################
  #get meta
  ################################################################################################################################
  
  #Data successfully downloaded
  
  observeEvent(if ((length(get_grouping(gset(), database = database())) > 0) | (length(get_grouping(gset(), database = database())) > 0)) TRUE,  {
    sendSweetAlert(
      session = session,
      title = "Success!!",
      text = "Dataset successfully selected! Now you can group the samples.",
      type = "success")
    
  })
  
  
  #Incorrect dataset selected
  
  observeEvent(if ((length(get_grouping(gset(), database = database())) == 0) | (length(get_grouping(gset(), database = database())) == 0)) TRUE, {
    sendSweetAlert(
      session = session,
      title = "Error!",
      text = "No valid dataset selected. Please try again.",
      type = "error")
    
  })
  
  
  #Info panel
  
  observeEvent(input$infopanel2, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "In this panel, you can classify the samples into groups. 
      These groups are required for the differential expression analysis.",
      type = "info"
    )
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
  
  #Make "meta" dataframe
  meta <- reactive({
    
    if (length(input$makegroups) > 0){
      
      if (input$makegroups == "Dataset"){
        
        if (length(input$groupselect) > 0){
          groups <- get_grouping(gset(), database = database())
          meta <- get_meta(gset = gset(), grouping_column = groups[,input$groupselect], database = database())
          return(meta)
        }
      }
      
      if (input$makegroups == "Manual grouping"){
        
        group2 <- cbind(input$owngroupsin, rep("group 2", length(input$owngroupsin)))
        
        group1 <- setdiff(samplelist(), input$owngroupsin)
        group1 <- cbind(group1, rep("group 1", length(group1)))
        
        meta <- as.data.frame(rbind(group1, group2))
        rownames(meta) <- NULL
        colnames(meta) <- c("Sample.ID", "Grouping")
        
        return(meta)
      }
    }
    
  })
  
  
  #make table of meta data
  output$grouping <- renderDataTable({
    
    meta()
    
  })
  
  
  #download table of meta data
  output$downloadmeta <- downloadHandler(
    filename = "meta data",
    content = function(file){
      write.table(meta(), file)
    }
  )
  
  
  observeEvent(input$meta.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel3")
    
    
  })
  
  
  observeEvent(input$meta.ok, {
    showTab("navbar", target = "panel3")
  })
  
  
  ################################################################################################################################
  #Pre-processing
  ################################################################################################################################
  
  #Success message
  
  observeEvent(input$meta.ok, {
    sendSweetAlert(
      session = session,
      title = "success!!",
      text = "Meta data successfully selected! Now you can choose between pre-processing options.",
      type = "success"
    )
  })
  
  
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

  
  outliers <- eventReactive(input$ann.ok, {
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
  
  
  
  #Read CEL files
  data1 <- eventReactive(input$ann.ok, {
    
    accession = accession()
    database = database()
    outliers = outliers()
    
    exdir = paste0("data_", accession)
    
    if (database == "GEO"){
      
      #Download data
      if (!file.exists(exdir)){
        
        getGEOSuppFiles(accession)
        tarfile = paste0(accession, "/", accession, "_RAW.tar")
        untar(tarfile, exdir = exdir)
        
      }
      
      
      celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", full.names = TRUE)
      
      if (!is.null(outliers)){
        
        outliers_select <- outliers[1]
        
        if (length(outliers) > 1){
          for (i in 1:(length(outliers)-1)){
            outliers_select <- paste(outliers_select, outliers[1+i], sep = "|")
          }
          
        }
        
        celfiles = celfiles[-grep(outliers_select, celfiles)]
      }
      
    }
    
    if (database == "ArrayExpress"){
      if (!file.exists(exdir)) {
        getAE(accession, type = "raw", extract = FALSE)
        unzip(paste0(accession, ".raw.1.zip"), exdir = exdir)
      }
      
      
      
      celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", full.names = TRUE)
      
      if (!is.null(outliers)){
        
        outliers_select <- outliers[1]
        
        if (length(outliers) > 1){
          for (i in 1:(length(outliers)-1)){
            outliers_select <- paste(outliers_select, outliers[1+i], sep = "|")
          }
          
        }
        
        celfiles = celfiles[-grep(outliers_select, celfiles)]
      }
      
    }
    
    if (database == "Upload CELs") {
      celfiles = unzip(input$uploadcelsin$datapath)
      
      if (!is.null(outliers)){
        
        outliers_select <- outliers[1]
        
        if (length(outliers) > 1){
          for (i in 1:(length(outliers)-1)){
            outliers_select <- paste(outliers_select, outliers[1+i], sep = "|")
          }
          
        }
        
        celfiles = celfiles[-grep(outliers_select, celfiles)]
      }
    }
    
    data1 = ReadAffy(filenames=celfiles)
    
    return(data1)
    
  })
  
  #RMA settings
  
  output$uibg <- renderUI({
    if(input$rma == "RMA"){
      awesomeCheckbox(
        inputId = "bg",
        label = "Background correction", 
        value = TRUE)
    }
  })
  
  output$uinorm <- renderUI({
    if(input$rma == "RMA"){
      awesomeCheckbox(
        inputId = "norm",
        label = "Quantile normalization", 
        value = TRUE)
    }
  })
  
  
  #Advanced settings
  
  output$uibgcorrectmethod <- renderUI({
    if(input$rma == "Advanced"){
      prettyRadioButtons(
        inputId = "bgcorrectmethod",
        label = "Background correction", 
        choices = c("rma", "mas", "none"),
        selected = "rma",
        status = "primary")
    }
  })
  
  output$uinormmethod <- renderUI({
    if(input$rma == "Advanced"){
      prettyRadioButtons(
        inputId = "normmethod",
        label = "Normalization", 
        choices = c("quantiles", "loess", "constant"),
        selected = "quantiles",
        status = "primary")
    }
  })
  
  output$uipmcorrectmethod <- renderUI({
    if(input$rma == "Advanced"){
      prettyRadioButtons(
        inputId = "pmcorrectmethod",
        label = "PM correction", 
        choices = c("pmonly", "mas", "substractmm"),
        selected = "pmonly",
        status = "primary")
    }
  })
  
  output$uisummarymethod <- renderUI({
    if(input$rma == "Advanced"){
      prettyRadioButtons(
        inputId = "summarymethod",
        label = "Summarization", 
        choices = c("medianpolish", "mas", "avgdiff"),
        selected = "medianpolish",
        status = "primary")
    }
  })
  

  
  #RMA background correction, normalization, and summarization
  data.expr <- eventReactive(input$ann.ok, {
    
    showModal(modalDialog("This might take a while. Please be patient.", 
                          title = "Data is being pre-processed....", 
                          footer=NULL))
    
    req(input$rma)
    
    if(input$rma == "RMA"){
    
      data.norm = rma(data1(), background = input$bg, normalize = input$norm)
    }
    
    if(input$rma == "Advanced"){
      req(input$bgcorrectmethod)
      req(input$normmethod)
      req(input$pmcorrectmethod)
      req(input$summarymethod)
      
      data.norm <- expresso(data1(), 
                            bgcorrect.method = input$bgcorrectmethod,
                            normalize.method = input$normmethod,
                            pmcorrect.method = input$pmcorrectmethod,
                            summary.method = input$summarymethod)
    }
    
    
    #Expression values
    data.expr <- exprs(data.norm)
    data.expr <- data.expr[rownames(data.expr) != "nonsense",]
      
    
    return(data.expr)
  })
  
  
  
  #logData

  
  normlogData <- eventReactive(input$ann.ok, {
    
    #Get raw expression data
    data.expr <- data.expr()
    
    #Get samples
    samples <- colnames(data.expr)
    
    #make dataframe with samples + intensities
    sampleNames = vector()
    logs = vector()
    for (i in 1:length(samples))
    {
      sampleNames = c(sampleNames,rep(samples[i],nrow(data.expr)))
      logs = c(logs,log2(data.expr[,samples[i]]))
    }
    
    logData = data.frame(logInt=logs,sampleName=sampleNames)
    
    
    logData <- logData %>% fuzzyjoin::fuzzy_inner_join(meta(), by = c("sampleName" = "Sample.ID"), match_fun = str_detect)
    
    logData <- arrange(logData, Grouping)
    logData$Grouping <- factor(logData$Grouping, levels = unique(logData$Grouping))
    logData$Sample.ID <- factor(logData$Sample.ID, levels = unique(logData$Sample.ID))
    
    return(logData)
  })
  
  output$normboxplot <- renderPlot(NULL)
  output$normhist <- renderPlot(NULL)
  
  observeEvent(if (length(normlogData()) > 0) TRUE, {
    #Normalized boxplot
    output$normboxplot <- renderPlot({
      
      logData <- normlogData()
      meta <- meta()
      
      #Make colours
      meta1 <- arrange(meta, Grouping)
      meta1$Grouping <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
      meta1$Sample.ID <- factor(meta1$Sample.ID, levels = unique(meta1$Sample.ID))
      
      myPalette <- colorsByFactor(experimentFactor = meta1[,2])
      
      samplecolours <- cbind(meta1, myPalette$plotColors)
      legendcolours <- inner_join(samplecolours, as.data.frame(myPalette$legendColors), by = c("myPalette$plotColors" = "myPalette$legendColors"))
      
      
      #Make Boxplot
      
      ggplot(logData, aes(x = Sample.ID, y = logInt, fill = Sample.ID)) +
        geom_boxplot() +
        scale_fill_manual(breaks = legendcolours[,1], 
                          label = legendcolours[,2],
                          values = samplecolours[,3]) +
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
        guides(fill=guide_legend(title="Grouping"))
      
    })

    
    
    
    
    #Normalized histogram
    output$normhist <- renderPlot({
      
      logData <- normlogData()
      meta <- meta()
      
      #Make colours
      meta1 <- arrange(meta, Grouping)
      meta1$Grouping <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
      meta1$Sample.ID <- factor(meta1$Sample.ID, levels = unique(meta1$Sample.ID))
      
      myPalette <- colorsByFactor(experimentFactor = meta1[,2])
      
      samplecolours <- cbind(meta1, myPalette$plotColors)
      legendcolours <- inner_join(samplecolours, as.data.frame(myPalette$legendColors), by = c("myPalette$plotColors" = "myPalette$legendColors"))
      
      
      #Make density plot
      
      ggplot(logData, aes(x = logInt, colour = Sample.ID)) +
        geom_density(size = 1) +
        scale_colour_manual(breaks = legendcolours[,1], 
                            label = legendcolours[,2],
                            values = samplecolours[,3]) +
        labs(title = "Density plot of normalized intensities",
             subtitle = "Distributions should be comparable between arrays") +
        xlab("Normalized log intensity") +
        ylab("Density") +
        theme_classic() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5),
              legend.background = element_rect(fill="#d3d3d3", colour = "black"),
              legend.title = element_text(face = "bold", hjust = 0.5)) +
        guides(colour=guide_legend(title="Grouping"))
      
    })
  })
  
  #After normalization....
  
  observeEvent(if (length(normlogData()) > 0) TRUE, {
    
    
    #get proceed button
    output$proceedann <- renderUI({
      actionBttn(inputId = "ann.proceed",
                 label = "Next",
                 style = "jelly",
                 color = "primary")
    })
    
    #Download data expr
    output$downloadexpr <- downloadHandler(
      filename = "data expr",
      content = function(file){
        write.table(data.expr(), file, sep = "\t")
      })
    
    removeModal()
    
    sendSweetAlert(
      session = session,
      title = "success!!",
      text = "Data successfully pre-processed! Please wait for the QC plots to be rendered.",
      type = "success"
    )
  })
  
  
  observeEvent(input$ann.proceed, {
    updateNavbarPage(session, "navbar",
                     selected = "panel4")
    
  })
  
  observeEvent(input$ann.proceed, {
    showTab("navbar", target = "panel4")
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
        pca.plot(data.PC(), meta(), pc.x(), pc.y())
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
        pca3d(data.PC(), meta(), x = pc.x(), y = pc.y(), z = pc.z())
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
      title = "success!!",
      text = "Now you can perform various statistical analyses",
      type = "success")
    
  })
    
  observeEvent(input$pca.ok, {
    updateNavbarPage(session, "navbar",
                     selected = "panel5")
    
  })
  
  observeEvent(input$pca.ok, {
    showTab("navbar", target = "panel5")
  })
  
  ################################################################################################################################
  #differential expression analysis
  ################################################################################################################################
  
  
  #get statistics
  
  top.table <- reactive({
    
    top.table <- diff_expr(data.expr(), meta(), comparisons = get_contrasts(meta()))
    
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
      choices = get_contrasts(meta()),
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
      choices = get_contrasts(meta()),
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
  #volcano plot
  ################################################################################################################################
  
  
  output$uivolcanocomp <- renderUI({
    pickerInput(
      inputId = "volcanocomp",
      label = "Comparison",
      choices = get_contrasts(meta()),
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
      choices = get_contrasts(meta()),
      multiple = TRUE
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
        genes <- toptableselection[[i]] %>% 
          filter(P.Value <= 0.05) %>%
          filter(logFC > 1 | logFC < -1)
        
        genelist[[i]] <- genes$Probeset.ID
      }
      
      
      ggvenn(genelist, set_name_size = 4, text_size = 3, stroke_size = 0.7)
      
    }
    
  })
  
}




