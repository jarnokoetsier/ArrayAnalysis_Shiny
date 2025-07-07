#==============================================================================#
# Name: FUNCTIONS.R
# Description: Functions of the ArrayAnalysis Shiny app
#==============================================================================#


#==============================================================================#
# firstup()
#==============================================================================#

# DESCRIPTION:
# Put first letter of character string to capital

# VARIABLES:
# x: vector of character strings

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#==============================================================================#
# getCELS()
#==============================================================================#

# DESCRIPTION:
# Unzip the celfiles and return the path to the unzipped folder

# VARIABLES:
# zippath: path to zip file
# shiny_upload: is it an object from shiny upload

getCELs <- function(zippath, shiny_upload = TRUE){
  tryCatch({
    if (shiny_upload){
      
      # Unzip
      unzip(zippath$datapath, 
            exdir = paste0(tempdir(), "/CELfiles"))
      
      # List CEL files
      celfiles <- list.files(tail(list.dirs(paste0(tempdir(), "/CELfiles")),1),
                             pattern = "CEL", 
                             full.names = TRUE)
    }
    if (!shiny_upload){
      
      # Unzip
      unzip(zippath, 
            exdir = paste0(tempdir(), "/CELfiles"))
      
      # List CEL files
      celfiles <- list.files(tail(list.dirs(paste0(tempdir(), "/CELfiles")),1),
                             pattern = "CEL", 
                             full.names = TRUE)
    }
    return(celfiles)
    
  }, error = function(cond){
    return(NULL)
  })
}

#==============================================================================#
# readCELs()
#==============================================================================#

# DESCRIPTION:
# Read CEL files

# VARIABLES:
# celfiles: vector with path to each of the celfiles
# rm: remove unzipped directory

readCELs <- function(celfiles, zippath, rm = FALSE){
  tryCatch({

    # List CEL files
    celfiles <- list.files(tail(list.dirs(paste0(tempdir(), "/CELfiles")),1),
                           pattern = "CEL", 
                           full.names = TRUE)
    
    if(length(celfiles) < 1){
      # Unzip
      unzip(zippath$datapath, 
            exdir = paste0(tempdir(), "/CELfiles"))
      
      # List CEL files
      celfiles <- list.files(tail(list.dirs(paste0(tempdir(), "/CELfiles")),1),
                             pattern = "CEL", 
                             full.names = TRUE) 
      
    }
    
    
    gxData <- affy::ReadAffy(filenames = celfiles)
    if(rm){unlink(paste0(tempdir(), "/CELfiles"), recursive = TRUE, force = TRUE)}
    return(gxData)
  }, error = function(cond){
    gxData <- oligo::read.celfiles(filenames = celfiles)
    if(rm){unlink(paste0(tempdir(), "/CELfiles"), recursive = TRUE, force = TRUE)}
    return(gxData)
  })
}


#==============================================================================#
# getMetaData()
#==============================================================================#

# DESCRIPTION:
# Read and prepare meta data

# VARIABLES:
# path: path to meta data file (.csv/.tsv or Series Matrix File)
# celfiles: path to celfiles (used to filter meta data)
# filetype: format of meta data file: ".tsv/.csv file" or "Series Matrix File"

getMetaData <- function(path, celfiles, filetype){
  
  if (filetype == ".tsv/.csv file"){
    
    # Read .tsv/.csv file
    if (stringr::str_detect(path, ".tsv")){
      metaData <- as.data.frame(data.table::fread(path, sep = "\t"))
    } else{
      metaData <- as.data.frame(data.table::fread(path))
    }
    
  }
  
  if (filetype == "Series Matrix File"){
    # Read Series Matrix file
    gse <- GEOquery::getGEO(filename=path)
    
    # Extract meta data
    metaData <- Biobase::phenoData(gse)@data
    
    # remove columns with same information for each sample
    remove_col <- NULL
    for (j in 1:ncol(metaData)){
      if ((length(unique(metaData[,j])) == 1)|(min(nchar(metaData[,j]), na.rm = TRUE) > 50)){
        remove_col <- c(remove_col,j)
      }
      
    }
    metaData <- metaData[,-remove_col]
  }
  metaData <- metaData[,!duplicated(colnames(metaData))]
  
  # get sample IDs
  CELsamples <- stringr::str_remove(basename(celfiles),"\\.CEL.*")
  
  # get column with samples IDs
  sumIDs <- rep(0, ncol(metaData))
  for (i in 1:ncol(metaData)){
    if (length(unique(metaData[,i])) == nrow(metaData)){
      sumIDs[i] <- sum(metaData[,i] %in% CELsamples)
    }
  }
  
  # If there are no common sample IDs, add them to the object (using fuzzy matching)
  if(max(sumIDs) == 0){
    CELsamples_alt <- stringr::str_remove(CELsamples, "_.*")
    
    # get column with samples IDs
    sumIDs <- rep(0, ncol(metaData))
    for (i in 1:ncol(metaData)){
      if (length(unique(metaData[,i])) == nrow(metaData)){
        sumIDs[i] <- sum(metaData[,i] %in% CELsamples_alt)
      }
    }
    
    # Add additional column with CEL names
    combineCELs <- data.frame(AltName = CELsamples_alt,
                              CELName = CELsamples)
    
    metaData_copy <- metaData
    colnames(metaData_copy)[which.max(sumIDs)] <- "y"
    metaData <- inner_join(metaData_copy,
                           combineCELs,
                           by = c("y" = "AltName"))
    metaSamples <- metaData[,ncol(metaData)]
    
    # Get common samples
    commonSamples <- intersect(CELsamples, metaSamples)
    
    # Filter meta data for common samples
    metaData_fil <- metaData[metaSamples %in% commonSamples,]
    
    # Remove duplicate samples
    metaData_fil <- metaData_fil[!duplicated(metaData_fil[,ncol(metaData)]),]
    
    # Set sample names as row names
    rownames(metaData_fil) <- metaData_fil[,ncol(metaData)]
    
  } else {
    metaSamples <- metaData[,which.max(sumIDs)]
    
    # Get common samples
    commonSamples <- intersect(CELsamples, metaSamples)
    
    # Filter meta data for common samples
    metaData_fil <- metaData[metaSamples %in% commonSamples,]
    
    # Remove duplicate samples
    metaData_fil <- metaData_fil[!duplicated(metaData_fil[,which.max(sumIDs)]),]
    
    # Set sample names as row names
    rownames(metaData_fil) <- metaData_fil[,which.max(sumIDs)]
  }
  colnames(metaData_fil) <- make.names(colnames(metaData_fil))
  return(metaData_fil)
}

#==============================================================================#
# autoGroup()
#==============================================================================#

# DESCRIPTION:
# Guess the experimental group from the meta data

# VARIABLES:
# metaData: meta data object (dataframe)

autoGroup <- function(metaData){
  tryCatch({
    term_col <- NULL
    len_col <- NULL
    for (i in 1:ncol(metaData)){
      if (length(grep("control|non|healthy|treat|wt|wildtype|wild-type", metaData[,i], ignore.case = TRUE)) > 0){
        term_col <- c(term_col,i)
      }
      
      if (((length(unique(metaData[,i])) > 1) & (length(unique(metaData[,i])) < nrow(metaData)))|(min(nchar(metaData[,i])) > 50)){
        len_col <- c(len_col, i)
      }
    }
    
    
    exp_col <- intersect(term_col, len_col)
    if (is.null(exp_col)){
      exp_col <- len_col
    }
    
    if (length(exp_col)>0){
      return(colnames(metaData)[exp_col[1]])
    } else {
      return(colnames(metaData)[1])
    }
  }, error = function(cond){
    return(colnames(metaData)[1])
  })
  
}

#==============================================================================#
# getOrganism()
#==============================================================================#

# DESCRIPTION:
# Get the organism from data object (AffyBatch or GeneFeaturSet)

# VARIABLES:
# gxData: data object (AffyBatch or GeneFeaturSet)

getOrganism <- function(gxData){
  
  organism <- "Homo sapiens"
  if (stringr::str_detect(gxData@annotation,regex("Cattle|Bos|taurus|Bt",ignore_case = T))){
    organism <- "Bos taurus"
  }
  if (stringr::str_detect(gxData@annotation,regex("nematode|Caenorhabditis|elegans|Ce",ignore_case = T))){
    organism <- "Caenorhabditis elegans"
  }
  if (stringr::str_detect(gxData@annotation,regex("Rat|Rattus|norvegicus|Rn",ignore_case = T))){
    organism <- "Rattus norvegicus"
  }
  if (stringr::str_detect(gxData@annotation,regex("mouse|Mus|musculus|Mm",ignore_case = T))){
    organism <- "Mus musculus"
  }
  if (stringr::str_detect(gxData@annotation,regex("human|Homo|sapiens|Hs",ignore_case = T))){
    organism <- "Homo sapiens"
  }
  return(organism)
}

#==============================================================================#
# uploadcdfenv()
#==============================================================================#

# DESCRIPTION:
# Add custom CDF environment

# VARIABLES:
# Data: data object (AffyBatch)
# cdf_path: path to CDF file

uploadcdfenv <- function(Data,cdf_path){
  
  #initial value
  CDFenv <- 0
  
  # recall which cdfName was added, in case no updated one is found (set back
  # even if it does not exist)
  presetCDF <- Data@cdfName
  
  #try to load cdf file
  try(CDFenv <- makecdfenv::make.cdf.env(filename = basename(cdf_path),cdf.path=dirname(cdf_path)),TRUE)
  
  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("Could not load custom CDF environment for this chip type - object kept as is")
  }
  
  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}

#==============================================================================#
# addUpdatedCDFenv()
#==============================================================================#

# DESCRIPTION:
# Add new CDF environment to the data object (AffyBatch)

# VARIABLES:
# Data: data object (AffyBatch)
# species: species/organism
# type: custom annotation type

addUpdatedCDFenv <- function(Data, species=NULL, type="ENSG") {
  # note: this function will add an updated cdf environment to the data object
  # and will overwrite a possible already loaded environment, unless no updated
  # cdf environment is found
  
  # developer's note: it may be of interest to find out whether available
  # species and types can be retrieved automatically from the brainarray website
  
  # Match the species to two letter symbols  
  spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
           "Sc","Sp","Ss")
  names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
                  "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
                  "Drosophila melanogaster","Gallus gallus","Homo sapiens",
                  "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
                  "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
  
  species <- spp[tolower(names(spp))==tolower(species)]
  
  
  #initial value
  CDFenv <- 0
  
  # recall which cdfName was added, in case no updated one is found (set back
  # even if it does not exist)
  presetCDF <- Data@cdfName
  
  #try to find updated cdf file of choice
  print(Data@cdfName<-paste(Data@annotation,species,type,sep="_"))
  suppressWarnings(try(CDFenv <- affy::getCdfInfo(Data),TRUE))
  #try without a version number
  print(Data@cdfName<-paste(gsub("v[0-9]$","",Data@annotation),species,type,sep="_"))
  suppressWarnings(try(CDFenv <- affy::getCdfInfo(Data),TRUE)) 
  
  #if it hasn't loaded, try to download
  if ((class(CDFenv)!="environment")) {
    install.packages(tolower(paste(Data@annotation,species,type,"cdf",sep="")),
                     repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- affy::getCdfInfo(Data),TRUE))
  }
  
  #if it hasn't loaded, try to download without version number
  if ((class(CDFenv)!="environment")) {
    install.packages(tolower(paste(gsub("v[0-9]$","",Data@annotation),species,type,"cdf",sep="")),
                     repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- affy::getCdfInfo(Data),TRUE))
  }
  
  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("Could not automatically retrieve CDF environment for this chip type - object kept as is")
  }
  
  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}

#==============================================================================#
# microarrayNormalization()
#==============================================================================#

# DESCRIPTION:
# Additional processing of raw microarray data (AffyBatch or GeneFeaturSet)

# VARIABLES:
# gxData: data object (AffyBatch or GeneFeaturSet)
# experimentFactor: factor of experimental group
# normMeth: normalization method ("RMA", "GCRMA", "PLIER")
# CDFtype: CDF type (for custom annotations only)
# species: species/organism (for custom annotations only)
# annotations: "Custom annotations"/"Upload annotation file"/"No annotations"
# perGroup_name: perform normalization per group
# annot_file_datapath: path to CDF annotation file

microarrayNormalization <- function(gxData,
                                    experimentFactor,
                                    normMeth,
                                    CDFtype,
                                    species,
                                    annotations,
                                    perGroup_name,
                                    annot_file_datapath){
  
  if (class(gxData) == "AffyBatch"){
    
    # Set species:
    
    # Match the species to two letter symbols  
    spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
             "Sc","Sp","Ss")
    names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
                    "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
                    "Drosophila melanogaster","Gallus gallus","Homo sapiens",
                    "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
                    "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
    
    species_abbr <- spp[tolower(names(spp))==tolower(species)]
    
    # Pre-processing settings:
    aType <- "PMMM"
    normMeth <- toupper(normMeth)
    ifelse(annotations =="Custom annotations",
           customCDF<-TRUE,customCDF<-FALSE)
    ifelse(annotations =="Upload annotation file",
           uploadCDF<-TRUE,uploadCDF<-FALSE)
    ifelse(perGroup_name== "Use all arrays",
           perGroup<-FALSE,perGroup<-TRUE)
    
    
    
    # If the "customCDF" option is chosen, apply copy the data in order not to 
    # change the original data object
    Data.copy <- gxData
    if(customCDF){
      print ("Changing CDF before pre-processing")
      Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
    }
    
    if(uploadCDF){
      print ("Changing CDF before pre-processing")
      Data.copy <- uploadcdfenv(Data.copy, annot_file_datapath)
    }
    
    print("Pre-processing is running")
    
    # Number of experimental groups
    nGroups <- 1
    if(perGroup) {
      nGroups <- max(1,length(levels(experimentFactor)))
      if(nGroups==1) warning("normalization per group requested,
                                 but no groups indicated in data set")
    }
    
    # If per group normalization required, or a method selected that does not 
    # return an ExpressionSet object, make a model of class ExpressionSet to 
    # paste real values in, use the relatively fast RMA method.
    # Note that binding of ExpressionSet objects is NOT possible
    if((nGroups>1)) {
      normData <- affy::rma(Data.copy)
      Biobase::exprs(normData)[] <- NA
    }
    
    # Perform normalization for each group
    for(group in 1:nGroups) {
      if(nGroups==1) {
        Data.tmp <- Data.copy
      } else {
        Data.tmp <- Data.copy[,experimentFactor==(levels(experimentFactor)[group])]
      }
      switch(normMeth,
             
             # MAS5 method: DOES NOT WORK!
             "MAS5" = {
               normData.tmp <- affy::mas5(Data.tmp)
             },
             
             # GCRMA method:
             "GCRMA" = {
               if(customCDF) {
                 #probe library needed, first try whether this has been intalled, otherwise do so
                 #probeLibrary <- tolower(paste(Data.tmp@annotation,species_abbr,CDFtype,"probe",sep=""))
                 probeLibrary <- tolower(paste(gsub("v[0-9]$","",Data.tmp@annotation),species_abbr,CDFtype,"probe",sep=""))
                 loaded <- suppressWarnings(
                   try(eval(parse("",-1,paste("library(",probeLibrary,")",
                                              sep=""))),TRUE))
                 if(class(loaded)=="try-error") {
                   install.packages(probeLibrary,
                                    repos="http://brainarray.mbni.med.umich.edu/bioc")
                 }
               }
               if(aType == "PMMM") ntype = "fullmodel"
               if(aType == "PMonly") ntype = "affinities" # good results if most of the genes are not expressed
               normData.tmp <- gcrma::gcrma(Data.tmp, type=ntype, fast = FALSE)
             },
             
             # RMA method:
             "RMA" = {
               normData.tmp <- affy::rma(Data.tmp)
             },
             
             # Plier method:
             "PLIER" = {
               if(aType == "PMMM") ntype = "together"
               if(aType == "PMonly") ntype = "pmonly"
               normData.tmp <- plier::justPlier(Data.tmp, normalize=TRUE,
                                                norm.type = ntype)
             }
      )
      
      # Return data if number of groups = 1
      if(nGroups==1) {
        normData <- normData.tmp
        
        # For MAS5 the data needs to be log-transformed first...
        if(normMeth=="MAS5") Biobase::exprs(normData)<-log2(exprs(normData))
        
      } else {
        
        # Paste the normalized expression values of the group into the 
        # ExpressionSet object (only if number of groups > 1)
        try(
          if(normMeth=="MAS5"){
            Biobase::exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- log2(exprs(normData.tmp))
          }else{
            Biobase::exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- Biobase::exprs(normData.tmp)
          },TRUE)
      }
      
      # Remove temporary files
      rm(normData.tmp, Data.tmp)
    }
    
    rm(Data.copy)
    
    # Remove CEL extension from the column names
    colnames(normData) <- stringr::str_remove(colnames(normData),"\\.CEL.*")
    return(normData)
  }
  
  if (class(gxData) == "GeneFeatureSet"){
    
    # Perform RMA normalization using the oligo packages if data object is
    # GeneFeatureSet
    normData <- oligo::rma(gxData)
    colnames(normData) <- stringr::str_remove(colnames(normData),"\\.CEL.*")
    return(normData)
  }
  
}


#==============================================================================#
# microarrayNormalization_processed()
#==============================================================================#

# DESCRIPTION:
# Additional processing of processed data (ExpressionSet)

# VARIABLES:
# gxData: data object (ExpressionSet)
# experimentFactor: factor of experimental group
# transMeth: transformation method (e.g., "None", "Log2-transformation")
# normMeth: normalization method (e.g., "None", "Quantile")
# perGroup_name: perform normalization per group?

microarrayNormalization_processed <- function(gxData,
                                              experimentFactor,
                                              transMeth,
                                              normMeth,
                                              perGroup_name){
  
  # Get expression values from data object
  m <- Biobase::exprs(gxData)
  
  # Perform log2-transformation (if selected)
  if (transMeth == "Log2-transformation"){
    m <- m[rowSums(m<=0)==0,]
    m <- log2(m+1)
  }
  
  # No additional normalization
  if(normMeth == "None"){
    gxData_final <- Biobase::ExpressionSet(assayData = m)
  }
  
  # Quantile normalization
  if(normMeth == "Quantile"){
    
    # Use all array for normalization
    if (perGroup_name == "Use all arrays"){
      m <- limma::normalizeQuantiles(m)
      gxData_final <- Biobase::ExpressionSet(assayData = m)
    } else{
      for (g in levels(experimentFactor)){
        m[,experimentFactor == g] <- limma::normalizeQuantiles(m[,experimentFactor == g])
      }
      # Perform normalization per experimental group
      gxData_final <- Biobase::ExpressionSet(assayData = m)
    }
  }
  return(gxData_final)
}

#==============================================================================#
# colorsByFactor()
#==============================================================================#

# DESCRIPTION:
# Make color palette

# VARIABLES:
# experimentFactor: factor that will be used for colouring

colorsByFactor <- function(experimentFactor) {
  
  #check whether a factor has been provided
  if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")
  
  if(length(levels(experimentFactor))==1) {
    #if there is only one group (or no groups are provided) take equally spread colors over the rainbow palette
    plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
    #set group legend color to white, as there is not a specific group color
    legendColors <- "white"
  } else {
    #compute the number of colors needed for each class
    tab.tmp <- table(experimentFactor)
    
    #set the two extreme colors for each class
    colors.light <- rainbow(length(levels(experimentFactor)),
                            s=1-sapply(tab.tmp,min,5)*.1)
    colors.dark <- rainbow(length(levels(experimentFactor)),
                           v=1-sapply(tab.tmp,min,5)*.14)
    
    #create the colors to plot, and colors for the legend (average one per experimental group)
    plotColors <- NULL
    legendColors <- NULL
    for(l in 1:length(levels(experimentFactor))) {
      colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
      tmpColors <- colorFun(tab.tmp[l])
      plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
      legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
    }
  }
  return(list(plotColors=plotColors,legendColors=legendColors))
  
}

#==============================================================================#
# geneBoxplot()
#==============================================================================#

geneBoxplot <- function(experimentFactor, 
                        normMatrix, 
                        sel_row,
                        legendColors,
                        groupOrder,
                        jitter = 0.1,
                        rnaseq = FALSE,
                        seed = 123){
  
  # Prepare expression data frame
  if(is.null(groupOrder)){groupOrder <- levels(experimentFactor)}
  if(is.null(legendColors)){legendColors <- colorsByFactor(experimentFactor)$legendColors}
  if(is.null(jitter)){jitter <- 0.1}
  plotExpr <- data.frame(
    logExpr = as.numeric(normMatrix[sel_row,]),
    Grouping = factor(experimentFactor, levels = groupOrder)
  )
  geneName <- rownames(normMatrix)[sel_row]
  
  # Make plot
  if (!rnaseq){
    p <- ggplot2::ggplot(plotExpr) +
      ggplot2::geom_boxplot(aes(x = Grouping, y = logExpr,
                                fill = Grouping), alpha = 0.1,
                            outlier.alpha = 0) +
      ggplot2::geom_point(aes(x = Grouping, y = logExpr, fill = Grouping),
                          position = ggplot2::position_jitter(width = jitter, seed = seed),
                           size = 4, shape = 21, color = "white") +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(expression(log[2] ~" intensity")) +
      ggplot2::ggtitle(geneName) +
      ggplot2::scale_color_manual(values = legendColors) +
      ggplot2::scale_fill_manual(values = legendColors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none",
                     axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 16, face = "bold"),
                     axis.title.y = element_text(size = 16, face = "bold"),
                     plot.title = element_text(size = 20, face = "bold", hjust= 0.5))
  }
  if (rnaseq){
    p <- ggplot2::ggplot(plotExpr) +
      ggplot2::geom_boxplot(aes(x = Grouping, y = logExpr,
                                fill = Grouping), alpha = 0.1,
                            outlier.alpha = 0) +
      ggplot2::geom_point(aes(x = Grouping, y = logExpr, fill = Grouping), 
                           position = ggplot2::position_jitter(width = jitter, seed = seed),
                           size = 4, shape = 21, color = "white") +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(expression(log[2] ~" expression")) +
      ggplot2::ggtitle(geneName) +
      ggplot2::scale_color_manual(values = legendColors) +
      ggplot2::scale_fill_manual(values = legendColors) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none",
                     axis.text.y = element_text(size = 12),
                     axis.text.x = element_text(size = 16, face = "bold"),
                     axis.title.y = element_text(size = 16, face = "bold"),
                     plot.title = element_text(size = 20, face = "bold", hjust= 0.5))
  }
  
  return(p)
  
}


#==============================================================================#
# getBoxplots()
#==============================================================================#

# DESCRIPTION:
# Make static boxplots

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normData: normalized expression data (ExpressionSet)

getBoxplots <- function(experimentFactor, 
                        legendColors,
                        normData, 
                        RNASeq = FALSE,
                        width = 1000,
                        height = 1414,
                        showAlways = TRUE){
  
  plotColors <- colorsByFactor(experimentFactor)$plotColors
  names(legendColors) <- levels(experimentFactor)
  
  for (i in 1:length(legendColors)){
    group_length <- sum(experimentFactor == levels(experimentFactor)[i])
    if (group_length > 1){
    colors <-  c(colorspace::lighten(legendColors[i], amount = rev(seq(0,0.5,length.out = round(group_length/2)))),
                 colorspace::darken(legendColors[i], amount = seq(0,0.5,length.out = group_length - round(group_length/2)+1)[-1])
    )
    } else{
      colors <- legendColors[i]
    }
    
    plotColors[experimentFactor == levels(experimentFactor)[i]] <- colors
  }
  
  if (isTRUE(RNASeq)){
    names(plotColors) <- colnames(normData)
  }
  if (!isTRUE(RNASeq)){
    names(plotColors) <- sampleNames(normData)
  }
  
  
  # Width of all the plots
  WIDTH <- ifelse(is.null(width), 1000, width)
  
  # Height of all the plots
  HEIGHT <- ifelse(is.null(height), 1414, height)
  
  # Point sizes for the plots
  POINTSIZE <- 24
  
  # Maximum number of arrays that can be computed
  MAXARRAY <- 41
  
  # Plot symbols
  plotSymbols <- 18-as.numeric(experimentFactor)
  
  # Legend symbols
  legendSymbols <- sort(plotSymbols, decreasing=TRUE)
  
  Type <- "Norm"
  
  if (!isTRUE(RNASeq)){
    tmain <- " "
    tmtext2 <- "Normalized log intensity\n\n\n"
    description <- NULL
    samples <- sampleNames(normData)
  }
  if (isTRUE(RNASeq)){
    tmain <- NULL
    description <- NULL
    tmtext2 <- "Normalized log counts\n\n\n"
    samples <- colnames(normData)
  }
  
  
  tryCatch({
    DataBoxplot<- tempfile(fileext='.png')
    png(file = DataBoxplot,width=WIDTH,height=HEIGHT, 
        pointsize=POINTSIZE)
    par(mar=c(5.1,4.1,0,2.1))
    par(oma=c(17,0,0,0), cex.axis=1)
    if (class(normData)[[1]] != "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, col=plotColors ,main=tmain,
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if (class(normData)[[1]] == "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, target = "core", col=plotColors,
                               main=" ", 
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if(length(levels(experimentFactor))>1){
      legend("topright", levels(experimentFactor),
             col=legendColors,fill=legendColors, cex = 0.7, bg = "white",
             bty = "o")
    }
    if(length(samples)<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }
    axis(1,at=1:length(samples),las=2,
         labels=samples, cex.axis=cexval)
    axis(2, cex.axis=0.7)
    mtext(tmtext2, side=2, cex=0.8)
    mtext(description, side=3,
          font=1, cex=0.7)
    dev.off()
    
    list(src = DataBoxplot,width = WIDTH,height = HEIGHT,
         alt = "This is alternate text")
    
    # set height and width in case of error
  }, error = function(cond){
    DataBoxplot<- tempfile(fileext='.png')
    png(file = DataBoxplot,width=1000,height=1414, 
        pointsize=POINTSIZE)
    par(mar=c(5.1,4.1,0,2.1))
    par(oma=c(17,0,0,0), cex.axis=1)
    if (class(normData)[[1]] != "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, col=plotColors ,main=tmain,
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if (class(normData)[[1]] == "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, target = "core", col=plotColors,
                               main=tmain, 
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if(length(levels(experimentFactor))>1){
      legend("topright", levels(experimentFactor),
             col=legendColors,fill=legendColors, cex = 0.7, bg = "white",
             bty = "o")
    }
    if(length(samples)<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }
    axis(1,at=1:length(samples),las=2,
         labels=samples, cex.axis=cexval)
    axis(2, cex.axis=0.7)
    mtext(tmtext2, side=2, cex=0.8)
    mtext(description, side=3,
          font=1, cex=0.7)
    dev.off()
    
    list(src = DataBoxplot,width = 1000,height = 1414,
         alt = "This is alternate text")
  })
}


#==============================================================================#
# getBoxplots_download()
#==============================================================================#

# DESCRIPTION:
# Make static boxplots

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normData: normalized expression data (ExpressionSet)

getBoxplots_download <- function(experimentFactor, 
                                 legendColors,
                                 normData, 
                                 RNASeq = FALSE,
                                 width = 1000,
                                 height = 1414){
  
  plotColors <- colorsByFactor(experimentFactor)$plotColors
  names(legendColors) <- levels(experimentFactor)
  
  for (i in 1:length(legendColors)){
    group_length <- sum(experimentFactor == levels(experimentFactor)[i])
    colors <-  c(colorspace::lighten(legendColors[i], amount = rev(seq(0,0.5,length.out = round(group_length/2)))),
                 colorspace::darken(legendColors[i], amount = seq(0,0.5,length.out = group_length - round(group_length/2)+1)[-1])
    )
    
    plotColors[experimentFactor == levels(experimentFactor)[i]] <- colors
  }
  
  if (isTRUE(RNASeq)){
    names(plotColors) <- colnames(normData)
  }
  if (!isTRUE(RNASeq)){
    names(plotColors) <- sampleNames(normData)
  }
  
  
  # Maximum number of arrays that can be computed
  MAXARRAY <- 41
  
  # Plot symbols
  plotSymbols <- 18-as.numeric(experimentFactor)
  
  # Legend symbols
  legendSymbols <- sort(plotSymbols, decreasing=TRUE)
  
  Type <- "Norm"
  
  if (!isTRUE(RNASeq)){
    tmain <- "Boxplot of normalized intensities"
    tmtext2 <- "Normalized log intensity\n\n\n"
    description <- "Distributions should be comparable between arrays\n"
    samples <- sampleNames(normData)
  }
  if (isTRUE(RNASeq)){
    tmain <- "Boxplot of normalized counts"
    tmtext2 <- "Normalized log counts\n\n\n"
    description <- "Distributions should be comparable between samples\n"
    samples <- colnames(normData)
  }
    #DataBoxplot<- tempfile(fileext='.png')
    #png(file = DataBoxplot,width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(mar=c(5.1,4.1,4.1,2.1))
    par(oma=c(17,0,0,0), cex.axis=1)
    if (class(normData)[[1]] != "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, col=plotColors ,main=tmain,
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if (class(normData)[[1]] == "GeneFeatureSet"){
      suppressWarnings(boxplot(normData, target = "core", col=plotColors,
                               main=tmain, 
                               axes=FALSE, pch = 20, cex=0.7))
    }
    if(length(levels(experimentFactor))>1){
      legend("topright", levels(experimentFactor),
             col=legendColors,fill=legendColors, cex = 0.7, bg = "white",
             bty = "o")
    }
    if(length(samples)<MAXARRAY){
      cexval <- 0.65
    }else{
      cexval <- 0.45
    }
    axis(1,at=1:length(samples),las=2,
         labels=samples, cex.axis=cexval)
    axis(2, cex.axis=0.7)
    mtext(tmtext2, side=2, cex=0.8)
    mtext(description, side=3,
          font=1, cex=0.7)
    #dev.off()
}

#==============================================================================#
# getDensityplots()
#==============================================================================#

# DESCRIPTION:
# Make interactive density plot

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normMatrix: normalized expression matrix

getDensityplots <- function(experimentFactor,
                            legendColors,
                            normMatrix, RNASeq = FALSE){
  
  # Prepare dataframe
  plot_df <- tidyr::gather(as.data.frame(normMatrix))
  plot_df$GeneID <- rep(rownames(normMatrix),ncol(normMatrix))
  plot_df$Group <- rep(experimentFactor, each = nrow(normMatrix))
  
  
  plotColors <- colorsByFactor(experimentFactor)$plotColors
  names(legendColors) <- levels(experimentFactor)
  
  for (i in 1:length(legendColors)){
    group_length <- sum(experimentFactor == levels(experimentFactor)[i])
    colors <-  c(colorspace::lighten(legendColors[i], amount = rev(seq(0,0.5,length.out = round(group_length/2)))),
                 colorspace::darken(legendColors[i], amount = seq(0,0.5,length.out = group_length - round(group_length/2)+1)[-1])
                 )
    
    plotColors[experimentFactor == levels(experimentFactor)[i]] <- colors
  }
  names(plotColors) <- colnames(normMatrix)

  
  if (!isTRUE(RNASeq)){
    xaxis_name <- "Normalized log<sub>2</sub> intensity"
  }
  if (isTRUE(RNASeq)){
    xaxis_name <- "Normalized log<sub>2</sub> counts"
  }
  
  
  #Make density plot
  densityPlot <- ggplot2::ggplot(plot_df, 
                                 ggplot2::aes(x = value, colour = key, shape = Group)) +
    ggplot2::geom_density(size = 1) +
    ggplot2::scale_colour_manual(values = plotColors) +
    ggplot2::geom_hline(yintercept = 0, color = "darkgrey", linewidth = 1.1) +
    #ggplot2::geom_vline(xintercept = 0, color = "darkgrey", linewidth = 1.1) +
    ggplot2::xlab(xaxis_name) +
    ggplot2::ylab("Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::guides(shape=guide_legend(title="Group",
                                       override.aes = list(
                                         colour = legendColors)),
                    colour = "none")
  
  p <-  plotly::ggplotly(densityPlot) %>% 
    layout(xaxis = list(   
             title=xaxis_name),   
           yaxis = list(   
             title='Density')) %>%
    config(
      toImageButtonOptions = list(
        format = "png",
        filename = "DensityPlot"
      )
    )

  return(p)
  
}

#==============================================================================#
# getDensityplots_static()
#==============================================================================#

# DESCRIPTION:
# Make interactive density plot

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normMatrix: normalized expression matrix

getDensityplots_static <- function(experimentFactor, legendColors,
                                   normMatrix, RNASeq = FALSE){
  
  # Prepare dataframe
  plot_df <- tidyr::gather(as.data.frame(normMatrix))
  plot_df$GeneID <- rep(rownames(normMatrix),ncol(normMatrix))
  plot_df$Group <- rep(experimentFactor, each = nrow(normMatrix))
  
  
  plotColors <- colorsByFactor(experimentFactor)$plotColors
  names(legendColors) <- levels(experimentFactor)
  
  for (i in 1:length(legendColors)){
    group_length <- sum(experimentFactor == levels(experimentFactor)[i])
    colors <-  c(colorspace::lighten(legendColors[i], amount = rev(seq(0,0.5,length.out = round(group_length/2)))),
                 colorspace::darken(legendColors[i], amount = seq(0,0.5,length.out = group_length - round(group_length/2)+1)[-1])
    )
    
    plotColors[experimentFactor == levels(experimentFactor)[i]] <- colors
  }
  names(plotColors) <- colnames(normMatrix)
  
  
  if (!isTRUE(RNASeq)){
    xaxis_name <- expression("Normalized"~ log[2] ~ "intensity")
  }
  if (isTRUE(RNASeq)){
    xaxis_name <- expression("Normalized"~ log[2] ~ "counts")
  }
  
  
  #Make density plot
  densityPlot <- ggplot2::ggplot(plot_df, 
                                 ggplot2::aes(x = value, colour = key, shape = Group)) +
    ggplot2::geom_density(size = 1) +
    ggplot2::scale_colour_manual(values = plotColors) +
    ggplot2::geom_hline(yintercept = 0, color = "darkgrey", linewidth = 1.1) +
    #ggplot2::geom_vline(xintercept = 0, color = "darkgrey", linewidth = 1.1) +
    ggplot2::xlab(xaxis_name) +
    ggplot2::ylab("Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank()) +
    ggplot2::guides(shape=guide_legend(title="Group",
                                       override.aes = list(
                                         colour = legendColors)),
                    colour = "none")
  
  return(densityPlot)
  
}

#==============================================================================#
# getHeatmap()
#==============================================================================#

# DESCRIPTION:
# Make interactive heatmap of sample-sample correlations

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normMatrix: normalized expression matrix
# clusterOption1: distance method (e.g., spearman, pearson, euclidean)
# clusterOption2: linkage method (e.g. ward.D2)
# theme: color theme of heatmap
getHeatmap <- function(experimentFactor,
                       legendColors,
                       normMatrix,
                       clusterOption1,
                       clusterOption2,
                       theme){
  
  # Make color palette
  #myPalette <- colorsByFactor(experimentFactor)
  
  # Legend colors
  #legendColors <- myPalette$legendColors
  #names(legendColors) <- levels(experimentFactor)
  
  #note: for computing array correlation, euclidean would not make sense
  #only use euclidean distance to compute the similarity of the correlation
  #vectors for the arrays
  COpt1 <- "pearson"
  if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
  crp <- cor(normMatrix, use="complete.obs", method=COpt1)
  
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
  sidecolors <- data.frame(experimentFactor)
  colnames(sidecolors) <- "Experimental group"
  
  # Select heatmap theme colour
  if (theme == "Default"){
    gradient = viridis(n = 256, alpha = 1, begin = 0, end = 1,
                       option = "viridis")
  }
  
  if (theme == "Yellow-red"){
    gradient = heat.colors(100)
  }
  
  if (theme == "Blues"){
    gradient = RColorBrewer::brewer.pal(9, "Blues")
  }
  if (theme == "Reds"){
    gradient = RColorBrewer::brewer.pal(9, "Reds")
  }
  
  # Make heatmap
  p <- heatmaply::heatmaply(crp, plot_method = "plotly", distfun = my.dist,
                 hclustfun = my.hclust, symm = TRUE, seriate = "mean",
                 titleX = FALSE, titleY = FALSE, key.title = NULL,
                 show_dendrogram = c(TRUE, FALSE), col_side_colors = sidecolors,
                 col_side_palette = legendColors, column_text_angle = 90,
                 colors = gradient)
  
  return(p)
}


#==============================================================================#
# getHeatmap_static()
#==============================================================================#

# DESCRIPTION:
# Make interactive heatmap of sample-sample correlations

# VARIABLES:
# experimentFactor: factor with experimental groups (will be used for colouring)
# normMatrix: normalized expression matrix
# clusterOption1: distance method (e.g., spearman, pearson, euclidean)
# clusterOption2: linkage method (e.g. ward.D2)
# theme: color theme of heatmap
getHeatmap_static <- function(experimentFactor,
                       legendColors,
                       normMatrix,
                       clusterOption1,
                       clusterOption2,
                       theme,
                       width,
                       height,
                       filetype,
                       file){
  
  reate_dend <- function(x, seriate, distfun, hclustfun, na.rm) {
    switch(seriate,
           "mean" = rowMeans(x, na.rm = na.rm),
           "none" = 1:nrow(x),
           "OLO" = {
             dist <- distfun(x)
             hc <- hclustfun(dist)
             dend <- as.dendrogram(hc)
             dend <- seriate_dendrogram(dend, dist, method = "OLO")
             dend
           },
           "GW" = {
             dist <- distfun(x)
             hc <- hclustfun(dist)
             dend <- as.dendrogram(hc)
             dend <- seriate_dendrogram(dend, dist, method = "GW")
             dend
           }
    )
  }
  
  #note: for computing array correlation, euclidean would not make sense
  #only use euclidean distance to compute the similarity of the correlation
  #vectors for the arrays
  COpt1 <- "pearson"
  if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
  crp <- cor(normMatrix, use="complete.obs", method=COpt1)
  
  minVal <- floor(min(crp)*100)/100
  maxVal <- ceiling(max(crp)*100)/100
  
  switch(tolower(clusterOption1),
         "pearson" = {
           my.dist <- function(x) bioDist::cor.dist(x, abs=FALSE)
         },
         "spearman" = {
           my.dist <- function(x) bioDist::spearman.dist(x, abs=FALSE)
         },
         "euclidean" = {
           my.dist <- function(x) bioDist::euc(x)
         }
  )
  
  # Perform clustering
  my.hclust <- function(d) hclust(d, method=clusterOption2)
  sidecolors <- data.frame(experimentFactor)
  colnames(sidecolors) <- "Experimental group"
  
  # Select heatmap theme colour
  if (theme == "Default"){
    gradient <- circlize::colorRamp2(
      c(minVal, mean(c(minVal, maxVal)), maxVal),
      viridis::viridis(n = 3, alpha = 1, begin = 0, end = 1,
                       option = "viridis") # 256
    ) 
  }
  if (theme == "Yellow-red"){
    gradient <- circlize::colorRamp2(
      c(minVal, mean(c(minVal, maxVal)), maxVal),
      heat.colors(3) # #100
    )
  }
  
  if (theme == "Blues"){
    gradient <- circlize::colorRamp2(
      c(minVal, mean(c(minVal, maxVal)), maxVal),
      RColorBrewer::brewer.pal(9, "Blues")[c(1,5,9)] #9
    )
  }
  if (theme == "Reds"){
    gradient <- circlize::colorRamp2(
      c(minVal, mean(c(minVal, maxVal)), maxVal),
      RColorBrewer::brewer.pal(9, "Reds")[c(1,5,9)] #9
    )
    
  }
  
  # Make dendrogram
  row_dend <- reorder(as.dendrogram(my.hclust(my.dist(crp))), rowMeans(crp, na.rm = TRUE))
  #row_dend <- as.dendrogram(my.hclust(my.dist(crp)))

  # Create annotation for columns (top)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    `Experimental group` = experimentFactor,
    col = list(`Experimental group` = legendColors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = grid::unit(1, "cm")
  )
  
  # Custom legends
  lgd2 <- ComplexHeatmap::Legend(title = "Correlation", 
                                 at = c(minVal, mean(c(minVal, maxVal)), maxVal),
                                 col_fun = gradient,
                                 legend_height = grid::unit(0.2, "npc"),
                                 legend_width = grid::unit(0.1, "npc"),
                                 grid_height = grid::unit(0.04, "npc"), 
                                 grid_width = grid::unit(0.16,"npc"),
                                 labels_gp = grid::gpar(fontsize = 16),     
                                 title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
                                 title_gap = grid::unit(5, "mm"))
  lgd1 <- ComplexHeatmap::Legend(title = "Experimental group",
                                 at = names(legendColors),
                                 legend_gp = grid::gpar(fill = legendColors),
                                 legend_height = grid::unit(0.4, "npc"),
                                 legend_width = grid::unit(0.1, "npc"),
                                 grid_height = grid::unit(0.04, "npc"),  
                                 grid_width = grid::unit(0.16,"npc"),
                                 labels_gp = grid::gpar(fontsize = 16),  
                                 title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
                                 title_gap = grid::unit(5, "mm"))
  packed_lgd <- ComplexHeatmap::packLegend(lgd1, lgd2, direction = "vertical",
                                           gap = grid::unit(10, "mm"))
  
  if (filetype == "PNG"){
    png(file,
        width = width * 3 + 3000,
        height = height * 3,
        pointsize = 24,
        res = 300)
  }
  if (filetype == "PDF"){
    pdf(file,
        width = width/160 + 6.25,
        height = height/160)
  }
  if (filetype == "TIF"){
    tiff(file,
        width = width * 3 + 3000,
        height = height * 3,
        pointsize = 24,
        res = 300)
  }
  
  # Create a layout with 2 columns: one for the heatmap, one for the legend
  grid::grid.newpage()
  grid::pushViewport(viewport(layout = grid::grid.layout(nrow = 1, ncol = 2, 
                                                         widths = grid::unit.c(grid::unit(1, "npc") - unit(0.4, "npc"), 
                                                                               grid::unit(0.2, "npc")))))
  
  # Plot the heatmap in the first column
  grid::pushViewport(grid::viewport(layout.pos.col = 1))
  ht <- ComplexHeatmap::Heatmap(
    crp,
    cluster_rows = rev(row_dend),
    cluster_columns = row_dend,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_dend_side = "right",
    top_annotation = ha,
    col = gradient,
    show_heatmap_legend = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    column_names_rot = 90,
    row_dend_width = grid::unit(0.2, "npc"),
    row_dend_gp = grid::gpar(lwd = 3)
  )
  ComplexHeatmap::draw(ht, newpage = FALSE)
  grid::popViewport()
  
  # Plot the packed legend in the second column
  grid::pushViewport(grid::viewport(layout.pos.col = 2))
  ComplexHeatmap::draw(packed_lgd, just = "left")
  grid::popViewport()
  
  dev.off()
}

#==============================================================================#
# plot_PCA()
#==============================================================================#

# DESCRIPTION:
# Make interactive PCA plot

# VARIABLES:
# PCA_data: PCA object retrieved from the prcomp function
# colorFactor: A factor the color the samples in the plot
# xpc: Principal component on horizontal axis
# ypc: Principal component on vertical axis
# zpc: Principal component on z-axis (for 3D plot only)
# xyz: Plot in 3D

plot_PCA <- function(PC_data, colorFactor, legendColors, xpc = 1, ypc = 2, zpc = 3, xyz = FALSE) {
  
  #get PCs
  PCA_df <- data.frame(x = as.data.frame(PC_data$x)[,xpc],
                       y = as.data.frame(PC_data$x)[,ypc],
                       z = as.data.frame(PC_data$x)[,zpc],
                       Group = colorFactor,
                       sampleID = rownames(PC_data$x))
  
  # # Make color palette
  # myPalette <- colorsByFactor(colorFactor)
  # 
  # 
  # # Legend colors
  # legendColors <- myPalette$legendColors
  # names(legendColors) <- levels(colorFactor)
  
  #calculate explained variance
  perc_expl <- round(((PC_data$sdev^2)/sum(PC_data$sdev^2))*100,2)
  
  # Make 2D plot
  if (!isTRUE(xyz)){
    if (length(unique(PCA_df$Group)) < 4){
      pca2d <- ggplot2::ggplot(data = PCA_df, 
                               ggplot2::aes(x = x, y = y, colour = Group, shape = Group, 
                                            text = paste("Sample:", sampleID))) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_colour_manual(values = legendColors) +
        ggplot2::xlab(paste0("PC",xpc, " (", perc_expl[xpc], "%)")) +
        ggplot2::ylab(paste0("PC",ypc, " (", perc_expl[ypc], "%)")) +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.title = ggplot2::element_text(face = "bold", size = 12),
                       legend.title = ggplot2::element_blank())
      
      p <- plotly::ggplotly(pca2d, tooltip = c("x", "y", "colour", "text"))
    } else{
      pca2d <- ggplot2::ggplot(data = PCA_df, 
                               ggplot2::aes(x = x, y = y, colour = Group, 
                                            text = paste("Sample:", sampleID))) +
        ggplot2::geom_point(size = 2) +
        ggplot2::scale_colour_manual(values = legendColors) +
        ggplot2::xlab(paste0("PC",xpc, " (", perc_expl[xpc], "%)")) +
        ggplot2::ylab(paste0("PC",ypc, " (", perc_expl[ypc], "%)")) +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5),
                       axis.title = ggplot2::element_text(face = "bold", size = 12),
                       legend.title = ggplot2::element_blank())
      
      p <- plotly::ggplotly(pca2d, tooltip = c("x", "y", "colour", "text"))
    }

    
    return(p)
  }
  # Make 3D plot
  if (isTRUE(xyz)){
    pca3d <- plotly::plot_ly(x=PCA_df$x, 
                             y=PCA_df$y, 
                             z=PCA_df$z, 
                             type="scatter3d", 
                             mode="markers", 
                             color=PCA_df$Group,
                             colors = legendColors,
                             text = paste('Sample:', PCA_df$sampleID, '<br>Group:', PCA_df$color))
    
    pca3d <- pca3d %>% plotly::layout(scene = list(xaxis = list(title = paste0('PC', xpc, " (", perc_expl[xpc], "%)")),
                                                   yaxis = list(title = paste0('PC', ypc, " (", perc_expl[ypc], "%)")),
                                                   zaxis = list(title = paste0('PC', zpc, " (", perc_expl[zpc], "%)"))))
    p <- pca3d #%>% 
      #plotly::layout(height = 600, width = 800)
    
    return(p)
  }
  
}

#==============================================================================#
# plot_PCA_static()
#==============================================================================#

# DESCRIPTION:
# Make interactive PCA plot

# VARIABLES:
# PCA_data: PCA object retrieved from the prcomp function
# colorFactor: A factor the color the samples in the plot
# xpc: Principal component on horizontal axis
# ypc: Principal component on vertical axis
# zpc: Principal component on z-axis (for 3D plot only)
# xyz: Plot in 3D

plot_PCA_static <- function(PC_data, colorFactor, legendColors, xpc = 1, ypc = 2) {
  
  #get PCs
  PCA_df <- data.frame(x = as.data.frame(PC_data$x)[,xpc],
                       y = as.data.frame(PC_data$x)[,ypc],
                       Group = colorFactor,
                       sampleID = rownames(PC_data$x))
  
  # # Make color palette
  # myPalette <- colorsByFactor(colorFactor)
  # 
  # 
  # # Legend colors
  # legendColors <- myPalette$legendColors
  # names(legendColors) <- levels(colorFactor)
  
  #calculate explained variance
  perc_expl <- round(((PC_data$sdev^2)/sum(PC_data$sdev^2))*100,2)
  
  # Make 2D plot
  
  if (length(unique(PCA_df$Group)) < 4){
    p <- ggplot2::ggplot(data = PCA_df, 
                             ggplot2::aes(x = x, y = y, colour = Group, shape = Group, 
                                          text = paste("Sample:", sampleID))) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = legendColors) +
      ggplot2::xlab(paste0("PC",xpc, " (", perc_expl[xpc], "%)")) +
      ggplot2::ylab(paste0("PC",ypc, " (", perc_expl[ypc], "%)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),
                     axis.title = ggplot2::element_text(face = "bold", size = 12),
                     legend.title = ggplot2::element_blank())
    
  } else{
    p <- ggplot2::ggplot(data = PCA_df, 
                             ggplot2::aes(x = x, y = y, colour = Group, 
                                          text = paste("Sample:", sampleID))) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = legendColors) +
      ggplot2::xlab(paste0("PC",xpc, " (", perc_expl[xpc], "%)")) +
      ggplot2::ylab(paste0("PC",ypc, " (", perc_expl[ypc], "%)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),
                     axis.title = ggplot2::element_text(face = "bold", size = 12),
                     legend.title = ggplot2::element_blank())
    
  }

  return(p)
}

#==============================================================================#
# makeComparisons()
#==============================================================================#

# DESCRIPTION:
# Make all possible comparisons from the experimental factor levels

# VARIABLES:
# ExpLevels: experimental factor levels

makeComparisons <- function(ExpLevels){
  comparisons1 <- NULL
  comparisons2 <- NULL
  for (m in 1:(length(ExpLevels)-1)){
    for (n in (m+1):length(ExpLevels)) {
      comparisons1 <- c(comparisons1, paste(ExpLevels[n], 
                                            ExpLevels[m], sep = " - "))
      comparisons2 <- c(comparisons2, paste(ExpLevels[m], 
                                            ExpLevels[n], sep = " - "))
      
    }
  }
  comparisons <- unique(c(comparisons1,comparisons2))
  return(comparisons)
}

#==============================================================================#
# selFilter()
#==============================================================================#

# DESCRIPTION:
# Select "filter" value for biomaRt annotations

# VARIABLES:
# ProbeAnnotation: selected probe annotation

selFilter <- function(ProbeAnnotation){
  
  selFilter <- "Entrez Gene ID"
  if (!is.null(ProbeAnnotation)){
    if (ProbeAnnotation == "ENSG"){
      selFilter <- "Ensembl Gene ID"
    }
  }
  return(selFilter)
}

#==============================================================================#
# getStatistics()
#==============================================================================#

# DESCRIPTION:
# Perform statistical analysis (limma) and optionally probeset annotation (biomaRt)

# VARIABLES:
# normMatrix: normalized intensity matrix
# metaData: meta data
# expFactor: column name(s) of the meta data that contain the experimental factor
# covGroups_num: continuous/numerical covariates
# covGroups_char: discrete/character covariates
# addAnnotation: add biomaRt annotations to the output table
# biomart_dataset: biomaRt dataset (e.g., hsapiens_gene_ensembl)
# biomart_attributes: biomaRt attributes (output IDs)
# biomart_filters: biomaRt filters (input IDs)

getStatistics <- function(normMatrix, 
                          metaData, 
                          expFactor, 
                          covGroups_num, 
                          covGroups_char, 
                          comparisons,
                          addAnnotation = TRUE,
                          biomart_dataset = "hsapiens_gene_ensembl",
                          biomart_attributes = c("Ensembl Gene ID",
                                                 "Entrez Gene ID",
                                                 "Gene Symbol/Name"),
                          biomart_filters = "Entrez Gene ID"){
  tryCatch({
    dataset <- NULL
    if (!is.null(biomart_dataset)){
      # Replace name of biomaRt filter
      if (biomart_filters %in% c("Ensembl Gene ID",
                                 "Entrez Gene ID",
                                 "Gene Symbol/Name")){
        biomart_filters <- tryCatch({
          switch(biomart_filters,
                 "Gene Symbol/Name" = "gene_name",
                 "Entrez Gene ID" = "entrezgene_id",
                 "Ensembl Gene ID" = "ensembl_gene_id",
          )
        }, error = function(cond){
          return(biomart_filters)
        })
      }
      
      # Replace name(s) of biomaRt attributes
      for (a in 1:length(biomart_attributes)){
        if (biomart_attributes[a] %in% c("Ensembl Gene ID",
                                         "Entrez Gene ID",
                                         "Gene Symbol/Name")){
          biomart_attributes[a] <- switch(biomart_attributes[a],
                                          "Gene Symbol/Name" = "gene_name",
                                          "Entrez Gene ID" = "entrezgene_id",
                                          "Ensembl Gene ID" = "ensembl_gene_id",
          )
        }
      }
    }
    
    
    # Get experiment factor
    if(length(expFactor) > 1){
      experimentFactor <- factor(apply(metaData[,expFactor], 1, paste, collapse = "_" ))
    } else{
      experimentFactor <- factor(metaData[,expFactor])
    }
    
    # Get covariates
    for (n in covGroups_num){
      metaData[,n] <- as.numeric(metaData[,n])
    }
    for (c in covGroups_char){
      metaData[,c] <- as.character(metaData[,c])
    }
    covariates <- c(covGroups_num, covGroups_char)
    
    # Make formula
    if (!is.null(covariates)){
      formula <- paste("~ 0", "experimentFactor", paste0("metaData$`",covariates,"`", collapse = "+"),sep = "+")
    } else {
      formula <- paste("~ 0", "experimentFactor",sep = "+")
    }
    
    # Make design matrix
    design <- model.matrix(as.formula(formula))
    colnames(design)[1:length(levels(experimentFactor))] <- make.names(levels(factor(experimentFactor)))
    colnames(design) <- make.names(colnames(design))
    
    #fit linear model
    data.fit <- limma::lmFit(normMatrix,design)
    
    # Perform statistical comparison for each of the selected comparison
    top_table <- list()
    for (i in comparisons){
      
      # Make contrast
      contrast.matrix <- limma::makeContrasts(contrasts = i,levels=design)
      
      # Fit contrasts
      data.fit.con <- limma::contrasts.fit(data.fit,contrast.matrix)
      
      # eBayes
      data.fit.eb <- limma::eBayes(data.fit.con)
      
      temp_toptable <- limma::topTable(data.fit.eb, 
                                       sort.by = "P",
                                       n = Inf,
                                       confint = TRUE)
      
      # Calculate SE of log2FC
      temp_toptable$SE <- (temp_toptable$CI.R - temp_toptable$CI.L)/3.92
      
      # get top table
      top_table[[i]] <- temp_toptable[,c("AveExpr",
                                         "logFC",
                                         "SE",
                                         "P.Value",
                                         "adj.P.Val")]
      top_table[[i]] <- cbind(rownames(top_table[[i]]), top_table[[i]])
      rownames(top_table[[i]]) <- NULL
      colnames(top_table[[i]]) <- c("GeneID", "meanExpr", "log2FC", "log2FC SE",
                                    "p-value", "adj. p-value")
    }
    message <- "Statistical analysis has been performed. 
    You can download the results as well as view them in interactive plots."
    
    # Add annotations to table if this option is selected
    if (addAnnotation == TRUE){
      
      # Change attribute and filter name
      if (biomart_dataset == "hsapiens_gene_ensembl"){
        biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                    "gene_name",
                                                    "hgnc_symbol")
        biomart_filters1 <- stringr::str_replace(biomart_filters,
                                                 "gene_name",
                                                 "hgnc_symbol")
      } else{
        biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                    "gene_name",
                                                    "external_gene_name")
        biomart_filters1 <- stringr::str_replace(biomart_filters,
                                                 "gene_name",
                                                 "external_gene_name")
      }
      
      
      # Do this for each top table in the list
      for (t in 1:length(top_table)){
        
        # Round numbers in top table
        for (n in 2:6){
          top_table[[t]][,n] <- signif(top_table[[t]][,n],3)
        }
        
        #Get annotations
        annotation_list <- tryCatch({
          ensembl <- biomaRt::useMart("ensembl")
          ensembl <- biomaRt::useDataset(biomart_dataset, mart=ensembl)
          annotations <- biomaRt::getBM(attributes=biomart_attributes1,
                                        filters = biomart_filters1,
                                        values = top_table[[t]]$GeneID,
                                        mart = ensembl)
          
          message <- "Statistical analysis has been performed. 
          Gene annotation was performed with biomaRt. You can download 
              the results as well as view them in interactive plots."
          dataset <- paste0(biomart_dataset, " (", searchDatasets(mart = ensembl, pattern = "hsapiens")$version, ")")
          list(annotations, message, dataset)
        },
        error = function(cond){
          
          # Load annotation package
          pkg <- switch(biomart_dataset,
                        "hsapiens_gene_ensembl" = "org.Hs.eg.db",
                        "btaurus_gene_ensembl" = "org.Bt.eg.db",
                        "celegans_gene_ensembl" = "org.Ce.eg.db",
                        "mmusculus_gene_ensembl" = "org.Mm.eg.db",
                        "rnorvegicus_gene_ensembl" = "org.Rn.eg.db"
          )
          
          if (!requireNamespace(pkg, quietly = TRUE))
            BiocManager::install(pkg, ask = FALSE)
          require(as.character(pkg), character.only = TRUE)
          
          
          # Change attribute names
          biomart_attributes2 <- biomart_attributes
          if (biomart_dataset == "hsapiens_gene_ensembl"){
            biomart_attributes2[biomart_attributes2 == "gene_name"] <- "SYMBOL"
          } else{
            biomart_attributes2[biomart_attributes2 == "gene_name"] <- "GENENAME"
          }
          biomart_attributes2[biomart_attributes2 == "entrezgene_id"] <- "ENTREZID"
          biomart_attributes2[biomart_attributes2 == "ensembl_gene_id"] <- "ENSEMBL"
          
          
          # Change filter names
          biomart_filters2 <- biomart_filters
          if (biomart_dataset == "hsapiens_gene_ensembl"){
            biomart_filters2[biomart_filters2 == "gene_name"] <- "SYMBOL"
          } else{
            biomart_filters2[biomart_filters2  == "gene_name"] <- "GENENAME"
          }
          biomart_filters2[biomart_filters2 == "entrezgene_id"] <- "ENTREZID"
          biomart_filters2[biomart_filters2 == "ensembl_gene_id"] <- "ENSEMBL"
          
          # Get annotations
          annotations <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                               columns = biomart_attributes2, 
                                               keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
          
          # Join with geneIDs
          annotations <- left_join(data.frame(GeneID = top_table[[t]]$GeneID),
                                   annotations,
                                   by = c("GeneID" = biomart_filters2))
          
          # Change column names back
          temp_col <- colnames(annotations)
          temp_col[temp_col == "GeneID"] <- biomart_filters1
          if (biomart_dataset == "hsapiens_gene_ensembl"){
            temp_col[temp_col == "SYMBOL"] <- "hgnc_symbol"
          } else{
            temp_col[temp_col == "SYMBOL"] <- "external_gene_name"
          }
          temp_col[temp_col == "ENTREZID"] <- "entrezgene_id"
          temp_col[temp_col == "ENSEMBL"] <- "ensembl_gene_id"
          colnames(annotations) <- temp_col
          message <- "Statistical analysis has been performed. 
          The Ensembl database was not available.
          So, the gene annotation was performed with the bioconductor annotation package. 
          You can download the results as well as view them in interactive plots."
          dataset <- paste0(pkg, " (", packageVersion(pkg),")")
          list(annotations, message, dataset)
        })
        annotations <- annotation_list[[1]]
        message <- annotation_list[[2]]
        dataset <- annotation_list[[3]]
        
        # Convert entrezgene id to character
        if("entrezgene_id" %in% biomart_attributes){
          annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
        }
        annotations[annotations == ""] <- NA
        annotations[annotations == " "] <- NA
        
        # Combine annotations with top table
        #annotations[,biomart_filters] <- as.character(annotations[,biomart_filters1])
        top_table_ann <- dplyr::left_join(top_table[[t]], annotations,
                                          by = c("GeneID" = biomart_filters1))
        
        colnames(top_table_ann)[colnames(top_table_ann) == "hgnc_symbol"] <- "gene_name"
        colnames(top_table_ann)[colnames(top_table_ann) == "external_gene_name"] <- "gene_name"
        
        # Make sure that there are no duplicate gene ids
        for (a in 1:(ncol(top_table_ann)-6)){
          temp1 <- unique(top_table_ann[,c(1,6+a)])
          temp1 <- temp1[!is.na(temp1[,2]),]
          temp_ann <- temp1 %>%
            dplyr::group_by(GeneID) %>%
            dplyr::summarise_at(colnames(top_table_ann)[6+a], function(x) paste(x,collapse = "; "))
          
          top_table[[t]] <- dplyr::left_join(top_table[[t]], temp_ann,
                                             by = c("GeneID" = "GeneID"))
        }
        
      }
      
    }
    top_table_list <- list(top_table, message, dataset)
    return(top_table_list)
  }, error = function(cond){
    NULL
  })
}

#==============================================================================#
# makelogFCHistogram()
#==============================================================================#

# DESCRIPTION:
# Make log2FC value histogram

# VARIABLES:
# logFC: vector of log2FCs

makelogFCHistogram <- function(logFC, color = "#d3d3d3", bins = 100, static = FALSE){

  if (is.null(color)){
    color <- "#d3d3d3"
  }
  
  if (is.null(bins)){
    bins <- 100
  }
  
  # Make edge color
  color_edge <- colorspace::darken(color, amount = 0.5)
  
  # Collect logFCs into data frame
  plotDF <- data.frame(Value = logFC)
  
  # Make plot
  p <- ggplot2::ggplot(data = plotDF, ggplot2::aes(x = Value)) +
    ggplot2::geom_histogram(bins = bins, colour = color_edge, fill = color) +
    ggplot2::labs(title = "logFC histogram") +
    ggplot2::xlab("logFC") +
    ggplot2::ylab("Count") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  
  # Return plot
  if (!static){
    return(plotly::ggplotly(p))
  } else{
    return(p)
  }
}

#==============================================================================#
# makePHistogram()
#==============================================================================#

# DESCRIPTION:
# Make P value histogram

# VARIABLES:
# P: vector of P values

makePHistogram <- function(P, color = "#d3d3d3", bins = 100, static = FALSE){
  
  if (is.null(color)){
    color <- "#d3d3d3"
  }
  
  if (is.null(bins)){
    bins <- 100
  }
  
  # Make edge color
  color_edge <- colorspace::darken(color, amount = 0.5)
  
  # Collect P values into data frame
  plotDF <- data.frame(Value = P)
  
  # Make plot
  p <- ggplot2::ggplot(data = plotDF, ggplot2::aes(x = Value)) +
    ggplot2::geom_histogram(bins = bins, colour = color_edge, fill = color) +
    ggplot2::labs(title = "P value histogram") +
    ggplot2::xlab("P value") +
    ggplot2::ylab("Count") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  
  # Return plot
  if (!static){
    return(plotly::ggplotly(p))
  } else{
    return(p)
  }
}

#==============================================================================#
# makeVolcano()
#==============================================================================#

# First, make a function for horizontal and vertical lines in plotly:

# Vertical line function
vline <- function(x = 0, color = "grey") {
  list(
    type = "line", 
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, dash = "dot")
  )
}

# Horizontal line function
hline <- function(y = 0, color = "grey") {
  list(
    type = "line", 
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color, dash = "dot")
  )
}

# DESCRIPTION:
# Make volcano plot

# VARIABLES:
# top_table: top table from the "getStatistics" function
# p: raw ("raw") or adjusted ("adj") P value
# p_threshold: P value threshold
# logFC_threshold: log2FC threshold

makeVolcano <- function(top_table, 
                        p = "raw", 
                        p_threshold = 0.05, 
                        logFC_threshold = 1,
                        unchanged_color = "darkgrey",
                        down_color = "blue",
                        up_color = "red"){
  
  if (is.null(unchanged_color)){
    color <- "darkgrey"
  }
  if (is.null(down_color)){
    color <- "blue"
  }
  if (is.null(up_color)){
    color <- "red"
  }
  
  plotDF <- data.frame(log2FC = top_table[,"log2FC"],
                       Pvalue = top_table[,"p-value"],
                       logP = -log10(top_table[,"p-value"]),
                       adjPvalue = top_table[,"adj. p-value"],
                       logadjP = -log10(top_table[,"adj. p-value"]),
                       GeneID = top_table[,"GeneID"])
  
  # Format pop-up text
  if (ncol(top_table) == 7){
    genes7 <- apply(top_table, 1, function(x) ifelse(nchar(x[7])>15,paste0(substr(x[7],1,15),", ..."), x[7]))
    plotDF$GeneID <- paste0("Gene ID: ", plotDF$GeneID, "\n", 
                            colnames(top_table)[7], ": ", genes7)
  }
  
  if (ncol(top_table) == 8){
    genes7 <- apply(top_table, 1, function(x) ifelse(nchar(x[7])>15,paste0(substr(x[7],1,15),", ..."), x[7]))
    genes8 <- apply(top_table, 1, function(x) ifelse(nchar(x[8])>15,paste0(substr(x[8],1,15),", ..."), x[8]))
    plotDF$GeneID <- paste0("Gene ID: ", plotDF$GeneID, "\n", 
                            colnames(top_table)[7], ": ", genes7, "\n", 
                            colnames(top_table)[8], ": ", genes8)
  }
  
  
  if (p == "raw"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$Pvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$Pvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$Pvalue <= p_threshold] <- "upregulated"
    
    
    volcano <- plotly::plot_ly() %>%
      # unchanged
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "unchanged",],
        x = ~log2FC,
        y = ~logP,
        type = "scatter",
        mode = "markers",
        marker = list(color = unchanged_color),
        hoverinfo = "none",
        name = "Unchanged"
      ) %>%
      # downregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "downregulated",],
        x = ~log2FC,
        y = ~logP,
        type = "scatter",
        mode = "markers",
        marker = list(color = down_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Downregulated"
      ) %>%
      # upregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "upregulated",],
        x = ~log2FC,
        y = ~logP,
        type = "scatter",
        mode = "markers",
        marker = list(color = up_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Upregulated"
      )%>%
      layout(xaxis = list(title = 'log<sub>2</sub>FC'), 
             yaxis = list(title = '-log<sub>10</sub> P value'),
             shapes = list(vline(logFC_threshold), 
                           vline(-1*logFC_threshold),
                           hline(-log10(p_threshold))),
             legend = list(traceorder = "alphabetical"))
  }
  
  
  
  
  
  if (p == "adj"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"
    
    
    volcano <- plotly::plot_ly() %>%
      # unchanged
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "unchanged",],
        x = ~log2FC,
        y = ~logadjP,
        type = "scatter",
        mode = "markers",
        marker = list(color = unchanged_color),
        hoverinfo = "none",
        name = "Unchanged"
      ) %>%
      # downregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "downregulated",],
        x = ~log2FC,
        y = ~logadjP,
        type = "scatter",
        mode = "markers",
        marker = list(color = down_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Downregulated"
      ) %>%
      # upregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "upregulated",],
        x = ~log2FC,
        y = ~logadjP,
        type = "scatter",
        mode = "markers",
        marker = list(color = up_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Upregulated"
      )%>%
      layout(xaxis = list(title = 'log<sub>2</sub>FC'), 
             yaxis = list(title = '-log<sub>10</sub> adj. P value'),
             shapes = list(vline(logFC_threshold), 
                           vline(-1*logFC_threshold),
                           hline(-log10(p_threshold))),
             legend = list(traceorder = "alphabetical"))
  }
  
  #p <- plotly::ggplotly(volcano, tooltip = "text")
  
  return(volcano)
  
  
}


#==============================================================================#
# makeVolcano_static()
#==============================================================================#

# DESCRIPTION:
# Make volcano plot (static)

# VARIABLES:
# top_table: top table from the "getStatistics" function
# p: raw ("raw") or adjusted ("adj") P value
# p_threshold: P value threshold
# logFC_threshold: log2FC threshold

makeVolcano_static <- function(top_table, 
                        p = "raw", 
                        p_threshold = 0.05, 
                        logFC_threshold = 1,
                        unchanged_color = "darkgrey",
                        down_color = "blue",
                        up_color = "red"){
  
  if (is.null(unchanged_color)){
    color <- "darkgrey"
  }
  if (is.null(down_color)){
    color <- "blue"
  }
  if (is.null(up_color)){
    color <- "red"
  }
  
  plotDF <- data.frame(log2FC = top_table[,"log2FC"],
                       Pvalue = top_table[,"p-value"],
                       logP = -log10(top_table[,"p-value"]),
                       adjPvalue = top_table[,"adj. p-value"],
                       GeneID = top_table[,"GeneID"])
  
  
  if (p == "raw"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$Pvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$Pvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$Pvalue <= p_threshold] <- "upregulated"
    
    
    volcano <- ggplot2::ggplot() +
      ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                           ggplot2::aes(x = log2FC, y = logP, color = "Unchanged")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                          ggplot2::aes(x = log2FC, y = logP, color = "Downregulated")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                          ggplot2::aes(x = log2FC, y = logP, color = "Upregulated")) +
      ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_vline(xintercept = -1*logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_vline(xintercept = logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                           c("Unchanged", "Downregulated", "Upregulated"))) +
      ggplot2::xlab(expression(log[2]~"FC")) +
      ggplot2::ylab(expression(-log[10]~"P value")) +
      ggplot2::labs(color = NULL) +
      ggplot2::theme_minimal()
    
   
  }
  
  if (p == "adj"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"
    
    
    volcano <- ggplot2::ggplot() +
      ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                           ggplot2::aes(x = log2FC, y = logP, color = "Unchanged")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                          ggplot2::aes(x = log2FC, y = logP, color = "Downregulated")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                          ggplot2::aes(x = log2FC, y = logP, color = "Upregulated")) +
      ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_vline(xintercept = -1*logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_vline(xintercept = logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                    c("Unchanged", "Downregulated", "Upregulated"))) +
      ggplot2::xlab(expression(log[2]~"FC")) +
      ggplot2::ylab(expression(-log[10]~"P value")) +
      ggplot2::labs(color = NULL) +
      ggplot2::theme_minimal()
    
  }
  
  return(volcano)
  
  
}

#==============================================================================#
# makeMAplot()
#==============================================================================#

# DESCRIPTION:
# Make MA plot

# VARIABLES:
# top_table: top table from the "getStatistics" function
# p: raw ("raw") or adjusted ("adj") P value
# p_threshold: P value threshold
# logFC_threshold: log2FC threshold

makeMAplot <- function(top_table, 
                       p = "raw", 
                       p_threshold = 0.05, 
                       logFC_threshold = 1,
                       unchanged_color = "darkgrey",
                       down_color = "blue",
                       up_color = "red",
                       RNAseq = TRUE){
  
  if (is.null(unchanged_color)){
    color <- "darkgrey"
  }
  if (is.null(down_color)){
    color <- "blue"
  }
  if (is.null(up_color)){
    color <- "red"
  }
  
  plotDF <- data.frame(log2FC = top_table[,"log2FC"],
                       Pvalue = top_table[,"p-value"],
                       logP = -log10(top_table[,"p-value"]),
                       adjPvalue = top_table[,"adj. p-value"],
                       meanExpr = top_table[,"meanExpr"],
                       GeneID = top_table[,"GeneID"])
  
  if (isTRUE(RNAseq)){
    plotDF$meanExpr <- log2(plotDF$meanExpr +1)
  }
  
  # Format pop-up text
  if (ncol(top_table) == 7){
    genes7 <- apply(top_table, 1, function(x) ifelse(nchar(x[7])>15,paste0(substr(x[7],1,15),", ..."), x[7]))
    plotDF$GeneID <- paste0("Gene ID: ", plotDF$GeneID, "\n", 
                            colnames(top_table)[7], ": ", genes7)
  }
  
  if (ncol(top_table) == 8){
    genes7 <- apply(top_table, 1, function(x) ifelse(nchar(x[7])>15,paste0(substr(x[7],1,15),", ..."), x[7]))
    genes8 <- apply(top_table, 1, function(x) ifelse(nchar(x[8])>15,paste0(substr(x[8],1,15),", ..."), x[8]))
    plotDF$GeneID <- paste0("Gene ID: ", plotDF$GeneID, "\n", 
                            colnames(top_table)[7], ": ", genes7, "\n", 
                            colnames(top_table)[8], ": ", genes8)
  }
  
  if (p == "raw"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$Pvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$Pvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$Pvalue <= p_threshold] <- "upregulated"
    
    
    MA <- plotly::plot_ly() %>%
      # unchanged
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "unchanged",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = unchanged_color),
        hoverinfo = "none",
        name = "Unchanged"
      ) %>%
      # downregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "downregulated",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = down_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Downregulated"
      ) %>%
      # upregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "upregulated",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = up_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Upregulated"
      )%>%
      layout(yaxis = list(title = 'log<sub>2</sub>FC'), 
             xaxis = list(title = 'Normalized log<sub>2</sub> expression'),
             shapes = list(hline(logFC_threshold), 
                           hline(-1*logFC_threshold)))
  }
  
  
  
  
  
  if (p == "adj"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"
    
    
    MA <- plotly::plot_ly() %>%
      # unchanged
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "unchanged",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = unchanged_color),
        hoverinfo = "none",
        name = "Unchanged"
      ) %>%
      # downregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "downregulated",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = down_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Downregulated"
      ) %>%
      # upregulated
      plotly::add_trace(
        data = plotDF[plotDF$Colour == "upregulated",],
        x = ~meanExpr,
        y = ~log2FC,
        type = "scatter",
        mode = "markers",
        marker = list(color = up_color),
        hoverinfo = "text",
        text = ~GeneID,
        name = "Upregulated"
      )%>%
      layout(yaxis = list(title = 'log<sub>2</sub>FC'), 
             xaxis = list(title = 'Normalized log<sub>2</sub> expression'),
             shapes = list(hline(logFC_threshold), 
                           hline(-1*logFC_threshold)))
  }
  
  return(MA)
  
  
}


#==============================================================================#
# makeMAplot_static()
#==============================================================================#

# DESCRIPTION:
# Make MA plot

# VARIABLES:
# top_table: top table from the "getStatistics" function
# p: raw ("raw") or adjusted ("adj") P value
# p_threshold: P value threshold
# logFC_threshold: log2FC threshold

makeMAplot_static <- function(top_table, 
                       p = "raw", 
                       p_threshold = 0.05, 
                       logFC_threshold = 1,
                       unchanged_color = "darkgrey",
                       down_color = "blue",
                       up_color = "red",
                       RNAseq = TRUE){
  
  if (is.null(unchanged_color)){
    color <- "darkgrey"
  }
  if (is.null(down_color)){
    color <- "blue"
  }
  if (is.null(up_color)){
    color <- "red"
  }
  
  plotDF <- data.frame(log2FC = top_table[,"log2FC"],
                       Pvalue = top_table[,"p-value"],
                       logP = -log10(top_table[,"p-value"]),
                       adjPvalue = top_table[,"adj. p-value"],
                       meanExpr = top_table[,"meanExpr"],
                       GeneID = top_table[,"GeneID"])
  
  if (isTRUE(RNAseq)){
    plotDF$meanExpr <- log2(plotDF$meanExpr +1)
  }
  
  if (p == "raw"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$Pvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$Pvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$Pvalue <= p_threshold] <- "upregulated"
    
    
    MA <- ggplot2::ggplot() +
      ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                           ggplot2::aes(x = meanExpr, y = log2FC, color = "Unchanged")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                          ggplot2::aes(x = meanExpr, y = log2FC, color = "Downregulated")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                          ggplot2::aes(x = meanExpr, y = log2FC, color = "Upregulated")) +
      ggplot2::geom_hline(yintercept = -1*logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_hline(yintercept = logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                    c("Unchanged", "Downregulated", "Upregulated"))) +
      ggplot2::xlab(expression("Normalized"~log[2]~"expression")) +
      ggplot2::ylab(expression(log[2]~"FC")) +
      ggplot2::labs(color = NULL) +
      ggplot2::theme_minimal()
  }
  
  
  
  
  
  if (p == "adj"){
    plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
    plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
    plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"
    
    
    
    MA <- ggplot2::ggplot() +
      ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                           ggplot2::aes(x = meanExpr, y = log2FC, color = "Unchanged")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                          ggplot2::aes(x = meanExpr, y = log2FC, color = "Downregulated")) +
      ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                          ggplot2::aes(x = meanExpr, y = log2FC, color = "Upregulated")) +
      ggplot2::geom_hline(yintercept = -1*logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_hline(yintercept = logFC_threshold, linetype = "dashed", color = "lightgrey") +
      ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                    c("Unchanged", "Downregulated", "Upregulated"))) +
      ggplot2::xlab(expression("Normalized"~log[2]~"expression")) +
      ggplot2::ylab(expression(log[2]~"FC")) +
      ggplot2::labs(color = NULL) +
      ggplot2::theme_minimal()
  }
  
  return(MA)
  
  
}
#==============================================================================#
# selID()
#==============================================================================#

# DESCRIPTION:
# Guess which gene ID can be used for the gene overrepresentation analysis (ORA)

# VARIABLES:
# ProbeAnnotation: selected ProbeAnnotation

selID <- function(ProbeAnnotation){
  
  selID <- "ENTREZID"
  if (!is.null(ProbeAnnotation)){
    if (ProbeAnnotation == "ENSG"){
      selID <- "ENSEMBL"
    }
  }
  return(selID)
}

#==============================================================================#
# ORA()
#==============================================================================#

# DESCRIPTION:
# Perform gene overrepresentation analysis (ORA)

# VARIABLES:
# top_table: top table from the "getStatistics" function
# geneset: which geneset to perform ORA on (GO-BP, GO-MF, GO-CC, or WikiPathways)
# geneID_col: Column that contains the gene IDs
# geneID_type: Which gene ID does the selected column contain 
#              (ENTREZID, ENSEMBL, SYMBOL)
# organism: organism (e.g., Homo sapiens, Mus musculus)
# rawadj: raw ("raw") or adjusted ("adj) P value
# updown: Use upregulated genes ("Upregulated genes only"), 
#         downregulated genes ("Upregulated genes only"), or both ("Both) 
#         for the ORA
# p_thres: P value treshold
# logFC_thres: logFC threshold

ORA <- function(top_table,
                geneset = "GO-BP",
                geneID_col = colnames(top_table)[1],
                geneID_type = "ENTREZID",
                organism = "Homo sapiens",
                updown = "Both",
                topN = FALSE,
                N = NULL,
                rawadj = "raw",
                p_thres = 0.05,
                logFC_thres = 0){
  
  tryCatch({
    
    
    #Required Bioconductor annotation packages:
    pkg <- switch(organism,
                  "Homo sapiens" = "org.Hs.eg.db",
                  "Bos taurus" = "org.Bt.eg.db",
                  "Caenorhabditis elegans" = "org.Ce.eg.db",
                  "Mus musculus" = "org.Mm.eg.db",
                  "Rattus norvegicus" = "org.Rn.eg.db"
    )
    
    
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE)
    require(as.character(pkg), character.only = TRUE)
    
    #..........................................................................#
    # Filter top table for background genes (in genesets)
    #..........................................................................#
    
    # get all background genes
    if (geneset == "WikiPathways"){
      load("Objects/gmt_WP_all.RData")
      gmt_all <- gmt_all[[organism]]
      bg_genes <- as.character(unique(gmt_all[,geneID_type]))
      
    } 
    if (geneset == "GO-BP"){
      bg_genes <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                        columns = c(geneID_type, "GO"), 
                                        keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
      
      bg_genes <- unique(bg_genes[bg_genes$ONTOLOGY=="BP",geneID_type])
      bg_genes <- bg_genes[!is.na(bg_genes)]
    }
    if (geneset == "GO-MF"){
      bg_genes <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                        columns = c(geneID_type, "GO"), 
                                        keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
      
      bg_genes <- unique(bg_genes[bg_genes$ONTOLOGY=="MF",geneID_type])
      bg_genes <- bg_genes[!is.na(bg_genes)]
    }
    if (geneset == "GO-CC"){
      bg_genes <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                        columns = c(geneID_type, "GO"), 
                                        keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
      
      bg_genes <- unique(bg_genes[bg_genes$ONTOLOGY=="CC",geneID_type])
      bg_genes <- bg_genes[!is.na(bg_genes)]
    }
    if (geneset == "KEGG"){
      gmt_KEGG <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                        columns = c(geneID_type, "PATH"), 
                                        keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
      gmt_KEGG <- gmt_KEGG[!is.na(gmt_KEGG$PATH),]
      
      bg_genes <- unique(gmt_KEGG[, geneID_type])
      
    }
    
    
    # If the gene IDs are not saved in the first column, some genes might be 
    # separated by a comma in the same row.
    if (which(colnames(top_table) == geneID_col) != 1){
      
      # Split all genes by a comma
      genes_all <- str_split(top_table[,geneID_col], "; ")
      
      # Combine all genes into a single vector
      single_name <- NULL
      combined_name <- NULL
      for (i in 1:length(genes_all)){
        single_name <- c(single_name, genes_all[[i]])
        combined_name <- c(combined_name, rep(top_table[i,geneID_col], length(genes_all[[i]])))
      }
      gene_names <- data.frame(single_name,
                               combined_name)
      gene_names <- gene_names[!is.na(single_name),]
      
      # Make separate row names for single gene name
      colnames(top_table)[colnames(top_table) == geneID_col] <- "ID"
      top_table_name <- left_join(top_table, gene_names, by = c("ID" = "combined_name"))
      colnames(top_table_name)[colnames(top_table_name) == "ID"] <- geneID_col
      
      # Filter top table
      top_table_fil <- unique(top_table_name[top_table_name[,"single_name"] %in% bg_genes,-which(colnames(top_table_name) == "single_name")])
      
    } else { # If all the gene IDs are in the first column:
      
      # Filter top table
      top_table_fil <- top_table[top_table[,geneID_col] %in% bg_genes,]
    }
    
    
    #..........................................................................#
    # Perform ORA on filtered top table
    #..........................................................................#
    top_table <- top_table_fil
    
    # Universe: background list of genes
    universe <- unlist(stringr::str_split(top_table[,geneID_col], "; "))
    
    # Filter top_table for up- or downregulated genes only
    if (updown == "Upregulated genes only"){
      top_table <- top_table[top_table[,"log2FC"] > 0,]
    }  
    if (updown == "Downregulated genes only"){
      top_table <- top_table[top_table[,"log2FC"] < 0,]
    }
    
    # Select differentially expressed gens (DEGs, geneList)
    
    if (isTRUE(topN)){
      top_table_sort <- arrange(top_table, `p-value`)
      geneList <- head(top_table,N)[,geneID_col]
      geneList <- unlist(stringr::str_split(geneList, "; "))
      
    } else{
      if (rawadj == "raw"){
        geneList <- top_table[((abs(top_table[,"log2FC"])) > logFC_thres) & 
                                (top_table[,"p-value"] < p_thres),
                              geneID_col]
        geneList <- unlist(stringr::str_split(geneList, "; "))
      }
      if (rawadj == "adj"){
        geneList <- top_table[((abs(top_table[,"log2FC"])) > logFC_thres) & 
                                (top_table[,"adj. p-value"] < p_thres),
                              geneID_col]
        geneList <- unlist(stringr::str_split(geneList, "; "))
      }
    }
    
    
    # perform ORA
    ORA_data <- NULL
    if (!is.null(geneList)){
      
      # GO-BP
      if (geneset == "GO-BP"){
        ORA_data <- clusterProfiler::enrichGO(
          gene = geneList,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "BP",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          universe = universe,
          qvalueCutoff = Inf,
          minGSSize = 10,
          maxGSSize = 500,
          readable = FALSE,
          pool = FALSE
        )
      }
      
      # GO-MF
      if (geneset == "GO-MF"){
        ORA_data <- clusterProfiler::enrichGO(
          gene = geneList,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "MF",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          universe = universe,
          qvalueCutoff = Inf,
          minGSSize = 10,
          maxGSSize = 500,
          readable = FALSE,
          pool = FALSE
        )
      }
      
      # GO-CC
      if (geneset == "GO-CC"){
        ORA_data <- clusterProfiler::enrichGO(
          gene = geneList,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "CC",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          universe = universe,
          qvalueCutoff = Inf,
          minGSSize = 10,
          maxGSSize = 500,
          readable = FALSE,
          pool = FALSE
        )
      }
      
      # WikiPathways
      if (geneset == "WikiPathways"){
        
        # load GMT file
        load("Objects/gmt_WP_all.RData")
        gmt_all <- gmt_all[[organism]]
        print(head(gmt_all))
        
        # Prepare GMT for analysis
        gmt <- unique(gmt_all[,c("name", "version", "wpid",
                                 "species", geneID_type)])
        print(head(gmt))
        path2gene <- gmt[,c("wpid", geneID_type)]
        path2name <- gmt[,c("wpid", "name")]
        
        # Perform ORA
        ORA_data <- clusterProfiler::enricher(
          gene = geneList,
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          universe = universe,
          minGSSize = 10,
          maxGSSize = 500,
          qvalueCutoff = Inf,
          gson = NULL,
          TERM2GENE = path2gene,
          TERM2NAME = path2name
        )
        
      }
      
      # KEGG
      if (geneset == "KEGG"){
        
        # Get correct organism name
        KEGG_org <- switch(organism,
                      "Homo sapiens" = "hsa",
                      "Bos taurus" = "bta",
                      "Caenorhabditis elegans" = "cel",
                      "Mus musculus" = "mmu",
                      "Rattus norvegicus" = "rno"
        )
        
        # Prepare GMT for analysis
        gmt <- unique(gmt_KEGG[,c("PATH", geneID_type)])
        path2gene <- gmt[,c("PATH", geneID_type)]
        colnames(path2gene) <- c("path", "gene")
        path2gene$path <- paste0(KEGG_org,path2gene$path)
        path2name <- read.table(url(paste0("https://rest.kegg.jp/list/pathway/", KEGG_org)), sep = "\t")
        colnames(path2name) <- c("path", "name")
        
        # Perform ORA
        ORA_data <- clusterProfiler::enricher(
          gene = geneList,
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          universe = universe,
          minGSSize = 10,
          maxGSSize = 500,
          qvalueCutoff = Inf,
          gson = NULL,
          TERM2GENE = path2gene,
          TERM2NAME = path2name
        )
      }
    }
    
    output <- ORA_data@result[,c("ID", "Description", "GeneRatio",
                                 "BgRatio", "pvalue", "p.adjust",
                                 "qvalue", "geneID", "Count")]
    
    output$pvalue <- signif(output$pvalue,3)
    output$p.adjust <- signif(output$p.adjust,3)
    output$qvalue <- signif(output$qvalue,3)
    
    rownames(output) <- NULL
    colnames(output) <- c("ID", "Description", "GeneRatio", "BgRatio",
                          "p-value", "adj. p-value", "q-value", "GeneIDs", "Count")
    output <- output[,c("ID", "Description", "GeneRatio", "BgRatio",
                        "p-value", "adj. p-value")]
    
    # make term names shorter
    output$Description[nchar(output$Description)>50] <- paste0(substring(output$Description[nchar(output$Description)>50],1,47),"...")
    
    ORA_data@result <- arrange(output, by = `p-value`)
    
    return(ORA_data)
  }, error = function(cond){
    NULL
  })
}

#==============================================================================#
# performGSEA()
#==============================================================================#

# DESCRIPTION:
# Perform geneset enrichment analysis (GSEA)

# VARIABLES:
# top_table: top table from the "getStatistics" function
# geneset: which geneset to perform ORA on (GO-BP, GO-MF, GO-CC, or WikiPathways)
# geneID_col: Column that contains the gene IDs
# geneID_type: Which gene ID does the selected column contain 
#              (ENTREZID, ENSEMBL, SYMBOL)
# organism: organism (e.g., Homo sapiens, Mus musculus)
# rawadj: raw ("raw") or adjusted ("adj) P value
# updown: Use upregulated genes ("Upregulated genes only"), 
#         downregulated genes ("Upregulated genes only"), or both ("Both) 
#         for the ORA
# p_thres: P value treshold
# logFC_thres: logFC threshold

performGSEA <- function(top_table,
                geneset = "GO-BP",
                geneID_col = colnames(top_table)[1],
                geneID_type = "ENTREZID",
                organism = "Homo sapiens",
                rankingVar = "logFC"){
  
  tryCatch({
    
    
    #Required Bioconductor annotation packages:
    pkg <- switch(organism,
                  "Homo sapiens" = "org.Hs.eg.db",
                  "Bos taurus" = "org.Bt.eg.db",
                  "Caenorhabditis elegans" = "org.Ce.eg.db",
                  "Mus musculus" = "org.Mm.eg.db",
                  "Rattus norvegicus" = "org.Rn.eg.db"
    )
    
    
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, ask = FALSE)
    require(as.character(pkg), character.only = TRUE)
    
    
    # If the gene IDs are not saved in the first column, some genes might be 
    # separated by a comma in the same row.
    if (which(colnames(top_table) == geneID_col) != 1){
      
      # Split all genes by a comma
      genes_all <- str_split(top_table[,geneID_col], "; ")
      
      # Combine all genes into a single vector
      single_name <- NULL
      combined_name <- NULL
      for (i in 1:length(genes_all)){
        single_name <- c(single_name, genes_all[[i]])
        combined_name <- c(combined_name, rep(top_table[i,geneID_col], length(genes_all[[i]])))
      }
      gene_names <- data.frame(single_name,
                               combined_name)
      gene_names <- gene_names[!is.na(single_name),]
      
      # Make separate row names for single gene name
      colnames(top_table)[colnames(top_table) == geneID_col] <- "ID"
      top_table_name <- left_join(top_table, gene_names, by = c("ID" = "combined_name"))
      colnames(top_table_name)[colnames(top_table_name) == "ID"] <- geneID_col
      top_table <- top_table_name
      top_table <- top_table[!is.na(top_table[,geneID_col]),]
    }
    
    
    #..........................................................................#
    # Perform GSEA on filtered top table
    #..........................................................................#
    
    # Make GSEA input vector
    if (rankingVar == "logFC"){
      temp <- arrange(top_table, by = `p-value`)
      temp <- temp[!duplicated(temp[,geneID_col]),]
      gsea_input <- sort(setNames(temp[,"log2FC"],
                                  temp[,geneID_col]), 
                         decreasing = TRUE)
    }
    if (rankingVar == "pvalue"){
      temp <- arrange(top_table, by = `p-value`)
      temp <- temp[!duplicated(temp[,geneID_col]),]
      gsea_input <- sort(setNames(-log10(temp[,"p-value"]),
                                  temp[,geneID_col]), 
                         decreasing = TRUE)
    }
    if (rankingVar == "signed_pvalue"){
      temp <- arrange(top_table, by = `p-value`)
      temp <- temp[!duplicated(temp[,geneID_col]),]
      gsea_input <- sort(setNames(sign(temp[,"log2FC"])*-log10(temp[,"p-value"]),
                                  temp[,geneID_col]), 
                         decreasing = TRUE)
    }
    
    
    # perform ORA
    GSEA_data <- NULL
    if (!is.null(gsea_input)){
      
      # GO-BP
      if (geneset == "GO-BP"){
        set.seed(123)
        GSEA_data <- clusterProfiler::gseGO(
          geneList = gsea_input,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "BP",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          minGSSize = 10,
          maxGSSize = 500
        )
      }
      
      # GO-MF
      if (geneset == "GO-MF"){
        set.seed(123)
        GSEA_data <- clusterProfiler::gseGO(
          geneList = gsea_input,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "MF",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          minGSSize = 10,
          maxGSSize = 500
        )
      }
      
      # GO-CC
      if (geneset == "GO-CC"){
        set.seed(123)
        GSEA_data <- clusterProfiler::gseGO(
          geneList = gsea_input,
          OrgDb = pkg,
          keyType = geneID_type,
          ont = "CC",
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          minGSSize = 10,
          maxGSSize = 500
        )
      }
      
      # WikiPathways
      if (geneset == "WikiPathways"){
        
        # load GMT file
        load("Objects/gmt_WP_all.RData")
        gmt_all <- gmt_all[[organism]]
        
        # Prepare GMT for analysis
        gmt <- unique(gmt_all[,c("name", "version", "wpid",
                                 "species", geneID_type)])
        path2gene <- gmt[,c("wpid", geneID_type)]
        path2name <- gmt[,c("wpid", "name")]
        
        # Perform GSEA
        set.seed(123)
        GSEA_data <- GSEA(
          geneList = gsea_input,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          TERM2GENE = path2gene,
          TERM2NAME = path2name
        )
        
        
      }
      
      # KEGG
      if (geneset == "KEGG"){
        
        # Get correct organism name
        KEGG_org <- switch(organism,
                           "Homo sapiens" = "hsa",
                           "Bos taurus" = "bta",
                           "Caenorhabditis elegans" = "cel",
                           "Mus musculus" = "mmu",
                           "Rattus norvegicus" = "rno"
        )
        
        # Prepare GMT for analysis
        gmt <- unique(gmt_KEGG[,c("PATH", geneID_type)])
        path2gene <- gmt[,c("PATH", geneID_type)]
        colnames(path2gene) <- c("path", "gene")
        path2gene$path <- paste0(KEGG_org,path2gene$path)
        path2name <- read.table(url(paste0("https://rest.kegg.jp/list/pathway/", KEGG_org)), sep = "\t")
        colnames(path2name) <- c("path", "name")
        
        # Perform GSEA
        set.seed(123)
        GSEA_data <- GSEA(
          geneList = gsea_input,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = Inf,
          pAdjustMethod = "fdr",
          TERM2GENE = path2gene,
          TERM2NAME = path2name
        )
      }
    }
    
    # Prepare output
    output <- GSEA_data@result[,c("ID", "Description", "NES", "pvalue", "p.adjust")]
    
    output$pvalue <- signif(output$pvalue,3)
    output$p.adjust <- signif(output$p.adjust,3)
    output$NES <- signif(output$NES,3)
    
    rownames(output) <- NULL
    colnames(output) <- c("ID", "Description", "NES",
                          "p-value", "adj. p-value")
    output <- output[,c("ID", "Description", "NES",
                        "p-value", "adj. p-value")]
    
    # make term names shorter
    output$Description[nchar(output$Description)>50] <- paste0(substring(output$Description[nchar(output$Description)>50],1,47),"...")
  
    # Order results
    GSEA_data@result <- arrange(output, by = `p-value`)
    
    # Return output
    return(GSEA_data)
  }, error = function(cond){
    NULL
  })
}

#==============================================================================#
# make_ORAgene_table()
#==============================================================================#

# DESCRIPTION:
# Make top table for the genes in the selected gene set (GO term/WikiPathways)

# VARIABLES:
# ORA_data: ORA data object (clusterProfiler object)
# top_table: top table from the "getStatistics" function
# geneID_col: column name of the column that contains the gene IDs that were
#             used in the gene overrepresentation analysis
# sel_row_ORA: selected row index of the ORA_data (which gene set?)

make_ORAgene_table <- function(ORA_data,
                               top_table,
                               geneID_col = colnames(top_table)[1],
                               sel_row_ORA){
  
  # If the gene IDs are not saved in the first column, some genes might be 
  # separated by a comma in the same row.
  if (which(colnames(top_table) == geneID_col) != 1){
    
    # Split all genes by a comma
    genes_all <- str_split(top_table[,geneID_col], "; ")
    
    # Combine all genes into a single vector
    single_name <- NULL
    combined_name <- NULL
    for (i in 1:length(genes_all)){
      single_name <- c(single_name, genes_all[[i]])
      combined_name <- c(combined_name, rep(top_table[i,geneID_col], length(genes_all[[i]])))
    }
    gene_names <- data.frame(single_name,
                             combined_name)
    gene_names <- gene_names[!is.na(single_name),]
    
    # Make separate row names for single gene name
    colnames(top_table)[colnames(top_table) == geneID_col] <- "ID"
    top_table_name <- left_join(top_table, gene_names, by = c("ID" = "combined_name"))
    colnames(top_table_name)[colnames(top_table_name) == "ID"] <- geneID_col
    
    # Get genes in GO/WP term
    termID <- ORA_data@result$ID[sel_row_ORA]
    termGenes <- ORA_data@geneSets[[termID]]
    
    # Filter top table
    top_table_fil <- unique(top_table_name[top_table_name[,"single_name"] %in% termGenes,-which(colnames(top_table_name) == "single_name")])
    
  } else { # If all the gene IDs are in the first column:
    
    # Get genes in GO/WP term
    termID <- ORA_data@result$ID[sel_row_ORA]
    termGenes <- ORA_data@geneSets[[termID]]
    
    # Filter top table
    top_table_fil <- top_table[top_table[,geneID_col] %in% termGenes,]
  }
  
  # Remove rownames
  rownames(top_table_fil) <- NULL
  
  # Return results
  return(top_table_fil)
}


#==============================================================================#
# makeORAplot()
#==============================================================================#

# DESCRIPTION:
# Make bar graph of the results from the gene overrepresentation anlysis (ORA)

# VARIABLES:
# ORA_data: ORA data object (clusterProfiler object)
# nSets: number of gene sets (i.e., GO terms/WikiPathways) to include in the plot

makeORAplot <- function(ORA_data, nSets, color = "Viridis", static = FALSE){
  
  # round the number of sets to an integer value
  nSets <- round(nSets)
  
  # Set the minimum number to 5
  if(nSets < 5){
    nSets <- 5
  }
  
  # Set the maximum number to 20
  if(nSets > 20){
    nSets <- 20
  }
  
  # Select color gradient
  gradient = ggplot2::scale_fill_viridis_c()
  
  if (color == "Yellow-red"){
    gradient = ggplot2::scale_fill_gradient(low = "yellow", high = "red")
  }
  if (color == "Blues"){
    gradient = ggplot2::scale_fill_gradient(low = "#DEEBF7", high = "#08519C")
  }
  if (color == "Reds"){
    gradient = ggplot2::scale_fill_gradient(low = "#FEE0D2", high = "#CB181D")
  }
  
  # Retrieve data from ORA object and make a data frame
  plotDF <- ORA_data@result
  plotDF$Name <- paste0(firstup(plotDF$Description), " (", plotDF$ID, ")")
  plotDF$Name <- factor(plotDF$Name, levels = rev(plotDF$Name))
  
  # Make the plot
  p <- ggplot2::ggplot(plotDF[1:nSets,], 
                       ggplot2::aes(text = paste0("ID: ",ID, "\n",
                                                  "Name: ", Description, "\n",
                                                  "Gene ratio: ", GeneRatio, "\n",
                                                  "Background ratio: ", BgRatio, "\n",
                                                  "p-value: ", `p-value`, "\n",
                                                  "FDR-adj. p-value: ", `adj. p-value`, "\n"))) +
    ggplot2::geom_bar(ggplot2::aes(x = -log10(`p-value`), 
                                   y = Name, 
                                   fill = -log10(`p-value`)),
                      position = ggplot2::position_dodge(), 
                      stat = "identity", 
                      color = "black",
                      size = 0.5) +
    ggplot2::xlab(expression(-log[10]~"p-value")) +
    ggplot2::ylab(NULL) +
    ggplot2::theme_minimal() +
    gradient +
    ggplot2::theme(legend.position = "none")
  
  
  if (static){
    return(p)
  }else{
    p1 <- plotly::ggplotly(p, tooltip = "text") %>% 
      #plotly::layout(height = 500, width = 1000) %>%
      layout(xaxis = list(title = '-log<sub>10</sub> p-value'))
    return(p1)
  }
  
  return(p1)
}


#==============================================================================#
# makeGSEAplot()
#==============================================================================#

# DESCRIPTION:
# Make bar graph of the results from the geneset enrichment analysis (GSEA)

# VARIABLES:
# GSEA_data: GSEA data object (clusterProfiler object)
# nSets: number of gene sets (i.e., GO terms/WikiPathways) to include in the plot

makeGSEAplot <- function(GSEA_data, nSets, color, static = FALSE){
  
  if (is.null(color)){
    color <- c("#000072", "white", "red")
  }
  if (is.null(nSets)){
    nSets <- 10
  }
  
  # round the number of sets to an integer value
  nSets <- round(nSets)
  
  # Set the minimum number to 5
  if(nSets < 5){
    nSets <- 5
  }
  
  # Set the maximum number to 20
  if(nSets > 20){
    nSets <- 20
  }
  
  # Retrieve data from ORA object and make a data frame
  plotDF <- GSEA_data@result
  plotDF$Name <- paste0(firstup(plotDF$Description), " (", plotDF$ID, ")")
  plotDF$Name <- factor(plotDF$Name, levels = rev(plotDF$Name))
  
  # Make the plot
  p <- ggplot2::ggplot(plotDF[1:nSets,], 
                       ggplot2::aes(text = paste0("ID: ",ID, "\n",
                                                  "Name: ", Description, "\n",
                                                  "NES: ", NES, "\n",
                                                  "p-value: ", `p-value`, "\n",
                                                  "FDR-adj. p-value: ", `adj. p-value`, "\n"))) +
    ggplot2::geom_bar(ggplot2::aes(x = -log10(`p-value`), 
                                   y = Name, 
                                   fill = NES),
                      position = ggplot2::position_dodge(), 
                      stat = "identity", 
                      color = "black",
                      size = 0.5) +
    ggplot2::xlab(expression(-log[10]~"p-value")) +
    ggplot2::ylab(NULL) +
    ggplot2::theme_minimal() +
    scale_fill_gradient2(low = color[1], mid = color[2], high = color[3], 
                         midpoint = 0) +
    ggplot2::theme(legend.position = "right")
  
  
  if (static){
    return(p)
  }else{
    p1 <- plotly::ggplotly(p, tooltip = "text") %>% 
      #plotly::layout(height = 500, width = 1000) %>%
      layout(xaxis = list(title = '-log<sub>10</sub> p-value'))
    return(p1)
  }

}

#==============================================================================#
# JI()
#==============================================================================#

# DESCRIPTION:
# Calculate Jaccard Index (similarity between two gene sets)

# VARIABLES:
# x: Gene set #1
# y: Gene set #2

JI <- function(x,y){length(intersect(x,y))/length(union(x,y))}

#==============================================================================#
# makeORAnetwork()
#==============================================================================#

# DESCRIPTION:
# Make network plot of the results from the gene overrepresentation anlysis (ORA)

# VARIABLES:
# ORA_data: ORA data object (clusterProfiler object)
# layout: graph layout
# nSets: number of gene sets (i.e., GO terms/WikiPathways) to include in the plot

makeORAnetwork <- function(ORA_data, layout, nSets, color, size = 5){
  
  # round the number of sets to an integer value
  nSets <- round(nSets)
  
  # Set the minimum number to 5
  if(nSets < 5){
    nSets <- 5
  }
  
  # Set the maximum number to 20
  if(nSets > 20){
    nSets <- 20
  }
  
  # Select color gradient
  gradient = ggplot2::scale_color_viridis_c()
  
  if (color == "Yellow-red"){
    gradient = ggplot2::scale_color_gradient(low = "yellow", high = "red")
  }
  if (color == "Blues"){
    gradient = ggplot2::scale_color_gradient(low = "#DEEBF7", high = "#08519C")
  }
  if (color == "Reds"){
    gradient = ggplot2::scale_color_gradient(low = "#FEE0D2", high = "#CB181D")
  }
  
  # Get results from ORA
  ORAresults <- ORA_data@result
  
  # Filter genes for top 10 terms
  ORAgenes <- ORA_data@geneSets
  ORAgenes_fil <- ORAgenes[ORAresults$ID[1:nSets]]
  
  # Make a matrix that shows pairwise Jaccard Index
  graph_matrix <- matrix(NA, nrow = length(ORAgenes_fil), ncol = length(ORAgenes_fil))
  for (i in 1:length(ORAgenes_fil)){
    graph_matrix[i,] <- unlist(lapply(ORAgenes_fil,function(x){JI(x,ORAgenes_fil[[i]])}))
  }
  rownames(graph_matrix) <- names(ORAgenes_fil)
  colnames(graph_matrix) <- names(ORAgenes_fil)
  
  # make a graph from this matrix
  g <- igraph::graph_from_adjacency_matrix(graph_matrix, 
                                           mode = "lower", 
                                           weighted = "Jaccard Index", 
                                           diag = FALSE)
  
  # Add -log10 p-value and GO name as vertex attributes
  rownames(ORAresults) <- ORAresults$ID
  V(g)$`-log10 p-value` <- (-log10(ORAresults[V(g),"p-value"]))
  V(g)$label <- firstup(ORAresults[V(g),"Description"])
  
  # make plot
  p <- ggraph::ggraph(g, 'igraph', algorithm = layout) +
    ggraph::geom_edge_link0(ggplot2::aes(width = `Jaccard Index`), edge_alpha = 0.1) + 
    ggraph::geom_node_point(ggplot2::aes(color = `-log10 p-value`), size = 7) + 
    ggraph::geom_node_text(ggplot2::aes(label = label), color = 'black', 
                           size = size, repel = TRUE) + 
    ggplot2::theme_void() +
    ggplot2::labs(color = expression(-log[10] ~ "p-value")) +
    gradient
    #ggplot2::scale_color_manual(values = gradient) 
    #ggplot2::scale_color_gradient(low = "#6BAED6", high = "#FB6A4A")
  
  return(p)
}

#==============================================================================#
# makeGSEAnetwork()
#==============================================================================#

# DESCRIPTION:
# Make network plot of the results from the geneset enrichment analysis (GSEA)

# VARIABLES:
# GSEA_data: GSEA data object (clusterProfiler object)
# layout: graph layout
# nSets: number of gene sets (i.e., GO terms/WikiPathways) to include in the plot

makeGSEAnetwork <- function(GSEA_data, layout, nSets, color, size = 5){
  
  # round the number of sets to an integer value
  nSets <- round(nSets)
  
  # Set the minimum number to 5
  if(nSets < 5){
    nSets <- 5
  }
  
  # Set the maximum number to 20
  if(nSets > 20){
    nSets <- 20
  }
  
  if (is.null(color)){
    color <- c("#000072", "white", "red")
  }
  
  # Get results from GSEA
  GSEAresults <- GSEA_data@result
  
  # Filter genes for top 10 terms
  GSEAgenes <- GSEA_data@geneSets
  GSEAgenes_fil <- GSEAgenes[GSEAresults$ID[1:nSets]]
  
  # Make a matrix that shows pairwise Jaccard Index
  graph_matrix <- matrix(NA, nrow = length(GSEAgenes_fil), ncol = length(GSEAgenes_fil))
  for (i in 1:length(GSEAgenes_fil)){
    graph_matrix[i,] <- unlist(lapply(GSEAgenes_fil,function(x){JI(x,GSEAgenes_fil[[i]])}))
  }
  rownames(graph_matrix) <- names(GSEAgenes_fil)
  colnames(graph_matrix) <- names(GSEAgenes_fil)
  
  # make a graph from this matrix
  g <- igraph::graph_from_adjacency_matrix(graph_matrix, 
                                           mode = "lower", 
                                           weighted = "Jaccard Index", 
                                           diag = FALSE)
  
  # Add -log10 p-value and GO name as vertex attributes
  rownames(GSEAresults) <- GSEAresults$ID
  V(g)$`-log10 p-value` <- (-log10(GSEAresults[V(g),"p-value"]))
  V(g)$`NES` <- GSEAresults[V(g),"NES"]
  V(g)$label <- firstup(GSEAresults[V(g),"Description"])
  
  # make plot
  p <- ggraph::ggraph(g, 'igraph', algorithm = layout) +
    ggraph::geom_edge_link0(ggplot2::aes(width = `Jaccard Index`), edge_alpha = 0.1) + 
    ggraph::geom_node_point(ggplot2::aes(color = NES), size = 7) + 
    ggraph::geom_node_text(ggplot2::aes(label = label), color = 'black', 
                           size = size, repel = TRUE) + 
    ggplot2::theme_void() +
    ggplot2::labs(color = "NES") +
    #ggplot2::scale_color_continuous() +
    ggplot2::scale_color_gradient2(low = color[1], mid = color[2], high = color[3])
  
  return(p)
}

#==============================================================================#
# make_unique_names()
#==============================================================================#
make_unique_names <- function(names_vector) {
  table_counts <- table(names_vector)
  name_counter <- setNames(rep(0, length(table_counts)), names(table_counts))
  
  sapply(names_vector, function(name) {
    name_counter[name] <<- name_counter[name] + 1
    if (table_counts[name] > 1) {
      paste0(name, "_", name_counter[name])
    } else {
      name
    }
  })
}

#==============================================================================#
# readExprMatrix()
#==============================================================================#

# DESCRIPTION:
# Read expression matrix

# VARIABLES:
# path: path to expression matrix

readExprMatrix <- function(path, checkDup = TRUE){
  
  # Read matrix
  exprMatrix <- as.matrix(data.table::fread(path))
  
  # Set rownames
  rownames(exprMatrix) <- exprMatrix[,1]
  exprMatrix <- exprMatrix[,-1]
  
  # Make numeric
  exprMatrix_num <- matrix(as.numeric(exprMatrix), nrow = nrow(exprMatrix), ncol = ncol(exprMatrix))
  
  # Set rownames and column names of numeric matrix
  
  if (checkDup){
    if (sum(duplicated(rownames(exprMatrix)))>0){
      rownames(exprMatrix_num) <- make_unique_names(rownames(exprMatrix)) 
    } else{
      rownames(exprMatrix_num) <- rownames(exprMatrix)
    }
  } else{
    rownames(exprMatrix_num) <- rownames(exprMatrix)
  }
  colnames(exprMatrix_num) <- colnames(exprMatrix)
  
  # Return expression matrix
  return(exprMatrix_num)
}

#==============================================================================#
# RNASeqNormalization()
#==============================================================================#

# DESCRIPTION:
# Perform filtering and normalization (DESeq2) of RNA-seq data

# VARIABLES:
# gxData: matrix of raw counts
# metaData: meta data
# filterThreshold: minimum number of counts for n = smallestGroupSize
# smallestGroupSize: size of smallest experimental group

RNASeqNormalization <- function(gxData,
                                metaData,
                                filterThres = 10,
                                smallestGroupSize){
  
  # Generate DESeqDataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(gxData),
                                        colData = metaData,
                                        design = as.formula("~ experimentFactor"))
  
  # Filtering
  keep <- rowSums(counts(dds) >= filterThres) >= smallestGroupSize
  dds <- dds[keep,]
  
  # estimate size factors (for normalization)
  dds <- BiocGenerics::estimateSizeFactors(dds)
  
  # Return normalized counts
  normalized_counts <- list()
  normalized_counts[["counts"]] <- log2(counts(dds, normalized=TRUE)+1)
  normalized_counts[["vst"]] <- assay(DESeq2::vst(dds, blind = FALSE))
  return(normalized_counts)
}

#==============================================================================#
# RNASeqNormalization_processed()
#==============================================================================#

# DESCRIPTION:
# Additional processing of processed data (ExpressionSet)

# VARIABLES:
# gxData: matrix of normalized counts
# experimentFactor: factor of experimental group
# transMeth: transformation method (e.g., "None", "Log2-transformation")
# normMeth: normalization method (e.g., "None", "Quantile")
# perGroup_name: perform normalization per group?

RNASeqNormalization_processed <- function(gxMatrix,
                                          experimentFactor,
                                          transMeth,
                                          normMeth,
                                          perGroup_name,
                                          filterThres){
  
  # Get expression values from data object
  m <- gxMatrix
  
  # Filtering
  smallestGroupSize <- min(table(experimentFactor))
  keep <- rowSums(m >= filterThres) >= smallestGroupSize
  m <- m[keep,]
  
  # Perform log2-transformation (if selected)
  if (transMeth == "Log2-transformation"){
    m <- m[rowSums(m<=0)==0,]
    m <- log2(m+1)
  }
  
  # No additional normalization
  if(normMeth == "None"){
    gxMatrix_final <-  m
  }
  
  # Quantile normalization
  if(normMeth == "Quantile"){
    
    # Use all array for normalization
    if (perGroup_name == "Use all samples"){
      m <- limma::normalizeQuantiles(m)
      gxMatrix_final <- m
    } else{
      for (g in levels(experimentFactor)){
        m[,experimentFactor == g] <- limma::normalizeQuantiles(m[,experimentFactor == g])
      }
      # Perform normalization per experimental group
      gxMatrix_final <- m
    }
  }
  return(gxMatrix_final)
}

#==============================================================================#
# getStatistics_RNASeq()
#==============================================================================#

# DESCRIPTION:
# Perform statistical analysis (DESeq2) and optionally gene annotation (biomaRt)

# VARIABLES:
# rawMatrix: matrix of raw counts
# metaData: meta data
# expFactor: column name(s) of the meta data that contain the experimental factor
# covGroups_num: continuous/numerical covariates
# covGroups_char: discrete/character covariates
# comparisons: for which which comparisons to calculate statistics
# filterThreshold: minimum number of counts for n = smallestGroupSize
# smallestGroupSize: size of smallest experimental group
# addAnnotation: add biomaRt annotations to the output table
# biomart_dataset: biomaRt dataset (e.g., hsapiens_gene_ensembl)
# biomart_attributes: biomaRt attributes (output IDs)
# biomart_filters: biomaRt filters (input IDs)

getStatistics_RNASeq <- function(rawMatrix, 
                                 metaData, 
                                 expFactor, 
                                 covGroups_num, 
                                 covGroups_char, 
                                 comparisons,
                                 filterThres = 10,
                                 smallestGroupSize,
                                 shrinkage = TRUE,
                                 addAnnotation = TRUE,
                                 biomart_dataset = "hsapiens_gene_ensembl",
                                 biomart_attributes = c("Ensembl Gene ID",
                                                        "Entrez Gene ID",
                                                        "Gene Symbol/Name"),
                                 biomart_filters = "Entrez Gene ID"){
  tryCatch({
  dataset <- NULL
  
  # Replace name of biomaRt filter
  biomart_filters <- tryCatch({
    switch(biomart_filters,
           "Gene Symbol/Name" = "gene_name",
           "Entrez Gene ID" = "entrezgene_id",
           "Ensembl Gene ID" = "ensembl_gene_id",
    )
  }, error = function(cond){
    return(biomart_filters)
  })
  
  # Replace name(s) of biomaRt attributes
  for (a in 1:length(biomart_attributes)){
    biomart_attributes[a] <- tryCatch({
      switch(biomart_attributes[a],
             "Gene Symbol/Name" = "gene_name",
             "Entrez Gene ID" = "entrezgene_id",
             "Ensembl Gene ID" = "ensembl_gene_id",
      )
    }, error = function(cond){
      return(biomart_attributes[a])
    })
  }
  
  # Get experiment factor
  if(length(expFactor) > 1){
    experimentFactor <- factor(make.names(apply(metaData[,expFactor], 1, paste, collapse = "_" )))
  } else{
    experimentFactor <- factor(make.names(metaData[,expFactor]))
  }
  
  # Get covariates
  for (n in covGroups_num){
    metaData[,n] <- as.numeric(metaData[,n])
  }
  for (c in covGroups_char){
    metaData[,c] <- factor(metaData[,c])
  }
  covariates <- c(covGroups_num, covGroups_char)
  
  # Make formula
  if (!is.null(covariates)){
    formula <- paste0("~ ", "experimentFactor + ", paste0(covariates, collapse = " + "))
  } else {
    formula <- paste0("~ ", "experimentFactor")
  }
  
  # make DESeqDataSet object
  sampleInfo <- cbind.data.frame(experimentFactor, metaData[,covariates])
  colnames(sampleInfo) <- c("experimentFactor", covariates)
  rownames(sampleInfo) <- rownames(metaData)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(rawMatrix),
                                        colData = sampleInfo,
                                        design = as.formula(formula))
  # Filtering
  keep <- rowSums(counts(dds) >= filterThres) >= smallestGroupSize
  dds <- dds[keep,]
  
  
  # Perform statistical comparison for each of the selected comparison
  top_table <- list()
  for (i in comparisons){
    
    # Change level
    referenceLevel <- make.names(stringr::str_split(i," - ")[[1]][2])
    dds$experimentFactor <- relevel(dds$experimentFactor, ref = referenceLevel)
    dds <- DESeq2::DESeq(dds)
    
    contrastName <- paste("experimentFactor",
                          make.names(stringr::str_split(i," - ")[[1]][1]),
                          "vs",
                          make.names(stringr::str_split(i," - ")[[1]][2]), sep = "_")
    
    if (shrinkage){
      resLFC <- as.data.frame(DESeq2::lfcShrink(dds, coef=contrastName, type="apeglm"))
      resLFC <- dplyr::arrange(resLFC,by = pvalue)
    } else{
      resLFC <- as.data.frame(DESeq2::results(dds, name=contrastName))
      resLFC <- dplyr::arrange(resLFC,by = pvalue)
    }
    
    top_table[[i]] <- resLFC[,c("baseMean",
                                "log2FoldChange",
                                "lfcSE",
                                "pvalue",
                                "padj")]
    top_table[[i]] <- cbind(rownames(top_table[[i]]), top_table[[i]])
    rownames(top_table[[i]]) <- NULL
    colnames(top_table[[i]]) <- c("GeneID", "meanExpr", "log2FC", "log2FC SE",
                                  "p-value", "adj. p-value")
  }
  
  message <- "Statistical analysis has been performed. 
    You can download the results as well as view them in interactive plots."
  
  # Add annotations to table if this option is selected
  if (addAnnotation == TRUE){
    
    # Change attribute and filter name
    if (biomart_dataset == "hsapiens_gene_ensembl"){
      biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                  "gene_name",
                                                  "hgnc_symbol")
      biomart_filters1 <- stringr::str_replace(biomart_filters,
                                               "gene_name",
                                               "hgnc_symbol")
    } else{
      biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                  "gene_name",
                                                  "external_gene_name")
      biomart_filters1 <- stringr::str_replace(biomart_filters,
                                               "gene_name",
                                               "external_gene_name")
    }
    
    
    # Do this for each top table in the list
    for (t in 1:length(top_table)){
      
      # Round numbers in top table
      for (n in 2:6){
        top_table[[t]][,n] <- signif(top_table[[t]][,n],3)
      }
      
      #Get annotations
      annotation_list <- tryCatch({
        ensembl <- biomaRt::useMart("ensembl")
        ensembl <- biomaRt::useDataset(biomart_dataset, mart=ensembl)
        annotations <- biomaRt::getBM(attributes=biomart_attributes1,
                                      filters = biomart_filters1,
                                      values = top_table[[t]]$GeneID,
                                      mart = ensembl)
        
        message <- "Statistical analysis has been performed. 
          Gene annotation was performed with biomaRt. You can download 
              the results as well as view them in interactive plots."
        dataset <- paste0(biomart_dataset, " (", searchDatasets(mart = ensembl, pattern = "hsapiens")$version, ")")
        list(annotations, message, dataset)
      },
      error = function(cond){
        
        # Load annotation package
        pkg <- switch(biomart_dataset,
                      "hsapiens_gene_ensembl" = "org.Hs.eg.db",
                      "btaurus_gene_ensembl" = "org.Bt.eg.db",
                      "celegans_gene_ensembl" = "org.Ce.eg.db",
                      "mmusculus_gene_ensembl" = "org.Mm.eg.db",
                      "rnorvegicus_gene_ensembl" = "org.Rn.eg.db"
        )
        
        if (!requireNamespace(pkg, quietly = TRUE))
          BiocManager::install(pkg, ask = FALSE)
        require(as.character(pkg), character.only = TRUE)
        
        
        # Change attribute names
        biomart_attributes2 <- biomart_attributes
        if (biomart_dataset == "hsapiens_gene_ensembl"){
          biomart_attributes2[biomart_attributes2 == "gene_name"] <- "SYMBOL"
        } else{
          biomart_attributes2[biomart_attributes2 == "gene_name"] <- "GENENAME"
        }
        biomart_attributes2[biomart_attributes2 == "entrezgene_id"] <- "ENTREZID"
        biomart_attributes2[biomart_attributes2 == "ensembl_gene_id"] <- "ENSEMBL"
        
        
        # Change filter names
        biomart_filters2 <- biomart_filters
        if (biomart_dataset == "hsapiens_gene_ensembl"){
          biomart_filters2[biomart_filters2 == "gene_name"] <- "SYMBOL"
        } else{
          biomart_filters2[biomart_filters2  == "gene_name"] <- "GENENAME"
        }
        biomart_filters2[biomart_filters2 == "entrezgene_id"] <- "ENTREZID"
        biomart_filters2[biomart_filters2 == "ensembl_gene_id"] <- "ENSEMBL"
        
        # Get annotations
        annotations <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                             columns = biomart_attributes2, 
                                             keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
        
        # Join with geneIDs
        annotations <- left_join(data.frame(GeneID = top_table[[t]]$GeneID),
                                 annotations,
                                 by = c("GeneID" = biomart_filters2))
        
        # Change column names back
        temp_col <- colnames(annotations)
        temp_col[temp_col == "GeneID"] <- biomart_filters1
        if (biomart_dataset == "hsapiens_gene_ensembl"){
          temp_col[temp_col == "SYMBOL"] <- "hgnc_symbol"
        } else{
          temp_col[temp_col == "SYMBOL"] <- "external_gene_name"
        }
        temp_col[temp_col == "ENTREZID"] <- "entrezgene_id"
        temp_col[temp_col == "ENSEMBL"] <- "ensembl_gene_id"
        colnames(annotations) <- temp_col
        message <- "Statistical analysis has been performed. 
          The Ensembl database was not available.
          So, the gene annotation was performed with the bioconductor annotation package. 
          You can download the results as well as view them in interactive plots."
        dataset <- paste0(pkg, " (", packageVersion(pkg),")")
        list(annotations, message, dataset)
      })
      annotations <- annotation_list[[1]]
      message <- annotation_list[[2]]
      dataset <- annotation_list[[3]]
      
      # Convert entrezgene id to character
      if("entrezgene_id" %in% biomart_attributes){
        annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
      }
      annotations[annotations == ""] <- NA
      annotations[annotations == " "] <- NA
      
      # Combine annotations with top table
      #annotations[,biomart_filters] <- as.character(annotations[,biomart_filters1])
      top_table_ann <- dplyr::left_join(top_table[[t]], annotations,
                                        by = c("GeneID" = biomart_filters1))
      
      colnames(top_table_ann)[colnames(top_table_ann) == "hgnc_symbol"] <- "gene_name"
      colnames(top_table_ann)[colnames(top_table_ann) == "external_gene_name"] <- "gene_name"
      
      # Make sure that there are no duplicate gene ids
      for (a in 1:(ncol(top_table_ann)-6)){
        temp1 <- unique(top_table_ann[,c(1,6+a)])
        temp1 <- temp1[!is.na(temp1[,2]),]
        temp_ann <- temp1 %>%
          dplyr::group_by(GeneID) %>%
          dplyr::summarise_at(colnames(top_table_ann)[6+a], function(x) paste(x,collapse = "; "))
        
        top_table[[t]] <- dplyr::left_join(top_table[[t]], temp_ann,
                                           by = c("GeneID" = "GeneID"))
      }
      
    }
    
  }
  top_table_list <- list(top_table, message, dataset)
  return(top_table_list)
  }, error = function(cond){
    NULL
  })
}



#==============================================================================#
# getStatistics_RNASeq_processed()
#==============================================================================#

# DESCRIPTION:
# Perform statistical analysis (limma) and optionally probeset annotation (biomaRt)

# VARIABLES:
# normMatrix: normalized intensity matrix
# metaData: meta data
# expFactor: column name(s) of the meta data that contain the experimental factor
# covGroups_num: continuous/numerical covariates
# covGroups_char: discrete/character covariates
# addAnnotation: add biomaRt annotations to the output table
# biomart_dataset: biomaRt dataset (e.g., hsapiens_gene_ensembl)
# biomart_attributes: biomaRt attributes (output IDs)
# biomart_filters: biomaRt filters (input IDs)

getStatistics_RNASeq_processed <- function(normMatrix, 
                                           metaData, 
                                           expFactor, 
                                           covGroups_num, 
                                           covGroups_char, 
                                           comparisons,
                                           addAnnotation = TRUE,
                                           biomart_dataset = "hsapiens_gene_ensembl",
                                           biomart_attributes = c("Ensembl Gene ID",
                                                                  "Entrez Gene ID",
                                                                  "Gene Symbol/Name"),
                                           biomart_filters = "Entrez Gene ID"){
  tryCatch({
    dataset <- NULL
    
    # Replace name of biomaRt filter
    biomart_filters <- tryCatch({
      switch(biomart_filters,
             "Gene Symbol/Name" = "gene_name",
             "Entrez Gene ID" = "entrezgene_id",
             "Ensembl Gene ID" = "ensembl_gene_id",
      )
    }, error = function(cond){
      return(biomart_filters)
    })
    
    # Replace name(s) of biomaRt attributes
    for (a in 1:length(biomart_attributes)){
      biomart_attributes[a] <- tryCatch({
        switch(biomart_attributes[a],
               "Gene Symbol/Name" = "gene_name",
               "Entrez Gene ID" = "entrezgene_id",
               "Ensembl Gene ID" = "ensembl_gene_id",
        )
      }, error = function(cond){
        return(biomart_attributes[a])
      })
    }
    
    # Get experiment factor
    if(length(expFactor) > 1){
      experimentFactor <- factor(make.names(apply(metaData[,expFactor], 1, paste, collapse = "_" )))
    } else{
      experimentFactor <- factor(make.names(metaData[,expFactor]))
    }
    
    # Get covariates
    for (n in covGroups_num){
      metaData[,n] <- as.numeric(metaData[,n])
    }
    for (c in covGroups_char){
      metaData[,c] <- as.character(metaData[,c])
    }
    covariates <- c(covGroups_num, covGroups_char)
    
    # Make formula
    if (!is.null(covariates)){
      formula <- paste("~ 0", "experimentFactor", paste0("metaData$`",covariates,"`", collapse = "+"),sep = "+")
    } else {
      formula <- paste("~ 0", "experimentFactor",sep = "+")
    }
    
    # Make design matrix
    design <- model.matrix(as.formula(formula))
    colnames(design)[1:length(levels(experimentFactor))] <- make.names(levels(factor(experimentFactor)))
    colnames(design) <- make.names(colnames(design))
    
    
    #fit linear model
    data.fit <- limma::lmFit(normMatrix,design)
    
    # Perform statistical comparison for each of the selected comparison
    top_table <- list()
    for (i in comparisons){
      
      # Make contrast
      contrast.matrix <- limma::makeContrasts(contrasts = i,levels=design)
      
      # Fit contrasts
      data.fit.con <- limma::contrasts.fit(data.fit,contrast.matrix)
      
      # eBayes (limma trend)
      data.fit.eb <- limma::eBayes(data.fit.con, trend = TRUE)
      
      temp_toptable <- limma::topTable(data.fit.eb, 
                                       sort.by = "P",
                                       n = Inf,
                                       confint = TRUE)
      
      # Calculate SE of log2FC
      temp_toptable$SE <- (temp_toptable$CI.R - temp_toptable$CI.L)/3.92
      
      # get top table
      top_table[[i]] <- temp_toptable[,c("AveExpr",
                                         "logFC",
                                         "SE",
                                         "P.Value",
                                         "adj.P.Val")]
      top_table[[i]] <- cbind(rownames(top_table[[i]]), top_table[[i]])
      rownames(top_table[[i]]) <- NULL
      colnames(top_table[[i]]) <- c("GeneID", "meanExpr", "log2FC", "log2FC SE",
                                    "p-value", "adj. p-value")
    }
    message <- "Statistical analysis has been performed. 
    You can download the results as well as view them in interactive plots."
    
    # Add annotations to table if this option is selected
    if (addAnnotation == TRUE){
      
      # Change attribute and filter name
      if (biomart_dataset == "hsapiens_gene_ensembl"){
        biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                    "gene_name",
                                                    "hgnc_symbol")
        biomart_filters1 <- stringr::str_replace(biomart_filters,
                                                 "gene_name",
                                                 "hgnc_symbol")
      } else{
        biomart_attributes1 <- stringr::str_replace(biomart_attributes,
                                                    "gene_name",
                                                    "external_gene_name")
        biomart_filters1 <- stringr::str_replace(biomart_filters,
                                                 "gene_name",
                                                 "external_gene_name")
      }
      
      
      # Do this for each top table in the list
      for (t in 1:length(top_table)){
        
        # Round numbers in top table
        for (n in 2:6){
          top_table[[t]][,n] <- signif(top_table[[t]][,n],3)
        }
        
        #Get annotations
        annotation_list <- tryCatch({
          ensembl <- biomaRt::useMart("ensembl")
          ensembl <- biomaRt::useDataset(biomart_dataset, mart=ensembl)
          annotations <- biomaRt::getBM(attributes=biomart_attributes1,
                                        filters = biomart_filters1,
                                        values = top_table[[t]]$GeneID,
                                        mart = ensembl)
          
          message <- "Statistical analysis has been performed. 
          Gene annotation was performed with biomaRt. You can download 
              the results as well as view them in interactive plots."
          dataset <- paste0(biomart_dataset, " (", searchDatasets(mart = ensembl, pattern = "hsapiens")$version, ")")
          list(annotations, message, dataset)
        },
        error = function(cond){
          
        # Load annotation package
          pkg <- switch(biomart_dataset,
                        "hsapiens_gene_ensembl" = "org.Hs.eg.db",
                        "btaurus_gene_ensembl" = "org.Bt.eg.db",
                        "celegans_gene_ensembl" = "org.Ce.eg.db",
                        "mmusculus_gene_ensembl" = "org.Mm.eg.db",
                        "rnorvegicus_gene_ensembl" = "org.Rn.eg.db"
          )
          
          if (!requireNamespace(pkg, quietly = TRUE))
            BiocManager::install(pkg, ask = FALSE)
          require(as.character(pkg), character.only = TRUE)
          
          
          # Change attribute names
          biomart_attributes2 <- biomart_attributes
          if (biomart_dataset == "hsapiens_gene_ensembl"){
            biomart_attributes2[biomart_attributes2 == "gene_name"] <- "SYMBOL"
          } else{
            biomart_attributes2[biomart_attributes2 == "gene_name"] <- "GENENAME"
          }
          biomart_attributes2[biomart_attributes2 == "entrezgene_id"] <- "ENTREZID"
          biomart_attributes2[biomart_attributes2 == "ensembl_gene_id"] <- "ENSEMBL"
          
          
          # Change filter names
          biomart_filters2 <- biomart_filters
          if (biomart_dataset == "hsapiens_gene_ensembl"){
            biomart_filters2[biomart_filters2 == "gene_name"] <- "SYMBOL"
          } else{
            biomart_filters2[biomart_filters2  == "gene_name"] <- "GENENAME"
          }
          biomart_filters2[biomart_filters2 == "entrezgene_id"] <- "ENTREZID"
          biomart_filters2[biomart_filters2 == "ensembl_gene_id"] <- "ENSEMBL"
        
          # Get annotations
          annotations <- AnnotationDbi::select(BiocGenerics::get(pkg), 
                                            columns = biomart_attributes2, 
                                            keys = AnnotationDbi::keys(BiocGenerics::get(pkg)))
          
          # Join with geneIDs
          annotations <- left_join(data.frame(GeneID = top_table[[t]]$GeneID),
                                   annotations,
                                   by = c("GeneID" = biomart_filters2))
          
          # Change column names back
          temp_col <- colnames(annotations)
          temp_col[temp_col == "GeneID"] <- biomart_filters1
          if (biomart_dataset == "hsapiens_gene_ensembl"){
          temp_col[temp_col == "SYMBOL"] <- "hgnc_symbol"
          } else{
            temp_col[temp_col == "SYMBOL"] <- "external_gene_name"
          }
          temp_col[temp_col == "ENTREZID"] <- "entrezgene_id"
          temp_col[temp_col == "ENSEMBL"] <- "ensembl_gene_id"
          colnames(annotations) <- temp_col
          message <- "Statistical analysis has been performed. 
          The Ensembl database was not available.
          So, the gene annotation was performed with the bioconductor annotation package. 
          You can download the results as well as view them in interactive plots."
          dataset <- paste0(pkg, " (", packageVersion(pkg),")")
          list(annotations, message, dataset)
        })
        annotations <- annotation_list[[1]]
        message <- annotation_list[[2]]
        dataset <- annotation_list[[3]]

        # Convert entrezgene id to character
        if("entrezgene_id" %in% biomart_attributes){
          annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
        }
        annotations[annotations == ""] <- NA
        annotations[annotations == " "] <- NA
        
        # Combine annotations with top table
        #annotations[,biomart_filters] <- as.character(annotations[,biomart_filters1])
        top_table_ann <- dplyr::left_join(top_table[[t]], annotations,
                                          by = c("GeneID" = biomart_filters1))
        
        colnames(top_table_ann)[colnames(top_table_ann) == "hgnc_symbol"] <- "gene_name"
        colnames(top_table_ann)[colnames(top_table_ann) == "external_gene_name"] <- "gene_name"
        
        # Make sure that there are no duplicate gene ids
        for (a in 1:(ncol(top_table_ann)-6)){
          temp1 <- unique(top_table_ann[,c(1,6+a)])
          temp1 <- temp1[!is.na(temp1[,2]),]
          temp_ann <- temp1 %>%
            dplyr::group_by(GeneID) %>%
            dplyr::summarise_at(colnames(top_table_ann)[6+a], function(x) paste(x,collapse = "; "))
          
          top_table[[t]] <- dplyr::left_join(top_table[[t]], temp_ann,
                                             by = c("GeneID" = "GeneID"))
        }
        
      }
      
    }
    top_table_list <- list(top_table, message, dataset)
    return(top_table_list)
  }, error = function(cond){
      NULL
  })
}


#==============================================================================#
# checkTransformation()
#==============================================================================#

# DESCRIPTION:
# Check if transformation is required

# VARIABLES:
# gxMatrix: expression matrix
# RNASeq data?

checkTransformation <- function(gxMatrix, RNASeq = TRUE){
  
  if (isTRUE(RNASeq)){
    gxMatrix <- gxMatrix[rowSums(gxMatrix<=0)==0,]
    r <- range(gxMatrix)
    diff <- r[2]-r[1]
    if (diff > 500){
      output <- "Count matrix"
    } else{
      output <- "log-transformed matrix"
    }
  } else{
    gxMatrix <- gxMatrix[rowSums(gxMatrix<=0)==0,]
    r <- range(gxMatrix)
    diff <- r[2]-r[1]
    if (diff > 500){
      output <- "Count matrix"
    } else{
      output <- "log-transformed matrix"
    }
  }
  return(output)
}

#==============================================================================#
# checkTransformation()
#==============================================================================#

# DESCRIPTION:
# Check if filtering is required

# VARIABLES:
# gxMatrix: expression matrix

checkFiltering <- function(gxMatrix){
  
  
  n0 <- sum(rowSums(gxMatrix==0)==ncol(gxMatrix))
  if (n0 > 0){
    output <- "No filtering"
  } else{
    output <- "Filtering"
  }
  
  return(output)
}

#==============================================================================#
# permute()
#==============================================================================#

permute <- function(vec) {
  if (length(vec) == 1) {
    return(list(vec))
  }
  
  result <- list()
  for (i in seq_along(vec)) {
    rest <- vec[-i]
    sub_permutations <- permute(rest)
    for (sub in sub_permutations) {
      result <- append(result, list(c(vec[i], sub)))
    }
  }
  return(result)
}