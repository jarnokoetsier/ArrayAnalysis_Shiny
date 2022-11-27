

################################################################################

#Install and load packages

################################################################################


#******************************************************************************#
#CRAN packages
#******************************************************************************#

#Required CRAN packages:
CRANpackages <- c("tidyverse", 
                  "fuzzyjoin", 
                  "plotly", 
                  "gplots", 
                  "heatmaply",
                  "ggvenn",
                  "DT",
                  "shiny",
                  "shinyWidgets",
                  "shinycssloaders")

#Install (if not yet installed) and load the required packages: 
for (pkg in CRANpackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}



#******************************************************************************#
#Bioconductor packages
#******************************************************************************#

#Required Bioconductor packages:
BiocPackages <- c("biomaRt",
                  "GEOquery",
                  "ArrayExpress",
                  "limma",
                  "makecdfenv",
                  "topGO",
                  "hugene10stprobeset.db",
                  "bioDist",
                  "gcrma",
                  "plier",
                  "clusterProfiler",
                  "enrichplot",
                  "KEGGREST")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (pkg in BiocPackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}


#******************************************************************************#
#Brainarray packages
#******************************************************************************#

#Modified affy package from Brainarray
if(!('affy' %in% installed.packages()[, "Package"])) {
  install.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/affy_1.68.0.tar.gz", 
                   type = "source", repos = NULL)
}
library(affy)


################################################################################

#uploadcdfenv

################################################################################

uploadcdfenv <- function(Data,cdf_path){
  #initial value
  CDFenv <- 0
  
  # recall which cdfName was added, in case no updated one is found (set back
  # even if it does not exist)
  presetCDF <- Data@cdfName
  
  #try to load cdf file
  try(CDFenv <- make.cdf.env(filename = basename(cdf_path),cdf.path=dirname(cdf_path)),TRUE)
  
  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("Could not load custom CDF environment for this chip type - object kept as is")
  }
  
  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}


################################################################################

#addUpdatedCDFenv

################################################################################

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
  suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  #try without a version number
  print(Data@cdfName<-paste(gsub("v[0-9]$","",Data@annotation),species,type,sep="_"))
  suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE)) 
  
  #if it hasn't loaded, try to download
  if ((class(CDFenv)!="environment")) {
    install.packages(tolower(paste(Data@annotation,species,type,"cdf",sep="")),
                     repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  }
  
  #if it hasn't loaded, try to download without version number
  if ((class(CDFenv)!="environment")) {
    install.packages(tolower(paste(gsub("v[0-9]$","",Data@annotation),species,type,"cdf",sep="")),
                     repos="http://brainarray.mbni.med.umich.edu/bioc")
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
  }
  
  if ((class(CDFenv)!="environment")) {
    Data@cdfName <- presetCDF
    warning("Could not automatically retrieve CDF environment for this chip type - object kept as is")
  }
  
  cat("current cdf environment loaded:",Data@cdfName,"\n")
  return(Data)
}


################################################################################

#colorsByFactor

################################################################################


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

################################################################################

#get_grouping

################################################################################

get_grouping <- function(gset, database = "GEO") {
  
  if (database == "GEO"){
    GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
    
    grouping <- as.data.frame(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)[1]
    
    for (l in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)){
      if (((length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) > 1) & 
          ((length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) < length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))){
        
        grouping <- cbind(grouping, as.data.frame(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)[l])
      }
      
    }
    
    grouping <- as.data.frame(grouping[,-1])
    #Get character columns only
    grouping <- grouping[, sapply(grouping, class) == 'character']
    
    grouping <- grouping %>%
      mutate(across(everything(), as.character))
  }
  
  if (database == "ArrayExpress"){
    
    grouping <- as.data.frame(gset@phenoData@data)[1]
    
    for (l in 1:length(gset@phenoData@data)){
      if (((length(unique(gset@phenoData@data[[l]]))) > 1) & 
          ((length(unique(gset@phenoData@data[[l]]))) < length(gset@phenoData@data[[l]]))){
        
        grouping <- cbind(grouping, as.data.frame(gset@phenoData@data)[l])
      }
      
    }
    
    grouping <- as.data.frame(grouping[,-1])
    
    grouping <- grouping %>%
      mutate(across(everything(), as.character))
  }
  
  
  #remove columns with the same information
  uni = c(1)
  
  tryCatch({
    suppressWarnings(
      
      for (i in 1:(ncol(grouping)-1)) {
        for (j in (i+1):ncol(grouping)) {
          if (all(str_detect(grouping[,i], grouping[,j]))) {
            if (length(unique(str_remove(grouping[,i], grouping[,j]))) == 1){
              uni <- c(uni, i)
            }
          }
          if (all(str_detect(grouping[,j], grouping[,i]))) {
            if (length(unique(str_remove(grouping[,j], grouping[,i]))) == 1) {
              if (length(unique(str_remove(grouping[,i], grouping[,j]))) != 1){
                uni <- c(uni, j)
              }
            }
          }
        }
      }
    )
  },
  error=function(cond) {
    NULL
  })
  
  
  
  uni <- unique(uni[-1])
  if (length(uni)>0){
    grouping <- grouping[, -uni]
  }
  
  #rownames(grouping) <- NULL
  
  #minimum group size
  minGroup <- function(x){
    groupSize <- NULL
    for (i in 1:length(unique(x))){
      groupSize[i] <- sum(x == unique(x)[i])
    }
    return(min(groupSize))
  }
  
  #Remove columns which have minimum group size of 1
  test <- apply(grouping,2,minGroup)
  grouping <- grouping[,test>1]
  
  
  
  colnames(grouping) <- str_remove_all(colnames(grouping), ":ch1")
  colnames(grouping) <- str_remove_all(colnames(grouping), "_ch1")
  colnames(grouping) <- str_replace_all(colnames(grouping), "_", " ")
  colnames(grouping) <- str_to_title(colnames(grouping))
  
  return(grouping)
  
}


################################################################################

#auto_group

################################################################################

auto_group <- function(groups, attempt = 1){
  
  
  participant_group <- NULL
  
  n_ch <- NULL
  for (l in 1:length(groups)){
    if ((length(grep("control|non|healthy|treat", groups[[l]], ignore.case = TRUE)) > 0) &
        (length(unique(groups[[l]]))) < length(groups[[l]])){
      
      n_ch <- rbind(n_ch, c(l,nchar(groups[[l]][1])))
    }
  }
  
  if (!is.null(n_ch)){
    if (attempt <= nrow(n_ch)) {
      y = 1
      repeat {
        
        if (y == attempt) {
          break
        }
        y = y + 1
        n_ch = n_ch[-which.min(n_ch[,2]),]
      }
      
      int1 <- n_ch[which.min(n_ch[,2]),1]
      participant_group <- as.data.frame(groups)[int1]
      print(noquote("Grouping has been done automatically by default. Please check whether grouping has occured correctly."))
      
    }
    
    if (attempt > nrow(n_ch)) {
      print(noquote("grouping information not found. Upload grouping data manually at grouping_column"))
    }
    
  }
  
  
  if (is.null(n_ch)) {
    nlevels = NULL
    for (i in 1:ncol(groups)){
      nlevels <- c(nlevels, length(unique(groups[,i])))
    }
    
    participant_group <- as.data.frame(groups)[which.min(nlevels)]
  }
  
  return(colnames(participant_group))
}


################################################################################

#get_meta

################################################################################


get_meta <- function(gset, grouping_column, pairing_column = NULL, file_name = NULL, database = "GEO") {
  
  #Get participant info
  
  if(database == "GEO"){
    GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
    
    participants <- rownames(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)
  
  }
  
  if(database == "ArrayExpress"){
    
    participants <- rownames(gset@phenoData@data)
    
    
  }
  
  
  #use grouping column
  if (length(grouping_column[1]) > 0){
    if (!is.vector(grouping_column)){
      if (ncol(grouping_column) > 1) {
        participant_group <- unite(grouping_column, col = grouping_column, sep = ".")
      }
      if (ncol(grouping_column) == 1){
        participant_group <- grouping_column
      } 
    }
    if (is.vector(grouping_column)){
      participant_group <- grouping_column
    }
    
  }
  
  
  
  #independ samples
  if (is.null(pairing_column) | length(pairing_column) < 1) {
    if (exists("participant_group")){
      meta <- as.data.frame(cbind(participants, participant_group))
      rownames(meta) <- NULL
      colnames(meta) <- c("Sample.ID", "Grouping")
      return(meta)
    }
  }
  
  
  
  #dependent samples
  if (!is.null(pairing_column) & length(pairing_column) >= 1){
    if (!is.vector(pairing_column)){
      if (ncol(pairing_column) > 1) {
        pairs <- unite(pairing_column, col = pairing_column, sep = ".")
      }
      if (ncol(pairing_column) == 1){
        pairs <- pairing_column
      } 
    }
    if (is.vector(pairing_column)){
      pairs <- pairing_column
    }
    
    meta <- as.data.frame(cbind(participants, participant_group, pairs))
    rownames(meta) <- NULL
    colnames(meta) <- c("Sample.ID", "Grouping", "Pairing")
    return(meta)
  }
  
  
  
  #file upload
  if (!is.null(file_name)){
    meta <- read_excel(file_name)
    return(meta)
  }
}




################################################################################

#pca.plot

################################################################################


pca.plot <- function(data.PC, meta1, hpc = 1, vpc = 2) {
  
  #get PCs
  PC <- cbind.data.frame(rownames(data.PC$x), data.PC$x)
  colnames(PC) <- c("cel_names", colnames(PC[,2:ncol(PC)]))
  rownames(PC) <- NULL
  
  #Match samples with grouping using meta
  PC <- meta1 %>% inner_join(PC, by = c("names" = "cel_names"))
  
  PC$Grouping <- factor(PC$Grouping, levels = unique(PC$Grouping))
  PC$Sample.ID <- factor(PC$Sample.ID, levels = unique(PC$Sample.ID))
  
  
  #Make colours
  experimentFactor <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
  
  myPalette <- colorsByFactor(experimentFactor)

  #calculate explained variance
  perc_expl1 <- round(((data.PC$sdev^2)/sum(data.PC$sdev^2))*100,2)
  
  
  #Make plot
  if (is.null(meta1$Pairing)){
    X = PC[,(3+hpc)]
    Y = PC[,(3+vpc)]
  }
  
  if (!is.null(meta1$Pairing)){
    X = PC[,(4+hpc)]
    Y = PC[,(4+vpc)]
  }

  
  
  pcaplot <-
    ggplot(data = PC, aes(x = X, y = Y, colour = Grouping, shape = Grouping, 
                          text = paste("Sample:", Sample.ID))) +
    geom_point(size = 2) +
    scale_colour_manual(values = myPalette$legendColors) +
    labs(title = "PCA plot",
         subtitle = "Similar samples are clustered together") +
    xlab(paste0("PC",hpc, " (", perc_expl1[hpc], "%)")) +
    ylab(paste0("PC",vpc, " (", perc_expl1[vpc], "%)")) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.title = element_blank())
  
  return(ggplotly(pcaplot, tooltip = c("x", "y", "colour", "text")))
  
  
}

################################################################################

#pca3d

################################################################################

pca3d <- function(data.PC, meta1, x = 1, y = 2, z = 3) {

  #get PCs
  PC <- cbind.data.frame(rownames(data.PC$x), data.PC$x)
  colnames(PC) <- c("cel_names", colnames(PC[,2:ncol(PC)]))
  rownames(PC) <- NULL
  
  #Match samples with grouping using meta
  PC <- meta1 %>% inner_join(PC, by = c("names" = "cel_names"))
  
  PC$Grouping <- factor(PC$Grouping, levels = unique(meta1$Grouping))
  PC$Sample.ID <- factor(PC$Sample.ID, levels = unique(meta1$Sample.ID))
  
  #Make colours
  experimentFactor <- factor(meta1$Grouping, levels = unique(meta1$Grouping))
  
  myPalette <- colorsByFactor(experimentFactor)

  
  #calculate explained variance
  perc_expl1 <- round(((data.PC$sdev^2)/sum(data.PC$sdev^2))*100,2)
  
  #Make plot
  
  if (is.null(meta1$Pairing)){
    X = PC[,(3+x)]
    Y = PC[,(3+y)]
    Z = PC[,(3+z)]
  }
  
  if (!is.null(meta1$Pairing)){
    X = PC[,(4+x)]
    Y = PC[,(4+y)]
    Z = PC[,(4+z)]
  }
  
  
  pca3d <- plot_ly(x=X, 
                   y=Y, 
                   z=Z, 
                   type="scatter3d", 
                   mode="markers", 
                   color=PC$Grouping,
                   colors = myPalette$legendColors,
                   text = paste('Sample:', PC$Sample.ID, '<br>Grouping:', PC$Grouping))
  
  pca3d <- pca3d %>% layout(scene = list(xaxis = list(title = paste0('PC', x, " (", perc_expl1[x], "%)")),
                                         yaxis = list(title = paste0('PC', y, " (", perc_expl1[y], "%)")),
                                         zaxis = list(title = paste0('PC', z, " (", perc_expl1[z], "%)"))))
  
  
}


################################################################################

#PCvariances

################################################################################

PCvariances <- function(data.PC, x = 1, y = 2, z= NULL) {
  
  perc_expl1 <- round(((data.PC$sdev^2)/sum(data.PC$sdev^2))*100,2)
  
  
  selected <- rep("no", length(perc_expl1))
  selected[x] <- "yes"
  selected[y] <- "yes"
  selected[z] <- "yes"
  
  
  percentages <- cbind.data.frame(seq(1:length(perc_expl1)), perc_expl1)
  percentages <- cbind.data.frame(percentages, selected)
  colnames(percentages) <- c("PC", "Percentage.explained", "Selected")
  
  
  
  pcahist <-
    ggplot(data = percentages, aes(x = PC, y = Percentage.explained, fill = Selected,
                                   text = paste0("PC", PC, "<br>", "Explained variance: ", Percentage.explained, "%"))) +
    geom_bar(stat = "identity", colour = "#696969") +
    scale_fill_manual(values = c("#d3d3d3", "#FF4233")) +
    labs(title = "Histogram of explained variances") +
    xlab("Components") +
    ylab("% of the variance explained") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")
  
  return(ggplotly(pcahist, tooltip = "text"))
  
}

################################################################################

#get_contrasts

################################################################################

get_contrasts <- function(meta) {
  
  levels1 <- unique(meta$Grouping)
  
    contrast <- paste(make.names(levels1)[2], make.names(levels1)[1], sep = " - ")
    
    for (m in 1:(length(levels1)-1)){
      for (n in (m+1):length(levels1)) {
        
        if (grepl("non|healthy|control", levels1[m])) {
          contrast <- c(contrast, paste(make.names(levels1)[n], 
                                        make.names(levels1)[m], sep = " - "))
        }
        
        if (!grepl("non|healthy|control", levels1[m])){
          contrast <- c(contrast, paste(make.names(levels1)[m], 
                                        make.names(levels1)[n], sep = " - "))
        }
      }
    }
    return(sort(contrast[-1]))
  
}
  


################################################################################

#auto_contrast

################################################################################

auto_contrasts <- function(meta) {
  
  levels1 <- unique(meta$Grouping)
  
  contrast <- paste(make.names(levels1)[2], make.names(levels1)[1], sep = " - ")
  
  for (m in 1:(length(levels1)-1)){
    for (n in (m+1):length(levels1)) {
      
      if (grepl("non|healthy|control", levels1[m])) {
        contrast <- c(contrast, paste(make.names(levels1)[n], 
                                      make.names(levels1)[m], sep = " - "))
      }
      
      if (!grepl("non|healthy|control", levels1[m])){
        contrast <- c(contrast, paste(make.names(levels1)[m], 
                                      make.names(levels1)[n], sep = " - "))
      }
    }
  }
  
  
  contrast <- contrast[-1]
  
  contrast <- str_replace_all(contrast, "_", "\\.") 
  
  contrast1 <- contrast
  
  
  for (k in 1:length(contrast)) {
    
    
    a <- strsplit(strsplit(contrast[k], split= " - ")[[1]][1], split = "\\.")
    a <- a[[1]]
    b <- strsplit(strsplit(contrast[k], split= " - ")[[1]][2], split = "\\.")
    b <- b[[1]]
    
    c <- setdiff(a,b)
    d <- setdiff(b,a)
    
    length_c <- length(c)
    length_d <- length(d)
    
    if (length_c > 1) {
      
      for (i in 1:(length_c-1)) {
        for (j in (i+1):(length_c)) {
          
          e <- paste(c[i], c[j], sep = ".")
          
          if (length(grep(e, make.names(unique(meta[,2])))) == length(match(grep(c[i], make.names(unique(meta[,2]))), grep(c[j], contrast)))) {
            length_c = length_c - 1
          }
        }
      }
    }
    
    
    if (length_d > 1) {
      
      for (i in 1:(length(d)-1)) {
        for (j in (i+1):(length(d))) {
          
          f <- paste(d[i], d[j], sep = ".")
          
          if (length(grep(f, make.names(unique(meta[,2])))) == length(match(grep(d[i], make.names(unique(meta[,2]))), grep(d[j], contrast)))) {
            length_d = length_d - 1
          }
          
          
        }
      }
    }
    
    if (length_c <= 1 & length_d <= 1) {
    }
    
    if (length_c > 1 | length_d > 1) {
      contrast1 <- contrast1[-which(contrast1 == contrast[k])]
    }
    
  }
  
  return(contrast1)
  
}

################################################################################

#diff_expr

################################################################################


diff_expr <- function(data.expr, meta, comparisons) {
  
  
  #get grouping variable in correct order
  
  ph = as.data.frame(colnames(data.expr))
  colnames(ph) <- "cel_names"
  
  ph1 <- ph %>% inner_join(meta, by = c("cel_names" = "names"), 
                           match_fun = str_detect)
  
  if (is.null(ph1$Pairing)){
    
    #model design
    groups <- ph1$Grouping
    levels <- unique(ph1$Grouping)
    
    f <- factor(groups,levels=levels)
    
    design <- model.matrix(~ 0 + f)
    colnames(design) <- make.names(levels)
    
    
    #fit linear model
    
    data.fit <- lmFit(data.expr,design)
    
    
    #get contrasts
    contrast <- comparisons
  }
    
  
  if (!is.null(ph1$Pairing)){
    
    #Model design
    fg <- factor(ph1$Grouping)
    fp <- factor(ph1$Pairing)
    
    length(levels(fg))
    
    design = model.matrix(~ 0 + fg + fp)
    
    levels1 <- levels(fg)
    
    remainder <- paste0("a",seq(1,(ncol(design) - length(comparisons) - 1)))
    colnames(design) <- make.names(c(levels1, remainder))
    
    
    #Fit linear model
    data.fit = lmFit(data.expr, design)
    
    #get contrasts
    contrast <- comparisons
    
  }
  

    #make top table for each comparison
    top.table <- list()
    
    for (i in contrast) {
      contrast.matrix <- makeContrasts(contrasts = i,levels=design)
      
      data.fit.con <- contrasts.fit(data.fit,contrast.matrix)
      
      data.fit.eb <- eBayes(data.fit.con)
      
      top.table[[i]] <- topTable(data.fit.eb, sort.by = "P", n = Inf)
      
      top.table[[i]] <- cbind(rownames(top.table[[i]]), top.table[[i]])
      
      rownames(top.table[[1]]) <- NULL
      
      colnames(top.table[[i]]) <- c("Probeset.ID", colnames(top.table[[i]])[-1])
      
    }
  
  return(top.table)
}



################################################################################

#plotvolcano

################################################################################

plotvolcano <- function(top.table, p = "raw", p.threshold = 0.05, logFC.threshold = 1){
  
  if (p == "raw"){
    top.table$colour[(top.table$logFC < logFC.threshold & top.table$logFC > (-1 * logFC.threshold)) | top.table$P.Value > p.threshold] = "darkgrey"
    top.table$colour[top.table$logFC < (-1 * logFC.threshold) & top.table$P.Value <= p.threshold] = "blue"
    top.table$colour[top.table$logFC >= logFC.threshold & top.table$P.Value <= p.threshold] = "red"
    
    
    volcano <- ggplot(top.table, aes(x = logFC, y = -log10(P.Value), text = paste0("Probeset ID: ", Probeset.ID))) +
      geom_point(colour = top.table$colour) +
      geom_hline(yintercept = -log10(p.threshold), color = "grey", linetype = "dotted", size = 0.5) +
      geom_vline(xintercept = c(-1 *logFC.threshold, logFC.threshold), color = "grey", linetype = "dotted", size = 0.5) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  
  if (p == "adj"){
    top.table$colour[(top.table$logFC < logFC.threshold & top.table$logFC > (-1 * logFC.threshold)) | top.table$adj.P.Val > p.threshold] = "darkgrey"
    top.table$colour[top.table$logFC < (-1 * logFC.threshold) & top.table$adj.P.Val <= p.threshold] = "blue"
    top.table$colour[top.table$logFC >= logFC.threshold & top.table$adj.P.Val <= p.threshold] = "red"
    
    
    volcano <- ggplot(top.table, aes(x = logFC, y = -log10(adj.P.Val), text = paste0("Probeset ID: ", Probeset.ID))) +
      geom_point(colour = top.table$colour) +
      geom_hline(yintercept = -log10(p.threshold), color = "grey", linetype = "dotted", size = 0.5) +
      geom_vline(xintercept = c(-1 *logFC.threshold, logFC.threshold), color = "grey", linetype = "dotted", size = 0.5) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  
  
  
  return(ggplotly(volcano))
  
  
}

################################################################################

#get_chiptype

################################################################################

get_chiptype = function(annotation){
  attributes <- c("affy_hc_g110", "affy_hg_focus", "affy_hg_u133a", "affy_hg_u133a_2", 
                  "affy_hg_u133b", "affy_hg_u133_plus_2", "affy_hg_u95a", "affy_hg_u95av2",
                  "affy_hg_u95b", "affy_hg_u95c", "affy_hg_u95d", "affy_hg_u95d", "affy_hg_u95e",
                  "affy_hta_2_0", "affy_huex_1_0_st_v2", "affy_hugenefl", "affy_hugene_1_0_st_v1",
                  "affy_hugene_2_0_st_v1", "affy_hugene_2_1_st_v1", "affy_primeview", "affy_u133_x3p")
  
  
  attributes1 = NULL
  attributes2 = NULL
  
  for (i in 1:length(attributes)){
    attributes1[i] <- str_remove(attributes[i], "affy_")
    attributes2[i] <- stringdist::stringdist(annotation, attributes1[i])
  }
  
  chiptype <- attributes[attributes2 == min(attributes2)]
  
  if (length(chiptype) > 1){
    chiptype <- NULL
  }
  
  return(chiptype)
}

