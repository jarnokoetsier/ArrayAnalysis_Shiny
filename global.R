#==============================================================================#
# Name: global.R
# Description: global environment of the ArrayAnalysis Shiny app
#==============================================================================#

#******************************************************************************#
# Load functions
#******************************************************************************#

source("FUNCTIONS.R")

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
                  "data.table",
                  "DT",
                  "shiny",
                  "shinyWidgets",
                  "shinycssloaders",
                  "shinybusy",
                  "purrr",
                  "prompter",
                  "htmlwidgets",
                  "igraph",
                  "ggraph")

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
                  "bioDist",
                  "gcrma",
                  "plier",
                  "clusterProfiler",
                  "affy",
                  "oligo")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (pkg in BiocPackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}


#******************************************************************************#
# Example metadata
#******************************************************************************#

exampleMeta <- as.data.frame(fread("Data/metaData_GSE6955.csv"))

#******************************************************************************#
# FilterList for biomaRt annotations
#******************************************************************************#

load("Objects/filterList_biomaRt.RData")

