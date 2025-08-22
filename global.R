#==============================================================================#
# Name: global.R
# Description: global environment of the ArrayAnalysis Shiny app
#==============================================================================#

# ArrayAnalysis version:
ArrayAnalysis_version <- "ArrayAnalysis version 0.1.1"
online = FALSE

#******************************************************************************#
# Load functions
#******************************************************************************#

options(repos=c(CRAN="https://cran.r-project.org"))
source("FUNCTIONS.R")

#******************************************************************************#
#CRAN packages
#******************************************************************************#

# Get required CRAN packages
CRANpackages <- read.table("Objects/CRANpackages.txt", header = TRUE)

# Install remotes package
if (!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes", ask = FALSE)
}
require("remotes", character.only = TRUE)

#Install (if not yet installed) and load the required packages: 
for (pkg in 1:nrow(CRANpackages)) {
  
  # Install latest version
  if (CRANpackages$direction[pkg] == "smaller"){
    
    # Install package if not installed
    if (!requireNamespace(CRANpackages$name[pkg], quietly = TRUE)){
      install.packages(CRANpackages$name[pkg], 
                       ask = FALSE,
                       repos = "https://cloud.r-project.org")
    } else {
      
      #Install package if package version is too low:
      if (packageVersion(CRANpackages$name[pkg]) < CRANpackages$version[pkg]){
        if (CRANpackages$name[pkg] %in% .packages()){
          detach(paste0("package:", CRANpackages$name[pkg]), unload = TRUE,
                 character.only = TRUE)
        }
        install.packages(CRANpackages$name[pkg], 
                         ask = FALSE,
                         repos = "https://cloud.r-project.org",
                         upgrade = "never")
      }
    }
  }
  
  # Install specific version (not needed yet)
  if (CRANpackages$direction[pkg] == "unequal"){
    
    # Install package if not installed
    if (!requireNamespace(CRANpackages$name[pkg], quietly = TRUE)){
      remotes::install_version(CRANpackages$name[pkg],
                               CRANpackages$version[pkg],
                               repos = "https://cloud.r-project.org")
    } else {
      
      #Install package if package version is not correct
      if (packageVersion(CRANpackages$name[pkg]) != CRANpackages$version[pkg]){
        if (CRANpackages$name[pkg] %in% .packages()){
          detach(paste0("package:", CRANpackages$name[pkg]), unload = TRUE,
                 character.only = TRUE)
        }
        remotes::install_version(CRANpackages$name[pkg],
                                 CRANpackages$version[pkg],
                                 repos = "https://cloud.r-project.org",
                                 upgrade = "never")
      }
    }
  }
  
  require(as.character(CRANpackages$name[pkg]), character.only = TRUE)
}

#******************************************************************************#
#Bioconductor packages
#******************************************************************************#

# Get required Bioconductor packages:
BiocPackages <- read.table("Objects/Bioconductorpackages.txt", header = TRUE)

for (pkg in 1:nrow(BiocPackages)) {
  if (!requireNamespace(BiocPackages$name[pkg], quietly = TRUE)){
    BiocManager::install(BiocPackages$name[pkg], ask = FALSE)
  } else{
    if (packageVersion(BiocPackages$name[pkg]) < BiocPackages$version[pkg]){
      BiocManager::install(BiocPackages$name[pkg], ask = FALSE)
    }
  }
  require(as.character(BiocPackages$name[pkg]), character.only = TRUE)
}


#******************************************************************************#
# Example metadata
#******************************************************************************#

exampleMeta <- as.data.frame(data.table::fread("Data/Microarray/metaData_GSE6955.csv"))

#******************************************************************************#
# FilterList for biomaRt annotations
#******************************************************************************#

load("Objects/filterList_biomaRt.RData")

