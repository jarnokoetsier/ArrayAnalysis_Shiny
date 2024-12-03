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

# Get required CRAN packages
CRANpackages <- read.table("Objects/CRANpackages.txt", header = TRUE)

# Install remotes package
if (!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes", ask = FALSE)
}
require("remotes", character.only = TRUE)

#Install (if not yet installed) and load the required packages: 
for (pkg in 1:nrow(CRANpackages)) {
  
  # Install package if not installed
  if (!requireNamespace(CRANpackages$name[pkg], quietly = TRUE)){
    remotes::install_version(CRANpackages$name[pkg], 
                              version = CRANpackages$version[pkg], 
                              repos = "http://cran.us.r-project.org")
  } else {
    
    #Install package if package version is incorrect:
    
    # Install if package version is smaller than needed
    if (CRANpackages$direction[pkg] == "smaller"){
      if (packageVersion(CRANpackages$name[pkg]) < CRANpackages$version[pkg]){
        remotes::install_version(CRANpackages$name[pkg],
                                 CRANpackages$version[pkg])
      }
    }
    # Install if package version is different than needed
    if (CRANpackages$direction[pkg] == "unequal"){
      if (packageVersion(CRANpackages$name[pkg]) != CRANpackages$version[pkg]){
        remotes::install_version(CRANpackages$name[pkg],
                                 CRANpackages$version[pkg])
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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (pkg in 1:nrow(BiocPackages)) {
  if (!requireNamespace(BiocPackages$name[pkg], quietly = TRUE)){
    BiocManager::install(BiocPackages$name[pkg], ask = FALSE, version = "3.20")
  } else{
    if (packageVersion(BiocPackages$name[pkg]) < BiocPackages$version[pkg]){
      BiocManager::install(BiocPackages$name[pkg], ask = FALSE, version = "3.20")
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

