# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Source global file
setwd("C:/Users/jarno/OneDrive/Documents/GitHub/ArrayAnalysis_Shiny")
source("global.R")

#download gmt via: https://data.wikipathways.org/current/gmt/

WPversion <- "20250610"
organisms <- c("Homo sapiens",
               "Bos taurus",
               "Caenorhabditis elegans",
               "Mus musculus",
               "Rattus norvegicus")

ensembl <- useMart("ensembl")
gmt_all <- list()
for (o in organisms){
  
  # Read GMT file
  gmt <- clusterProfiler::read.gmt.wp(paste0("PrepareFiles/wikipathways-", WPversion, "-gmt-",
                                             stringr::str_replace(o, " ", "_"),
                                             ".gmt"))
  
  # Select biomaRt dataset
  biomaRt_dataset <- switch(o,
                            "Homo sapiens" = "hsapiens_gene_ensembl" ,
                            "Bos taurus" = "btaurus_gene_ensembl",
                            "Caenorhabditis elegans" = "celegans_gene_ensembl",
                            "Mus musculus" = "mmusculus_gene_ensembl",
                            "Rattus norvegicus" = "rnorvegicus_gene_ensembl"
  )
  
  # Get annotations
  ensembl_dataset <- useDataset(biomaRt_dataset, mart=ensembl)
  
  if (o == "Homo_sapiens"){
    annotations <- getBM(attributes=c("ensembl_gene_id",
                                      "entrezgene_id",
                                      "hgnc_symbol"), 
                         filters = "entrezgene_id",
                         values = gmt$gene[gmt$species == o],
                         mart = ensembl_dataset)
  } else{
    annotations <- getBM(attributes=c("ensembl_gene_id",
                                      "entrezgene_id",
                                      "external_gene_name"), 
                         filters = "entrezgene_id",
                         values = gmt$gene[gmt$species == o],
                         mart = ensembl_dataset)
  }

  
  annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
  temp <- left_join(gmt[gmt$species == o, ], annotations, by = c("gene" = "entrezgene_id"))
  colnames(temp) <- c("name", "version", "wpid",
            "species", "ENTREZID", "ENSEMBL", "SYMBOL")
  
  gmt_all[[o]] <- temp
}




test <- listAttributes(ensembl_dataset)

# Save file
save(gmt_all, file = "Objects/gmt_WP_all.RData")
save(WPversion, file = "Objects/WPversion.RData")


sapply(c("Homo_sapiens","Bos_taurus"), function(x) {
  switch(x,
         "Homo_sapiens" = "hsapiens_gene_ensembl" ,
         "Bos_taurus" = "btaurus_gene_ensembl",
         "Caenorhabditis_elegans" = "celegans_gene_ensembl",
         "Mus_musculus" = "mmusculus_gene_ensembl",
         "Rattus_norvegicus" = "rnorvegicus_gene_ensembl"
  )
})
