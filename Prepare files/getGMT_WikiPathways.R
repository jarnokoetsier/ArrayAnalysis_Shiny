# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Source global file
setwd("C:/Users/jarno/OneDrive/Documents/GitHub/ArrayAnalysis_Shiny")
source("global.R")

organisms <- c("Homo_sapiens",
               "Bos_taurus",
               "Caenorhabditis_elegans",
               "Mus_musculus",
               "Rattus_norvegicus")

ensembl <- useMart("ensembl")
gmt_all <- list()
for (o in organisms){
  
  # Read GMT file
  gmt <- clusterProfiler::read.gmt.wp(paste0("D:/ArrayAnalysis_Shiny/Prepare files/wikipathways-20231210-gmt-",o,".gmt"))
  
  # Select biomaRt dataset
  biomaRt_dataset <- switch(o,
                            "Homo_sapiens" = "hsapiens_gene_ensembl" ,
                            "Bos_taurus" = "btaurus_gene_ensembl",
                            "Caenorhabditis_elegans" = "celegans_gene_ensembl",
                            "Mus_musculus" = "mmusculus_gene_ensembl",
                            "Rattus_norvegicus" = "rnorvegicus_gene_ensembl"
  )
  
  # Get annotations
  ensembl_dataset <- useDataset(biomaRt_dataset, mart=ensembl)
  
  if (o == "Homo_sapiens"){
    annotations <- getBM(attributes=c("ensembl_gene_id",
                                      "entrezgene_id",
                                      "hgnc_symbol"), 
                         filters = "entrezgene_id",
                         values = gmt$gene,
                         mart = ensembl_dataset)
  } else{
    annotations <- getBM(attributes=c("ensembl_gene_id",
                                      "entrezgene_id",
                                      "external_gene_name"), 
                         filters = "entrezgene_id",
                         values = gmt$gene,
                         mart = ensembl_dataset)
  }

  
  annotations$entrezgene_id <- as.character(annotations$entrezgene_id)
  temp <- left_join(gmt, annotations, by = c("gene" = "entrezgene_id"))
  colnames(temp) <- c("name", "version", "wpid",
            "species", "ENTREZID", "ENSEMBL", "SYMBOL")
  
  gmt_all[[o]] <- temp
}




test <- listAttributes(ensembl_dataset)

# Save file
save(gmt_all, file = "Objectsgmt_WP_all.RData")


sapply(c("Homo_sapiens","Bos_taurus"), function(x) {
  switch(x,
         "Homo_sapiens" = "hsapiens_gene_ensembl" ,
         "Bos_taurus" = "btaurus_gene_ensembl",
         "Caenorhabditis_elegans" = "celegans_gene_ensembl",
         "Mus_musculus" = "mmusculus_gene_ensembl",
         "Rattus_norvegicus" = "rnorvegicus_gene_ensembl"
  )
})
