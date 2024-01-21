# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Source global file
source("D:/ArrayAnalysis_Shiny/global.R")

# Get filter sets for each organism
sets <- c("hsapiens_gene_ensembl" ,
          "btaurus_gene_ensembl",
          "celegans_gene_ensembl",
          "mmusculus_gene_ensembl",
          "rnorvegicus_gene_ensembl")
ensembl <- useMart("ensembl")
filterList <- list()
for (i in 1:length(sets)){
  ensembl1 <- useDataset(sets[i], mart=ensembl)
  all_filters <- biomaRt::listFilters(ensembl1)
  filterList[[i]] <- c("ensembl_gene_id",
                       "entrezgene_id",
                       setdiff(all_filters$name[str_detect(all_filters$name,"affy")],
                               all_filters$name[str_detect(all_filters$name,"with")]))
}

names(filterList) <- sets

# Save file
save(filterList, file = "filterList_biomaRt.RData")
