library(data.table)
library(dplyr)
library(MAGMA.Celltyping)

#run based on R-4.3.0

rm(list = ls())
gc()

#### Cell type Enrichment ####
### note: file path do not be too long in this analysis, otherwise it will pose an error ###
CellTypeDataset_dir <- "E:/04.Specific_Cell_Types/ScRNA_seq_dataset/" #modify into your path
MAGMA_result_dir <- "E:/04.Specific_Cell_Types/MGAMA_results/" #modify into your path

phenos <- "Anxiety2021"
CellTypeDataset <- "ctd_BlueLake2018_FrontalCortexOnly"

### move MAGMA results to fixed name path (if not do so, it may pose an error that can not find MAGMA results) ###
for(i in 1:length(phenos))
{
    dir.create(paste0(MAGMA_result_dir,"MAGMA_Files/",phenos[i],".35UP.10DOWN/"),recursive = TRUE)
}

for(i in 1:length(phenos))
{
  file.rename(from = paste0(MAGMA_result_dir,phenos[i],".35UP.10DOWN.genes.out"),to = paste0(MAGMA_result_dir,"MAGMA_Files/",phenos[i],".35UP.10DOWN/",phenos[i],".35UP.10DOWN.genes.out"))
  file.rename(from = paste0(MAGMA_result_dir,phenos[i],".35UP.10DOWN.genes.raw"),to = paste0(MAGMA_result_dir,"MAGMA_Files/",phenos[i],".35UP.10DOWN/",phenos[i],".35UP.10DOWN.genes.raw"))
}

for(i in 1:length(phenos))
{
  for(j in 1:length(CellTypeDataset))
  {
    CellTypeDataset <- readRDS(paste0(CellTypeDataset_dir,CellTypeDataset[j],".rds"))
    MAGMA_result_dir_modified <- paste0(MAGMA_result_dir,"MAGMA_Files/",phenos[i],".35UP.10DOWN")
    #### Linear mode ####
    MAGMA.Celltyping::calculate_celltype_associations(
      ctd = CellTypeDataset,
      magma_dir = MAGMA_result_dir_modified,
      EnrichmentMode = "Linear",
      ctd_species = "hsapiens",
      force_new=T)

    #### top 10% mode ####
    MAGMA.Celltyping::calculate_celltype_associations(
      ctd = CellTypeDataset,
      magma_dir = MAGMA_result_dir_modified,
      EnrichmentMode = "Top 10%",
      ctd_species = "hsapiens",
      force_new=T)
  }
}
