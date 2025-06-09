# -------------------------------- Packages ------------------------------- ####
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(magrittr)
library(CARD)
library(dplyr)
library(RColorBrewer)


# ---------------------------------- CARD --------------------------------- ####
Diagnosis <- "OSCC"

# Load Reference
if (Diagnosis == "HNSCC") {
  ref <- readRDS("./CC_annotations/SingleCellExperiment_Refs/ref_HNSCC.rds")
}   else if (Diagnosis == "OSCC") {
  ref <- readRDS("./CC_annotations/SingleCellExperiment_Refs/ref_OSCC.rds")
}   else if (Diagnosis == "OPSCC") {
  ref <- readRDS("./CC_annotations/SingleCellExperiment_Refs/ref_OPSCC.rds")
}

# Load Data
if (Diagnosis == "HNSCC_nanostring"){
  merged <- readRDS("./JF_RDS_files/final_objects/QL_nanostring_single_seed_1.rds")
}   else if (Diagnosis == "HNSCC") {
  merged <- readRDS("./Integrated_RDS/HNSCC/HNSCC_Integrated.rds")
}   else if (Diagnosis == "OSCC") {
  merged <- readRDS("./JF_RDS_files/final_objects/GSE208253_OSCC_GSE220978_OSCC_Zenodo_OSCC_integrate_seed_1.rds")
}   else if (Diagnosis == "OPSCC") {
  merged <- readRDS("./JF_RDS_files/final_objects/UQ698BB9E_OPSCC_Zenodo_OPSCC_integrate_seed_1.rds")
}

sample_list <- SplitObject(merged, split.by = "sample_id")
image_name <- unique(merged@meta.data$sample_id)

# CARD function
for (i in seq_along(sample_list)) {
  sample <- sample_list[[i]]
  SCT_count <- sample@assays$SCT@data
  
  sample@images$slice1 <- sample@images[[image_name[[i]]]]
  spatial_location <- as.data.frame(sample@images$slice1@boundaries$centroids@coords)
  rownames(spatial_location) <- colnames(SCT_count)
  spatial_location <- spatial_location[order(match(rownames(spatial_location), colnames(SCT_count))), , drop = FALSE]
  
  if (!identical(colnames(SCT_count), rownames(spatial_location))){
    #ensure that count matrix contains the same spots as the location matrix 
    SCT_count <- SCT_count[, colnames(SCT_count) %in% rownames(spatial_location)]
  }
  #match order
  spatial_location <- spatial_location[order(match(rownames(spatial_location), colnames(SCT_count))), , drop = FALSE]
  
  sc_count <- ref@assays@data$logcounts
  sc_meta <- ref@colData@listData
  print(all(rownames(sc_meta) == colnames(sc_count)))
  
  #deconvolute using CARD
  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = SCT_count,
    spatial_location = spatial_location,
    ct.varname = "celltype",
    ct.select = unique(sc_meta$celltype),
    sample.varname = NULL,
    minCountGene = 100,
    minCountSpot = 5) 
  
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  proportions <- as.data.frame(CARD_obj@Proportion_CARD)
  
  proportions$dominant_celltype <- apply(proportions, 1, function(x) names(proportions)[which.max(x)])
  
  sample <- Seurat::AddMetaData(sample,proportions)
  
  # Recreate original format
  sample@images[[image_name[[i]]]] <- sample@images$slice1
  
  sample@images$slice1 <- NULL
  
  sample_list[[i]] <- sample
}
