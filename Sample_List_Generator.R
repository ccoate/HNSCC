# H5 Write Function ####

setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/OSCC/")

filter_matrix <- Read10X("./filtered_feature_bc_matrix/")

write10xCounts("./filtered_feature_bc_matrix.h5", filter_matrix, type = "HDF5",
               genome = "grch38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))

# ---------------- RDS List generator - Chandler Coate ---------------- ####

setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/")
metadata <- read.csv("./Sample_Metadata.csv", header = TRUE)
dataset_pull <- subset(metadata, Diagnosis == "TSCC")

sample_names <- dataset_pull$Meta_Name

base_path <- ("./merge_test/data//")
  
# Load in samples ####
load_spatial_sample <- function(sample_name) {
  sample_path <- paste0(base_path, sample_name)
  img <- Read10X_Image(paste0(sample_path, "/spatial"), image.name = "tissue_hires_image.png")
  tissue <- Load10X_Spatial(
    data.dir = sample_path,
    assay = "Spatial",
    image = img
  )
  # Set the lowres scale factor to match hires
  tissue@images$slice1@scale.factors$lowres = tissue@images$slice1@scale.factors$hires
  # Keep only spots with counts
  zero_counts <- colSums(tissue[["Spatial"]]$counts) == 0
  if (sum(zero_counts) > 0) {
    tissue <- tissue[, !zero_counts]
  }
  return(tissue)
}

seurat_list <- lapply(seq_along(sample_names), function(i) {
  obj <- load_spatial_sample(sample_names[i])
  obj$orig.ident <- as.character(i)
  if (exists("dataset_pull")) {
    if (nrow(dataset_pull) >= i) {
      obj$tissue_origin <- dataset_pull$Diagnosis
      obj$HPV_status <- dataset_pull$HPVstatus
      obj$Sample_ID <- dataset_pull$Study_Sample_ID
      obj$Treatment_Status <- dataset_pull$Treatment_Status
      obj$Dataset <- dataset_pull$Dataset
    } else {
      warning(paste("No metadata associated for sample ", i))
    }
  }
  return(obj)
})

saveRDS(seurat_list, file = "./merge_test/TSCC/seurat_list")

