library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# --------------------------- General References -------------------------- ####
#setup reference
ref <- celldex::HumanPrimaryCellAtlasData()

#create single cell experiment from seurat
scrdataTSCC <- as.SingleCellExperiment(merged)
#run singleR
scroutputTSCC <- SingleR(test = scrdataTSCC, ref = ref, assay.type.test = 1, labels = ref$label.main)
#port back
merged[["SingleR.labels"]] <- scroutputTSCC$labels
#check it
DimPlot(merged, group.by = c("seurat_clusters","SingleR.labels"), label = T, repel = T)

# --------------------------- Specific References -------------------------- ####
# ------------- Merged Reference Builder
setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/OSCC/")

base_path <- "/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/OSCC/"

sample_names <- paste0("sample_", 1:9)

load_ref_sample <- function(sample_name) {
  sample_path <- paste0(base_path, sample_name)
  tissue <- Read10X_h5(paste0(sample_path, "/filtered_feature_bc_matrix.h5"))
  tissue <- CreateSeuratObject(counts = tissue)
  return(tissue)
}

load_celltype <- function(sample_name) {
  celltype_path <- paste0(base_path, sample_name, "/celltype_metadata.csv")
  cell_ids <- read_csv(celltype_path)
}

ref_list <- lapply(seq_along(sample_names), function(i) {
  obj <- load_ref_sample(sample_names[i])
  cell_ids <- load_celltype(sample_names[i])
  annotated_cells <- cell_ids$cell
  obj <- subset(obj, cells = annotated_cells)
  obj$celltype <- cell_ids$celltype_classification
  return(obj)
})

ref_merge <- merge(ref_list[[1]], ref_list[-1], add.cell.ids = sample_names)
ref_merge <- JoinLayers(ref_merge)
ref_merge <- NormalizeData(ref_merge)
ref_merge <- as.SingleCellExperiment(ref_merge)


scrdataTSCC <- as.SingleCellExperiment(merged)
scroutputTSCC <- SingleR(test = scrdataTSCC, ref = ref_merge, assay.type.test = 1, labels = ref_merge$celltype)

merged[["SingleR.labels"]] <- scroutputTSCC$labels
DimPlot(merged, group.by = c("seurat_clusters","SingleR.labels"), label = T, repel = T)

saveRDS(ref_merge, file = file.path("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/ref_merge_OSCC"))

# ------------- Single Reference
setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/OSCC/")

# Build object to use as reference
ref_obj <- Read10X_h5("./filtered_feature_bc_matrix.h5")
ref_obj <- CreateSeuratObject(counts = ref_obj)

# Add annotations and subset for only annotated cells
cell_ids <- read_csv("./celltype_metadata-P1.csv")
annotated_cells <- cell_ids$cell
ref_obj <- subset(ref_obj, cells = annotated_cells)
ref_obj$celltype <- cell_ids$celltype_classification
ref_obj <- NormalizeData(ref_obj)

# Transform reference to single cell experiment
ref_OSCC <- as.SingleCellExperiment(ref_obj)

# Create single cell experiment from query Seurat
scrdataTSCC <- as.SingleCellExperiment(merged)
scroutputTSCC <- SingleR(test = scrdataTSCC, ref = ref_OSCC, assay.type.test = 1, labels = ref_OSCC$celltype)
#port back
merged[["SingleR.labels"]] <- scroutputTSCC$labels
#check it
DimPlot(merged, group.by = c("seurat_clusters","SingleR.labels"), label = T, repel = T)

# ------------- OPSCC
setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/OPSCC")
cell_ids <- read_csv("./celltype_metadata-P1.csv")
annotated_cells <- cell_ids$cell
ref_obj <- subset(ref_obj, cells = annotated_cells)
ref_obj$celltype <- cell_ids$celltype_classification
ref_obj <- NormalizeData(ref_obj)

# ------------- HNSCC Ref
setwd("/Users/chandlercoate/Desktop/Ferrara_Lab/cSCC_Project/annotations/HNSCC/GSE181919")
HNSCC_counts <- read.delim("GSE181919_UMI_counts.txt.gz", check.names = FALSE)
HNSCC_metadata <- read.delim("GSE181919_Barcode_metadata.txt.gz")
rownames(HNSCC_metadata) <- colnames(HNSCC_counts)
HNSCC_metadata <- subset(HNSCC_metadata, !tissue.type %in% c("LP", "LN"))
common_cells <- intersect(colnames(HNSCC_counts), rownames(HNSCC_metadata))
HNSCC_counts <- HNSCC_counts[, common_cells]
