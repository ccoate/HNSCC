library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(patchwork)
library(BiocManager)
library(cluster)

# rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)

setwd("/labs/kwferrar/chandler/head_neck")

# -------------------Functions------------------- ####
#Silhouette Function Josie
optimize_clustering_silhouette <- function(seurat_obj,
                                           range_pc = c(10, seq(20, 70, 10)),
                                           range_kn = c(10, seq(20, 40, 10)),
                                           range_res = 0.4) {
  message("Running PCA with max PCs: ", max(range_pc))
  seurat_obj <- RunPCA(seurat_obj, npcs = max(range_pc), verbose = FALSE)
  best_score <- -Inf
  best_params <- list(pc = NA, k = NA, res = NA)
  results <- data.frame(pc = numeric(), k = numeric(), res = numeric(),
                        silhouette = numeric(), n_clusters = numeric())
  for (p in range_pc) {
    for (k in range_kn) {
      for (r in range_res) {
        message(sprintf("Testing: PCs = %d | k = %d | res = %.2f", p, k, r))
        obj <- FindNeighbors(seurat_obj, dims = 1:p, k.param = k, verbose = FALSE)
        obj <- FindClusters(obj, resolution = r, algorithm = 4, verbose = FALSE)
        obj <- RunUMAP(obj, dims = 1:p, verbose = FALSE)
        dist.matrix <- dist(Embeddings(obj, reduction = "umap"))
        sil <- silhouette(as.numeric(as.factor(obj$seurat_clusters)), dist.matrix)
        avg_sil <- mean(summary(sil)$clus.avg.widths)
        n_clust <- length(unique(obj$seurat_clusters))
        results <- rbind(results, data.frame(pc = p, k = k, res = r,
                                             silhouette = avg_sil,
                                             n_clusters = n_clust))
        if (avg_sil > best_score) {
          best_score <- avg_sil
          best_params <- list(pc = p, k = k, res = r, silhouette = avg_sil)
        }
      }
    }
  }
  # Apply best parameters to the original object
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:best_params$pc, k.param = best_params$k, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = best_params$res, algorithm = 4, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:best_params$pc, verbose = FALSE)
  return(list(
    best_params = best_params,
    results = results,
    seurat_object = seurat_obj
  ))
}

# Correlation functions
tarcorr <- function(corset){
  #do the correlation
  cormat <- round(cor(as.matrix(corset),method = "pearson"), 2)
  #melt it
  melted_cormat <- reshape2::melt(cormat)
  upper_tri <- get_upper_tri(cormat)
  #melt again
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position.inside = c(0.6, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(position = "right")
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
#reordering
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

get_metadata_by_datasets <- function(dataset_names, metadata) {
  subset(metadata, Dataset %in% dataset_names)
}

# ----------Load in cGESPs_prior and GOI---------- ####
source("./scripts/cGESPs_5.12.25.R")

# ---------------------Load RDS------------------- ####
sample_list <- readRDS("./HNSCC_merged/seurat_list")

for (i in seq_along(sample_list)) {
  sample_list[[i]] <- SCTransform(sample_list[[i]], assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)
}

# -------------------Integration------------------ ####
features <- SelectIntegrationFeatures(sample_list, nfeatures = 3000, verbose = TRUE) #select a good feature number, 2000 or 3000

sample_list <- PrepSCTIntegration(sample_list, anchor.features = features, verbose = TRUE)

anchors <- FindIntegrationAnchors(sample_list, normalization.method = "SCT",
                                  anchor.features = features,
                                  verbose = TRUE)
merged <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                        verbose = TRUE)

set.seed(1234)

merged <- RunPCA(merged, npcs = 100, assay = "SCT", features = VariableFeatures(merged))

opt <- optimize_clustering_silhouette(merged)
merged <- opt$seurat_object
Idents(merged) <- "seurat_clusters"
best_params <- opt$best_params

print(best_params)

output_dir <- "./HNSCC_merged/outs/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save RDS after Clustering
saveRDS(merged, file = file.path(output_dir, "HNSCC_Integrated_RDS"))

# FindMarkers
merged <- PrepSCTFindMarkers(merged, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(merged, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)

# Save RDS after FindAllMarkers
write.csv(markers, file = file.path(output_dir, "HNSCC_AllMarkers.csv"), row.names = FALSE)

#calculating silhouette scores of above clustering
dist.matrix <- dist(Embeddings(merged,reduction = "umap"))
sil <- silhouette(as.numeric(as.factor(merged$seurat_clusters)), dist.matrix)
avg_sil <- mean(summary(sil)$clus.avg.widths)

#Sanity check
print(paste0("Dims: ", dim(merged)))
print(paste0("Silhouette: ", avg_sil))

# Plots ####

# ──────── UMAP Plot ────────
um <- DimPlot(merged, label = TRUE) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  ggtitle(paste0("Merged Clustering: ", 
                 length(sample_list_1_4), 
                 " Samples", 
                 "\n",
                 "| PCs: ", best_params$pc,
                 "| k: ", best_params$k,
                 "| res: ", best_params$res,
                 "| silhouette: ", round(avg_sil, 2))
)
# print(um)
# Save UMAP plot
umap_file <- file.path(output_dir, paste0("UMAP_Merged.png"))
ggsave(umap_file, plot = um, width = 8, height = 6, dpi = 300, bg = "transparent")
rm(um)
gc()

# ──────── Heatmaps ────────
cGESPs_exp <- markers[markers$gene %in% cGESPs_prior & !(duplicated(markers$gene)), ]

hm_cGESP <- DoHeatmap(
  merged,
  features = cGESPs_exp$gene,
  disp.min = -2,
  disp.max = 2,
  slot = "scale.data"
) +
  ggtitle(paste0("Merged Optimized Clustering: ", length(sample_list_1_4), " Samples")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  guides(color = "none")
# Save heatmap
heatmap_file <- file.path(output_dir, "Heatmap_cGESP_Merged.png")
ggsave(heatmap_file, plot = hm_cGESP, width = 8, height = 6, dpi = 300, bg = "transparent")
rm(hm_cGESP)
gc()

top_markers <- subset(markers, markers$pct.1 >= 0.8)
top_markers <- subset(top_markers, top_markers$avg_log2FC > 1.5)
cGESP_80 <- subset(top_markers, gene %in% cGESPs_prior)

hm_cGESP80 <- DoHeatmap(merged, 
              features = cGESP_80$gene,
              disp.min = -2,
              disp.max = 2,
              slot = "scale.data",
              draw.lines = FALSE,
              angle = 0,
              size = 3
) +
  ggtitle(paste0("Expression of cGESPs")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  guides(color = "none")

heatmap_file_cGESP80 <- file.path(output_dir, "Heatmap_cGESP_80.png")
ggsave(heatmap_file_cGESP80, plot = hm_cGESP80, width = 8, height = 6, dpi = 300, bg = "transparent")
rm(hm_cGESP80)

comparisons <- c("seurat_clusters", "orig.ident")
for (i in seq_along(comparisons)) {
  hm <- DoHeatmap(
    merged,
    features = GOI,
    disp.min = -2,
    disp.max = 2,
    slot = "scale.data",
    group.by = comparisons[[i]],
    draw.lines = FALSE
  ) +
    ggtitle(paste0("Expression by ", comparisons[[i]])) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    guides(color = "none")
  heatmap_file_sample <- file.path(output_dir, paste0("Heatmap_", comparisons[[i]], ".png"))
  ggsave(heatmap_file_sample, plot = hm, width = 8, height = 6, dpi = 300, bg = "transparent")
  rm(hm) 
}

# ────── Violin ──────
expr_values <- FetchData(object = merged, vars = GOI)
y_max_val <- max(expr_values, na.rm = TRUE)

for (i in seq_along(comparisons)) {
  vln_list <- lapply(GOI, function(gene) {
    VlnPlot(merged, group.by = comparisons[[i]], features = gene, pt.size = FALSE) +
      theme(legend.position = "none",
            plot.title = element_text(size = 12),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10)) +
      scale_y_continuous(limits = c(0, y_max_val)) +
      ggtitle(gene, )
  })
  vln_combined <- wrap_plots(vln_list, ncol = 4) +
    plot_annotation(title = paste0("Violin Plot: ", comparisons[[i]]))
  # print(vln_combined)
  # Save violin plot
  violin_file <- file.path(output_dir, paste0("Violin_", comparisons[[i]], ".png"))
  ggsave(violin_file, plot = vln_combined, width = 13, height = 7.5, dpi = 300, bg = "transparent")  # Adjust size for grid
  rm(vln_combined)
  gc()
}

# ────── Correlation ──────
corset <- FetchData(merged, vars = GOI)
cp <- tarcorr(corset)
corrplot_file <- file.path(output_dir, paste0("Corrplot_Merged.png"))
ggsave(corrplot_file, plot = cp, width = 8, height = 6, dpi = 300, bg = "transparent")
rm(corrplot_file, cp)

# ────── Spatial Dim ──────
slices <- merged@images
SD_list <- lapply(names(slices), function(slice_name) {
  SpatialDimPlot(merged, 
                 group.by = "seurat_clusters", 
                 images = slice_name,
                 pt.size.factor = 2) +
    ggtitle(slice_name)
})
SD_combined <- wrap_plots(SD_list, ncol = 4)
# Save the plot
SD_file <- file.path(output_dir, "Spatial_Dim_Merged.png")
ggsave(SD_file, plot = SD_combined, width = 13, height = 7.5, dpi = 300, bg = "transparent")
# Clean up
rm(SD_combined)
gc()

