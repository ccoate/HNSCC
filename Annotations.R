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

# ---------------------------- Cluster profiling -------------------------- ####
library(clusterProfiler)
library(org.Hs.eg.db)   # Use org.Mm.eg.db for mouse genes
library(enrichplot)     # For visualization
library(ggplot2)

gene_symbols <- markers$gene 
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
markers.df <- markers
markers.df <- merge(markers.df, gene_entrez, by.x = "gene", by.y = "SYMBOL")

gene_list <- markers.df$avg_log2FC  # Ranked fold changes
names(gene_list) <- markers.df$ENTREZID  # Assign Entrez IDs as names

marker_list <- split(markers.df, markers.df$cluster)
names(marker_list) <- paste0("cluster_", names(marker_list))

# Loop over clusters and assign each ranked gene list to gene_list_<cluster>
gene_lists <- list()
for (i in seq_along(marker_list)) {
  cluster_id <- names(marker_list)[i]
  df <- marker_list[[i]]
  
  # Filter for rows that have valid ENTREZIDs and logFCs
  df <- df[!is.na(df$ENTREZID) & !is.na(df$avg_log2FC), ]
  
  # Create ranked gene list
  gene_list <- df$avg_log2FC
  names(gene_list) <- df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Add to the gene_lists list
  gene_lists[[cluster_id]] <- gene_list
}

#GO
dotplot_list <- list()
for (i in seq_along(gene_lists)) {
  cluster_id <- names(gene_lists)[i]
  gene_vector <- gene_lists[[i]]
  gene_ids <- names(gene_vector)  # ENTREZ IDs
  
  go_results <- enrichGO(
    gene = gene_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # Only add to the list if enrichment returned results
  if (nrow(go_results@result) > 0) {
    dot <- dotplot(go_results, showCategory = 5) +
      ggtitle(paste0("Cluster ", cluster_id, " GO Enrichment")) +
      theme(plot.title = element_text(size = 10),
            axis.text.y = element_text(size = 8),
            )
    
    dotplot_list[[i]] <- dot
  } else {
    message(paste("No enrichment found for cluster", cluster_id))
  }
}

# Combine all dotplots together
go_combined_plot <- wrap_plots(dotplot_list, ncol = 2) +
  plot_annotation(title = "GO Enrichment Across Clusters")

# Print or save
print(go_combined_plot)
ggsave("GO_Enrichment_AllClusters.png", plot = go_combined_plot, width = 12, height = 10, dpi = 300)

vln_list <- lapply(GOIs, function(gene) {
  VlnPlot(sample_list[[i]], features = gene, pt.size = 0.2) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 4)) +
    ggtitle(gene)
})
vln_combined <- wrap_plots(vln_list, ncol = 4) +
  plot_annotation(title = paste("Sample", i, "| Violin Plot of Genes of Interest"))


go_results <- enrichGO(
  gene = names(gene_list),  
  OrgDb = org.Hs.eg.db,     
  keyType = "ENTREZID",      
  ont = "BP",                # "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component)
  pAdjustMethod = "BH",      
  pvalueCutoff = 0.05,       
  qvalueCutoff = 0.05        
)
dotplot(go_results, showCategory = 10)  # Top 20 terms

#KEGG
KEGGplot_list <- list()
for (i in seq_along(gene_lists)) {
  cluster_id <- names(gene_lists)[i]
  gene_vector <- gene_lists[[i]]
  gene_ids <- names(gene_vector)  # ENTREZ IDs
  
  KEGG_results <- enrichKEGG(
    gene = gene_ids,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # Only add to the list if enrichment returned results
  if (nrow(KEGG_results@result) > 0) {
    KEGG <- dotplot(KEGG_results, showCategory = 5) +
      ggtitle(paste0("Cluster ", cluster_id, " KEGG Enrichment")) +
      theme(plot.title = element_text(size = 10),
            axis.text.y = element_text(size = 8),
      )
    
    KEGGplot_list[[i]] <- KEGG
  } else {
    message(paste("No enrichment found for cluster", cluster_id))
  }
}

# Combine all dotplots together
KEGG_combined_plot <- wrap_plots(KEGGplot_list, ncol = 2) +
  plot_annotation(title = "KEGG Enrichment Across Clusters")

# Print or save
print(KEGG_combined_plot)
ggsave("KEGG_Enrichment_AllClusters.png", plot = KEGG_combined_plot, width = 12, height = 10, dpi = 300)

