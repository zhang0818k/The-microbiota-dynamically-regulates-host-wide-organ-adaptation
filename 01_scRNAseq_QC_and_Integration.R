# ==============================================================================
# SCRIPT: scRNA-seq Analysis Pipeline for Publication
# DATE: September 2, 2025
#
# DESCRIPTION:
# This script details the quality control, integration, and clustering of
# single-cell RNA sequencing data for the Mesentery tissue. The workflow
# includes:
#   1. Initial assessment of raw, un-integrated data for batch effects.
#   2. Per-sample QC: empty droplet removal, low-quality cell filtering,
#      and doublet detection.
#   3. Data integration using Harmony.
#   4. Unsupervised clustering and visualization.
#   5. Identification of cluster marker genes.
#
# ==============================================================================


#
# Section 0: Setup and Environment
# ------------------------------------------------------------------------------

# --- 0.1: Load Required Libraries
# Ensure all necessary packages are installed before running.
# install.packages(c("Seurat", "ggplot2", "harmony", "DropletUtils", "DoubletFinder", "BiocParallel"))

library(Seurat)
library(ggplot2)
library(harmony)
library(DropletUtils)
library(DoubletFinder)
library(BiocParallel) # For emptyDrops
library(dplyr)       # For data manipulation


# --- 0.2: Define Parameters and Set Working Directory
set.seed(123) # for reproducibility
options(future.globals.maxSize = 8000 * 1024^2) # Increase memory limit if needed

# --- Define paths and tissue name
TISSUE_NAME <- "Mesentery"
BASE_DIR <- "/home/zk_amas/GM_SPF/processed_data"
TISSUE_DIR <- file.path(BASE_DIR, TISSUE_NAME)

# --- Set working directory
setwd(TISSUE_DIR)

# --- Define custom color palette (if available)
# my36colors <- ... # Define your color vector here. For this example, we'll let Seurat choose.


#
# Section 1: Data Loading and Initial QC (Pre-Integration)
# ------------------------------------------------------------------------------
# In this section, we load the raw data for each sample and merge them without
# integration to visualize potential batch effects.

# --- 1.1: Load list of Seurat objects
# This file should contain a named list of Seurat objects, one for each sample.
seurat_list <- readRDS("~/GM_SPF/processed_data/mesentery_seurat_objects.rds")

# --- 1.2: Merge raw data for initial assessment
raw.merged <- merge(x = seurat_list[[1]],
                    y = unlist(seurat_list[2:length(seurat_list)]),
                    add.cell.ids = names(seurat_list))

# --- 1.3: Standard analysis on un-integrated data
raw.merged <- NormalizeData(raw.merged)
raw.merged <- FindVariableFeatures(raw.merged, selection.method = "vst", nfeatures = 2000)
raw.merged <- ScaleData(raw.merged)
raw.merged <- RunPCA(raw.merged, npcs = 30)
raw.merged <- RunUMAP(raw.merged, dims = 1:30)

# --- 1.4: Visualize batch effects
# UMAP plot split by sample (`orig.ident`) shows if cells cluster by sample
# of origin, indicating a strong batch effect that needs correction.
pdf("01_Raw_UMAP_BatchEffect_Check.pdf", width = 16, height = 6)
DimPlot(raw.merged, reduction = "umap", group.by = "orig.ident", raster = FALSE)
dev.off()

# Violin plots of QC metrics for each sample
pdf("02_Raw_VlnPlot_QC_Metrics.pdf", width = 12, height = 6)
VlnPlot(raw.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "orig.ident", pt.size = 0)
dev.off()

# Clean up large object before proceeding
rm(raw.merged)
gc()


#
# Section 2: Per-Sample Filtering and Quality Control
# ------------------------------------------------------------------------------
# We process each sample individually to remove empty droplets, filter cells
# based on QC metrics, and identify potential doublets.

# --- 2.1: Initialize a list to store cleaned Seurat objects
seurat_list_qc <- list()

# --- 2.2: Loop through each sample for QC
for (sample_name in names(seurat_list)) {
  
  cat(paste("\n--- Processing sample:", sample_name, "---\n"))
  
  seurat_obj <- seurat_list[[sample_name]]
  
  # Step A: Remove Empty Droplets using emptyDrops
  cat("Step A: Removing empty droplets...\n")
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  lower_bound <- sort(seurat_obj$nCount_RNA)[2] # Second least RNA count
  empty_droplets <- emptyDrops(counts, lower = lower_bound, BPPARAM = MulticoreParam())
  
  # Filter out cells identified as empty droplets
  # Note: some cells might be 'NA' in FDR, keep those (they are not empty)
  is.cell <- empty_droplets$FDR <= 0.01
  seurat_obj <- seurat_obj[, which(is.cell | is.na(is.cell))]
  
  # Step B: Filter Low-Quality Cells
  cat("Step B: Filtering low-quality cells...\n")
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > 300 &
                         nFeature_RNA < 6000 &
                         nCount_RNA < 30000 &
                         percent.mt < 15)
  
  # Step C: Detect and Annotate Doublets using DoubletFinder
  cat("Step C: Detecting doublets...\n")
  
  # Pre-process data for DoubletFinder
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # Needed for homotypic proportion
  
  # Find optimal pK value
  sweep.res <- paramSweep_v3(seurat_obj, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pK <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  
  # Estimate homotypic doublet proportion
  homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)
  
  # Calculate expected doublet rate (adjust this formula based on your platform)
  # General rule of thumb: ~0.8% per 1000 cells captured.
  # Formula used in original script: N_cells * 8 * 1e-6
  nExp_poi <- round((ncol(seurat_obj) * 8 * 1e-6) * ncol(seurat_obj))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(seurat_obj,
                                 PCs = 1:30,
                                 pN = 0.25,
                                 pK = optimal_pK,
                                 nExp = nExp_poi.adj,
                                 reuse.pANN = FALSE,
                                 sct = FALSE)
  
  # Standardize metadata column name for doublets
  df_colname <- colnames(seurat_obj@meta.data)[grep("^DF.classifications", colnames(seurat_obj@meta.data))]
  colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == df_colname] <- "DF.classification"
  
  # Store the cleaned object
  seurat_list_qc[[sample_name]] <- seurat_obj
}


#
# Section 3: Data Integration, Clustering, and Visualization
# ------------------------------------------------------------------------------
# Merge the cleaned samples, remove doublets, integrate using Harmony, and
# perform final clustering and UMAP visualization.

# --- 3.1: Merge the quality-controlled Seurat objects
all.merged <- merge(x = seurat_list_qc[[1]],
                    y = unlist(seurat_list_qc[2:length(seurat_list_qc)]),
                    add.cell.ids = names(seurat_list_qc))

# --- 3.2: Remove identified doublets
# In DoubletFinder, cells are classified as 'Singlet' or 'Doublet'.
all.merged <- subset(all.merged, subset = DF.classification == 'Singlet')

# --- 3.3: Perform data integration with Harmony
cat("\n--- Integrating data with Harmony ---\n")
all.merged <- NormalizeData(all.merged, verbose = FALSE)
all.merged <- FindVariableFeatures(all.merged, selection.method = "vst", nfeatures = 2000)
all.merged <- ScaleData(all.merged, verbose = FALSE)
all.merged <- RunPCA(all.merged, features = VariableFeatures(all.merged), npcs = 50, verbose = FALSE)

# Run Harmony to correct for batch effects from 'orig.ident'
all.merged <- RunHarmony(all.merged, group.by.vars = "orig.ident")

# --- 3.4: Post-integration clustering and UMAP
cat("--- Performing post-integration clustering and visualization ---\n")
all.merged <- RunUMAP(all.merged, reduction = "harmony", dims = 1:50)
all.merged <- FindNeighbors(all.merged, reduction = "harmony", dims = 1:50)
all.merged <- FindClusters(all.merged, resolution = 0.3) # Using a representative resolution

# --- 3.5: Save the final integrated Seurat object
saveRDS(all.merged, file = paste0(TISSUE_NAME, "_integrated_processed.rds"))


# --- 3.6: Visualize integrated results
pdf(paste0("03_", TISSUE_NAME, "_Integrated_UMAP_by_Cluster.pdf"), width = 8, height = 7)
DimPlot(all.merged, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf(paste0("04_", TISSUE_NAME, "_Integrated_UMAP_by_Sample.pdf"), width = 16, height = 6)
DimPlot(all.merged, reduction = "umap", group.by = "orig.ident")
dev.off()


#
# Section 4: Cluster Marker Gene Identification
# ------------------------------------------------------------------------------

cat("\n--- Finding markers for all clusters ---\n")

# Set the identity to the desired clustering resolution
DefaultAssay(all.merged) <- "RNA"
Idents(all.merged) <- all.merged$seurat_clusters # Or RNA_snn_res.0.3 if that's the column name

# Find all markers (positive markers only)
all.markers <- FindAllMarkers(all.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers to a file
write.csv(all.markers, file = paste0(TISSUE_NAME, "_cluster_markers_res0.3.csv"), row.names = FALSE)

cat("\n--- Analysis complete! ---\n")