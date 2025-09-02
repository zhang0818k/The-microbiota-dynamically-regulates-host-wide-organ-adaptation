# ==============================================================================
# SCRIPT: Differential Gene Expression (DEG) Analysis Pipeline
# DATE: September 2, 2025
#
# DESCRIPTION:
# This script performs differential gene expression analysis across multiple
# tissues and experimental conditions (GF, SPF, nABX, ExGF).
#
# WORKFLOW:
#   1. Loads processed and annotated Seurat objects for multiple tissues.
#   2. For each tissue, it filters for cell types with sufficient representation
#      across all experimental conditions.
#   3. It then iterates through predefined comparisons (e.g., GF vs. SPF) for
#      each filtered cell type.
#   4. The script runs Seurat's FindMarkers, filters for significant DEGs,
#      and saves the consolidated results into a CSV file for each tissue.
#
# ==============================================================================


#
# Section 0: Setup and Environment
# ------------------------------------------------------------------------------

# --- 0.1: Load Required Libraries
library(Seurat)
library(dplyr)
library(purrr)

# --- 0.2: Define Parameters and Paths
# Base directory where tissue-specific subfolders are located
BASE_DIR <- "/home/GM_SPF/"

# Directory to save the final DEG results
RESULTS_DIR <- "/home/GM_SPF/DEG"
# Create the results directory if it doesn't exist
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}

# List of tissues to be analyzed
TISSUE_NAMES <- c("aorta", "brain", "cecum", "mesentery", "colon", "duodenum", "heart", "ileum", "jejunum", "liver", "lung", "muscle", "pancreas", "PBMC", "spleen", "kidney", "stomach")

# Minimum number of cells required in each condition for a cell type to be included in DEG analysis
# (e.g., >2 in the original script is equivalent to >= 3)
MIN_CELLS_PER_CONDITION <- 3

# Define the experimental comparisons to be performed
# 'compare': Name of the comparison
# 'ident1'/'ident2': Suffixes for the conditions being compared
COMPARISON_TABLE <- data.frame(
  compare = c('GFvsSPF', 'nABXvsSPF', 'ExGFvsGF'),
  ident1  = c('_GF', '_nABX', '_ExGF'),
  ident2  = c('_SPF', '_SPF', '_GF')
)


#
# Section 1: Data Loading
# ------------------------------------------------------------------------------

#' Reads a processed Seurat object for a given tissue.
#'
#' @param tissue_name The name of the tissue (e.g., "colon").
#' @return A Seurat object if the file is found, otherwise NULL.
read_seurat <- function(tissue_name) {
  tissue_dir <- file.path(BASE_DIR, tissue_name)
  
  # **MODIFIED**: Looks for files ending with '_integrated_processed.rds'
  seurat_file <- list.files(
    path = tissue_dir,
    pattern = "_integrated_processed\\.rds$",
    full.names = TRUE
  )[1]
  
  if (is.na(seurat_file) || length(seurat_file) == 0) {
    warning(paste("No '*_integrated_processed.rds' file found in", tissue_dir))
    return(NULL)
  }
  
  cat("Loading data for:", tissue_name, "\n")
  seurat_obj <- readRDS(seurat_file)
  return(seurat_obj)
}


#
# Section 2: DEG Analysis Function
# ------------------------------------------------------------------------------

#' Performs DEG analysis for a single tissue across multiple comparisons and cell types.
#'
#' @param seurat_obj A processed and annotated Seurat object.
#' @param tissue_name The name of the tissue being analyzed.
#' @return This function writes a CSV file with DEG results and does not return an object.
perform_deg_analysis <- function(seurat_obj, tissue_name) {
  if (is.null(seurat_obj)) {
    warning(paste("Skipping DEG analysis for", tissue_name, "due to missing data."))
    return()
  }
  
  cat(paste("\n--- Starting DEG analysis for:", tissue_name, "---\n"))
  
  # --- 2.1: Prepare Seurat object for comparisons
  # Create a new metadata column combining cell annotation and experimental status
  seurat_obj$cell.type.condition <- paste(seurat_obj$annotation, seurat_obj$Status, sep = "_")
  Idents(seurat_obj) <- seurat_obj$cell.type.condition
  
  # --- 2.2: Filter cell types based on cell counts
  # Count cells for each annotation within each status
  cell_counts <- as.data.frame.matrix(table(seurat_obj$annotation, seurat_obj$Status))
  
  # Identify cell types that have at least MIN_CELLS_PER_CONDITION in all four conditions
  celltypes_to_analyze <- row.names(cell_counts)[
    cell_counts$ExGF >= MIN_CELLS_PER_CONDITION &
      cell_counts$GF   >= MIN_CELLS_PER_CONDITION &
      cell_counts$nABX >= MIN_CELLS_PER_CONDITION &
      cell_counts$SPF  >= MIN_CELLS_PER_CONDITION
  ]
  
  if (length(celltypes_to_analyze) == 0) {
    message(sprintf("No cell types in %s met the minimum cell count criteria for all conditions. Skipping.", tissue_name))
    return()
  }
  
  message(sprintf("Found %d cell types to analyze in %s.", length(celltypes_to_analyze), tissue_name))
  
  # --- 2.3: Iterate through comparisons and cell types to find DEGs
  all_deg_results <- list() # Store results in a list for efficiency
  
  for (i in 1:nrow(COMPARISON_TABLE)) {
    comparison_name <- COMPARISON_TABLE$compare[i]
    
    for (cell_type in celltypes_to_analyze) {
      ident1 <- paste0(cell_type, COMPARISON_TABLE$ident1[i])
      ident2 <- paste0(cell_type, COMPARISON_TABLE$ident2[i])
      
      # Check if both comparison groups exist in the data
      if (ident1 %in% levels(Idents(seurat_obj)) && ident2 %in% levels(Idents(seurat_obj))) {
        tryCatch({
          # Run FindMarkers
          deg_markers <- FindMarkers(
            seurat_obj,
            ident.1 = ident1,
            ident.2 = ident2,
            logfc.threshold = 0.25, # A lower threshold to capture more before filtering
            verbose = FALSE
          )
          
          # Filter for significant DEGs and add metadata
          if (nrow(deg_markers) > 0) {
            deg_markers <- deg_markers %>%
              filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>% # Apply logFC filter here
              filter(!startsWith(rownames(.), "mt-")) # More robust mitochondrial gene filtering
            
            if (nrow(deg_markers) > 0) {
              deg_markers$Gene_symbol <- rownames(deg_markers)
              deg_markers$Cell_type <- cell_type
              deg_markers$Change <- ifelse(deg_markers$avg_log2FC > 0, 'Up', 'Down')
              deg_markers$Tissue <- tissue_name
              deg_markers$Comparison <- comparison_name
              all_deg_results <- append(all_deg_results, list(deg_markers))
            }
          }
        }, error = function(e) {
          message(sprintf("Skipping comparison for %s: %s", ident1, e$message))
        })
      }
    }
  }
  
  # --- 2.4: Consolidate and save results
  if (length(all_deg_results) > 0) {
    final_deg_df <- bind_rows(all_deg_results)
    output_file <- file.path(RESULTS_DIR, paste0('DEG_results_', tissue_name, '.csv'))
    write.csv(final_deg_df, output_file, row.names = FALSE)
    message(sprintf("SUCCESS: Saved DEG results for %s. Found %d total DEGs.", tissue_name, nrow(final_deg_df)))
  } else {
    message(sprintf("No significant DEGs found for %s after filtering.", tissue_name))
  }
}


#
# Section 3: Main Execution
# ------------------------------------------------------------------------------
cat("\n--- Starting analysis for all tissues ---\n")

# Use purrr::map to read all Seurat objects into a named list
seurat_list <- map(TISSUE_NAMES, read_seurat)
names(seurat_list) <- TISSUE_NAMES

# Use purrr::walk2 to apply the DEG analysis function to each object and its name
walk2(seurat_list, names(seurat_list), perform_deg_analysis)

cat("\n--- All analyses complete. ---\n")
