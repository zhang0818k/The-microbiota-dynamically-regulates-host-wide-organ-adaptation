# ==============================================================================
# SCRIPT: Pathway Enrichment Analysis Pipeline
# DATE: September 2, 2025
#
# DESCRIPTION:
# This script performs Gene Ontology (GO) enrichment analysis for Biological
# Processes (BP) on differentially expressed genes (DEGs). It processes DEGs
# from multiple tissues, cell types, and experimental comparisons. The workflow
# is designed to be efficient and reproducible by functionalizing the core
# analysis steps.
#
# WORKFLOW:
#   1. Loads the merged DEG results from the previous analysis step.
#   2. Defines a reusable function to perform GO enrichment on any given set of genes.
#   3. Splits the DEG data into groups based on the experimental comparison and
#      regulation direction (Up/Down).
#   4. Applies the enrichment function to each group of DEGs.
#   5. Combines all significant enrichment results into a single data frame.
#   6. Saves the final, filtered results to a CSV file.
#
# PREREQUISITES:
#   - Requires 'clusterProfiler' and 'org.Mm.eg.db' packages.
#   - Expects 'merged_DEG_all_tissues.csv' to be present in the specified directory.
#
# ==============================================================================


#
# Section 0: Setup and Environment
# ------------------------------------------------------------------------------

# --- 0.1: Load Required Libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(purrr)

# --- 0.2: Define Paths and Set Working Directory
# Directory where DEG results are stored and where pathway results will be saved.
RESULTS_DIR <- "/home/GM_SPF/Pathway"
DEG_FILE_PATH <- "/home/GM_SPF/DEG/merged_DEG_all_tissues.csv"

# Create the results directory if it doesn't exist and set as working directory
if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}
setwd(RESULTS_DIR)


#
# Section 1: Data Preparation
# ------------------------------------------------------------------------------
cat("--- 1. Preparing DEG data for analysis ---\n")

# --- 1.1: Load and filter the merged DEG data
deg_data <- read.csv(DEG_FILE_PATH) %>%
   # Create a unique identifier for each tissue/cell type combination
  mutate(tissue_cell_type = paste0(Tissue, "_", Cell_type))


#
# Section 2: Define Enrichment Analysis Function
# ------------------------------------------------------------------------------

#' Performs GO enrichment analysis on a subset of DEGs.
#'
#' @param deg_subset A data frame containing DEGs for a single group.
#'        Must contain columns: 'Gene_symbol', 'Tissue', 'tissue_cell_type', 'Cell_type'.
#' @return A data frame with significant enrichment results, or NULL if none are found.
perform_enrichment <- function(deg_subset) {
  
  # Ensure the input is not empty
  if (nrow(deg_subset) == 0) return(NULL)
  
  tryCatch({
    # --- Convert gene symbols to Entrez IDs
    gene_ids <- bitr(
      unique(deg_subset$Gene_symbol),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db
    )
    
    if (nrow(gene_ids) == 0) return(NULL)
    
    # --- Perform GO enrichment for Biological Process
    enrich_result <- enrichGO(
      gene          = gene_ids$ENTREZID,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05, # Use a less strict cutoff here, will filter on q-value later
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    
    # --- Format and return results if any are found
    if (!is.null(enrich_result) && nrow(enrich_result) > 0) {
      result_df <- as.data.frame(enrich_result) %>%
        mutate(
          Tissue = unique(deg_subset$Tissue),
          tissue_cell_type = unique(deg_subset$tissue_cell_type),
          Cell_type = unique(deg_subset$Cell_type)
        )
      return(result_df)
    }
    
    return(NULL)
    
  }, error = function(e) {
    # If an error occurs (e.g., no genes map), print a message and return NULL
    message(sprintf("Could not perform enrichment for group: %s. Error: %s", 
                    unique(deg_subset$tissue_cell_type), e$message))
    return(NULL)
  })
}


#
# Section 3: Run Analysis and Consolidate Results
# ------------------------------------------------------------------------------
cat("--- 2. Running enrichment analysis for all DEG groups ---\n")

# --- 3.1: Group data and apply the enrichment function to each group
# This splits the data by Comparison, Change, and tissue_cell_type,
# then runs the analysis on each small chunk of data.
all_enrichment_results <- deg_data %>%
  group_by(Comparison, Change, tissue_cell_type) %>%
  # nest() creates a list-column containing a small data frame for each group
  nest() %>%
  # pmap runs perform_enrichment on each row of the grouped data
  mutate(
    enrich_results = pmap(
      list(data, Comparison, Change), 
      ~ perform_enrichment(..1) %>%
        { if (!is.null(.)) mutate(., Comparison = ..2, Change = tolower(..3)) else . }
    )
  )

# --- 3.2: Combine all results into a single data frame
final_results <- all_enrichment_results %>%
  select(enrich_results) %>%
  unnest(cols = c(enrich_results))

# --- 3.3: Final filtering and save the results
if (nrow(final_results) > 0) {
  final_results_filtered <- final_results %>%
    filter(qvalue < 0.01) %>%
    # Reorder columns for clarity
    select(
      Comparison, Change, Tissue, Cell_type, tissue_cell_type,
      ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
      geneID, Count
    )
  
  write.csv(final_results_filtered, "Pathway_DEG_data_final.csv", row.names = FALSE)
  
  message(sprintf(
    "--- 3. Analysis complete. Found %d significant pathways. Results saved to 'Pathway_DEG_data_final.csv' ---", 
    nrow(final_results_filtered)
  ))
} else {
  message("--- 3. Analysis complete. No significant pathways found after final filtering. ---")
}
