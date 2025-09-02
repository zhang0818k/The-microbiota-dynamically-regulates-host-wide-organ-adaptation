# ==============================================================================
# SCRIPT: Recovery Score Calculation Pipeline
# AUTHOR: [Your Name]
# DATE: September 2, 2025
#
# DESCRIPTION:
# This script calculates "recovery scores" for both cell proportions and
# differentially expressed genes (DEGs) following different treatments.
# The goal is to quantify the extent to which changes observed in Germ-Free (GF)
# mice compared to SPF mice are reversed by antibiotic treatment (nABX) or
# co-housing (ExGF).
#
# WORKFLOW:
#   1.  **Cell Proportion Recovery**:
#       - Loads cell proportion change data.
#       - Calculates recovery scores for nABX and ExGF treatments.
#       - Classifies each cell type based on its recovery pattern.
#       - Saves detailed and summary results.
#   2.  **DEG Recovery**:
#       - Loads all tissue-specific DEG files.
#       - Calculates recovery scores at the individual gene level.
#       - Summarizes the scores for each cell type.
#       - Classifies each cell type based on the average recovery of its DEGs.
#       - Saves detailed and summary results.
#
# RECOVERY SCORE LOGIC:
#   - nABX_recovery = log2FC(nABXvsSPF) / log2FC(GFvsSPF)
#   - ExGF_recovery = log2FC(ExGFvsGF) / -log2FC(GFvsSPF)
#   A score > 0.5 is considered "reversible".
#
# ==============================================================================


#
# Section 0: Setup and Environment
# ------------------------------------------------------------------------------

# --- 0.1: Load Required Libraries
library(dplyr)
library(purrr)
library(readr)
library(tidyr)

# --- 0.2: Define Input/Output Paths
# Directory containing input files
PROP_DATA_PATH <- "/home/GM_SPF/prop_results/all_tissues_full_results.csv"
DEG_DATA_DIR <- "/home/GM_SPF/DEG/"

# Directory to save the output results
RESULTS_DIR <- "/home/GM_SPF/DEG/recovery_scores"
if (!dir.exists(RESULTS_GUDIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE)
}


#
# Section 1: Cell Proportion Recovery Analysis
# ------------------------------------------------------------------------------
cat("--- 1. Starting Cell Proportion Recovery Analysis ---\n")

# --- 1.1: Load and process cell proportion data
all_prop_results <- read.csv(PROP_DATA_PATH)

# --- 1.2: Calculate recovery scores and classify
recovery_data_prop <- all_prop_results %>%
  # Focus on the relevant comparisons
  filter(comparison %in% c("GFvsSPF", "nABXvsSPF", "ExGFvsGF")) %>%
  select(celltype, tissue, comparison, log2FC, FDR, is_significant) %>%
  # Pivot the data to have one row per cell type, with columns for each comparison
  pivot_wider(
    names_from = comparison,
    values_from = c(log2FC, FDR, is_significant)
  ) %>%
  # Calculate recovery scores
  mutate(
    min_change_threshold = 0.5, # Minimum |log2FC| in GFvsSPF to be considered
    
    # nABX recovery: Ratio of nABX effect to GF effect
    nABX_recovery_score = log2FC_nABXvsSPF / log2FC_GFvsSPF,
    
    # ExGF recovery: Ratio of ExGF reversal effect to GF effect
    ExGF_recovery_score = log2FC_ExGFvsGF / (-log2FC_GFvsSPF),
    
    # Set scores to NA if the initial change in GFvsSPF was not significant or too small
    nABX_recovery_score = ifelse(
      abs(log2FC_GFvsSPF) < min_change_threshold | !is_significant_GFvsSPF,
      NA, nABX_recovery_score
    ),
    ExGF_recovery_score = ifelse(
      abs(log2FC_GFvsSPF) < min_change_threshold | !is_significant_GFvsSPF,
      NA, ExGF_recovery_score
    ),
    
    # Classify recovery based on the scores
    recovery_class = case_when(
      abs(log2FC_GFvsSPF) < min_change_threshold ~ "No change",
      !is.na(nABX_recovery_score) & nABX_recovery_score > 0.5 & 
        !is.na(ExGF_recovery_score) & ExGF_recovery_score > 0.5 ~ "Both reversible",
      !is.na(nABX_recovery_score) & nABX_recovery_score > 0.5 ~ "nABX reversible",
      !is.na(ExGF_recovery_score) & ExGF_recovery_score > 0.5 ~ "ExGF reversible",
      TRUE ~ "Not reversible"
    ),
    tissue_cell = paste(tissue, celltype, sep = ": ")
  ) %>%
  # Keep only rows where scores could be calculated
  filter(!is.na(nABX_recovery_score) | !is.na(ExGF_recovery_score))

# --- 1.3: Calculate summary statistics for proportion recovery
recovery_stats_prop <- recovery_data_prop %>%
  group_by(recovery_class) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

# --- 1.4: Save results
write.csv(recovery_data_prop, file.path(RESULTS_DIR, "cell_proportion_recovery_details.csv"), row.names = FALSE)
write.csv(recovery_stats_prop, file.path(RESULTS_DIR, "cell_proportion_recovery_summary.csv"), row.names = FALSE)

cat("--- Cell Proportion Recovery Analysis Complete. Results saved. ---\n\n")


#
# Section 2: DEG Recovery Analysis
# ------------------------------------------------------------------------------
cat("--- 2. Starting DEG Recovery Analysis ---\n")

# --- 2.1: Load and merge all DEG data files
deg_files <- list.files(
  path = DEG_DATA_DIR,
  pattern = "DEG_.+\\.csv",
  full.names = TRUE
)
# Exclude the already merged file to avoid duplication
deg_files <- grep("merged_DEG_all_tissues.csv", deg_files, invert = TRUE, value = TRUE)

all_deg_data <- map_dfr(deg_files, read_csv)

# --- 2.2: Calculate recovery scores at the individual gene level
deg_recovery_details <- all_deg_data %>%
  # Filter for significant genes to reduce noise
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  select(Gene_symbol, Cell_type, Tissue, Comparison, avg_log2FC) %>%
  # Pivot to get comparisons in columns
  pivot_wider(
    names_from = Comparison,
    values_from = avg_log2FC
  ) %>%
  # Calculate recovery scores for each gene
  mutate(
    min_change_threshold = 0.25,
    nABX_recovery = nABXvsSPF / GFvsSPF,
    ExGF_recovery = ExGFvsGF / (-GFvsSPF),
    
    # Handle cases where the initial change (GFvsSPF) is too small
    nABX_recovery = ifelse(abs(GFvsSPF) < min_change_threshold, NA, nABX_recovery),
    ExGF_recovery = ifelse(abs(GFvsSPF) < min_change_threshold, NA, ExGF_recovery)
  ) %>%
  # Remove genes where scores couldn't be calculated
  filter(!is.na(nABX_recovery) | !is.na(ExGF_recovery))

# --- 2.3: Summarize recovery at the cell type level
recovery_data_deg <- deg_recovery_details %>%
  group_by(Cell_type, Tissue) %>%
  summarise(
    n_deg = n(),
    # Use mean of scores across all DEGs in that cell type
    nABX_recovery_score = mean(nABX_recovery, na.rm = TRUE),
    ExGF_recovery_score = mean(ExGF_recovery, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Classify recovery based on the average scores
  mutate(
    recovery_class = case_when(
      n_deg < 2 ~ "No change", # Too few DEGs to make a call
      !is.na(nABX_recovery_score) & nABX_recovery_score > 0.5 & 
        !is.na(ExGF_recovery_score) & ExGF_recovery_score > 0.5 ~ "Both reversible",
      !is.na(nABX_recovery_score) & nABX_recovery_score > 0.5 ~ "nABX reversible",
      !is.na(ExGF_recovery_score) & ExGF_recovery_score > 0.5 ~ "ExGF reversible",
      TRUE ~ "Not reversible"
    ),
    tissue_cell = paste(Tissue, Cell_type, sep = ": ")
  )

# --- 2.4: Calculate summary statistics for DEG recovery
recovery_stats_deg <- recovery_data_deg %>%
  group_by(recovery_class) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(percent = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n))

# --- 2.5: Save results
write.csv(deg_recovery_details, file.path(RESULTS_DIR, "deg_recovery_gene_level_details.csv"), row.names = FALSE)
write.csv(recovery_data_deg, file.path(RESULTS_DIR, "deg_recovery_cell_type_summary.csv"), row.names = FALSE)
write.csv(recovery_stats_deg, file.path(RESULTS_DIR, "deg_recovery_class_summary.csv"), row.names = FALSE)

cat("--- DEG Recovery Analysis Complete. Results saved. ---\n")
cat("\n--- Full pipeline complete. ---\n")
