###############################################################
## Generates every figure of this repository, in order
##
## Usage:  Rscript run_all.R        (from the root of the repository)
##
## Before running it, edit the four data paths at the top of
## R/00_paths.R.
###############################################################

scripts <- c(
  ## ---- Figure 8 ----
  "R/01_fig8AB_CRC_UMAP_composition.R",          # 8A, 8B
  "R/03_fig8C_suppl8B_CRC_per_patient_scores.R", # 8C y S8B
  "R/04_fig8D_CRC_CD4_vs_CD8_delta.R",           # 8D   (uses the table from 03)
  "R/05_fig8E_suppl8A_CRC_CD4_allgenes_grid.R",  # 8E y S8A
  "R/07a_gene_level_concordance_table.R",        # table for 8F
  "R/07b_fig8F_direction_heatmap.R",             # 8F   (uses the table from 07a)
  ## ---- Suppl. Figure 7 (AOM/DSS) ----
  "R/08_suppl7CD_AOM_UMAP_composition.R",        # S7C, S7D
  "R/09_suppl7E_AOM_CD4_allgenes_grid.R",        # S7E
  ## ---- Suppl. Figure 8 (CRC) ----
  "R/06_suppl8C_CRC_CD8_allgenes_4columns.R",    # S8C
  ## ---- Suppl. Figure 9 (colitis and cross-dataset summary) ----
  "R/10_suppl9A_UC_gene_heatmaps.R",             # S9A
  "R/11a_proteomics_mLN_module_delta.R",         # proteomics table
  "R/11b_crossdataset_table.R",                  # cross-dataset summary table
  "R/11c_suppl9B_crossdataset_dotmatrix.R"       # S9B
)

for (s in scripts) {
  cat("\n==========", s, "==========\n")
  source(s, echo = FALSE)
  rm(list = setdiff(ls(), c("scripts", "s")))   # each script starts clean
  invisible(gc())
}
## record the exact environment of the run (R and package versions)
writeLines(capture.output(sessionInfo()), file.path("tables", "sessionInfo.txt"))

cat("\nDone. Figures in figures/, tables in tables/ (environment in tables/sessionInfo.txt).\n")
