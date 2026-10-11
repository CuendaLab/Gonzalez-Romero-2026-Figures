###############################################################
## Subset of T cells from TISSUE only (rectal biopsy, R), excluding
## PBMC. From here on the functional-programme analysis uses only
## funcionales en colitis se hace unicamente sobre este objeto.
##
## Se recalculan PCA/UMAP sobre el subset (no se reusa el embedding
## tissue cells, and the UMAP is recomputed (the previous one, built
## with blood+tissue together, is dominated by the blood/tissue axis).
##
## Tissue clusters (from the tissue_composition table):
##   CD4: T1, T2, T5, T7 (T7 = Treg, still within CD4)
##   CD8: T10, T14
## MAIT (T15) and gdT (T13, T16) were almost exclusively from blood and
## are absent (or marginal) in tissue with this subclustering.
###############################################################

library(Seurat)
library(ggplot2)

INPUT_RDS  <- file.path("GSE125527_data", "GSE125527_Tcells_umap.rds")
OUTPUT_RDS <- file.path("GSE125527_data", "GSE125527_Tcells_tissue.rds")
OUTPUT_DIR <- file.path("output", "tcell_umap_tissue")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

Tcells <- readRDS(INPUT_RDS)

Tissue <- subset(Tcells, subset = tissue_assignment == "R")
Tissue$Tcell_cluster <- droplevels(Tissue$Tcell_cluster)
Tissue$lineage <- droplevels(Tissue$lineage)

cat("Celulas T de tejido:", ncol(Tissue), "\n")
print(table(Tissue$Tcell_cluster, Tissue$lineage))
print(table(Tissue$disease_assignment))

## --- Recompute PCA/UMAP on tissue cells only ----------------------
Tissue <- FindVariableFeatures(Tissue, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Tissue <- ScaleData(Tissue, verbose = FALSE)
Tissue <- RunPCA(Tissue, npcs = 30, verbose = FALSE)
Tissue <- RunUMAP(Tissue, dims = 1:30, verbose = FALSE)

saveRDS(Tissue, OUTPUT_RDS)

## --- Plots -----------------------------------------------------------
## A few residual MAIT/gd cells remain in tissue (clusters that are
## almost all blood but hold a handful of tissue cells) -- they are
## still plotted, but are far too small to draw conclusions from.
lineage_cols <- c(CD4 = "#3B82C4", CD8 = "#D9534F", MAIT = "#F0AD4E", gdT = "#5CB85C")

p_cluster <- DimPlot(Tissue, group.by = "Tcell_cluster", label = TRUE, repel = TRUE) +
  ggtitle("Subclusters T de tejido")

p_lineage <- DimPlot(Tissue, group.by = "lineage", cols = lineage_cols, label = TRUE, repel = TRUE) +
  ggtitle("Linaje (CD4 / CD8 / MAIT / gdT) - solo tejido")

p_disease <- DimPlot(Tissue, group.by = "disease_assignment",
                      cols = c(healthy = "#4C9F70", diseased = "#B84C4C")) +
  ggtitle("Sano vs. enfermo - solo tejido")

combo <- p_cluster + p_lineage + p_disease

ggsave(file.path(OUTPUT_DIR, "Tcell_tissue_UMAP.png"), combo, width = 18, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Tcell_tissue_UMAP.pdf"), combo, width = 18, height = 6)

cat("Saved:", OUTPUT_RDS, "and plots in", OUTPUT_DIR, "\n")
