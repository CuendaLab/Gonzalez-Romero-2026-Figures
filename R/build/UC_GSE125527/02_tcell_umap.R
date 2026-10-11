###############################################################
## UMAP of the T cells (GSE125527), coloured by:
##  - subcluster original (T1-T17, Smillie et al.)
##  - marker-based lineage assignment (CD4 / CD8 / MAIT / gdT)
##  - condicion (sano / enfermo)
## The per-subcluster lineage table comes from
## code/01_cd4_cd8_marker_check.R (output/cd4_cd8_classification/...).
###############################################################

library(Seurat)
library(ggplot2)

INPUT_RDS  <- file.path("GSE125527_data", "GSE125527_seurat_raw.rds")
CALL_CSV   <- file.path("output", "cd4_cd8_classification", "Tcell_cluster_marker_pct_positive.csv")
OUTPUT_DIR <- file.path("output", "tcell_umap")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

UMAP_RDS <- file.path("GSE125527_data", "GSE125527_Tcells_umap.rds")

## --- Lineage per subcluster (CD4/CD8/MAIT/gdT) -- always used, whether
## recalculado el UMAP o no -----------------------------------------------
calls <- read.csv(CALL_CSV, stringsAsFactors = FALSE)
lineage_map <- setNames(calls$call, calls$cluster)

if (file.exists(UMAP_RDS)) {
  Tcells <- readRDS(UMAP_RDS)
} else {
  seu <- readRDS(INPUT_RDS)
  Tcells <- subset(seu, subset = celltype == "T")

  Tcells$lineage <- factor(unname(lineage_map[Tcells$Tcell_cluster]),
                            levels = c("CD4", "CD8", "MAIT", "gdT"))
  Tcells$Tcell_cluster <- factor(Tcells$Tcell_cluster, levels = paste0("T", 1:17))

  ## --- Pipeline estandar hasta UMAP -----------------------------------
  Tcells <- NormalizeData(Tcells, verbose = FALSE)
  Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  Tcells <- ScaleData(Tcells, verbose = FALSE)
  Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)
  Tcells <- RunUMAP(Tcells, dims = 1:30, verbose = FALSE)

  saveRDS(Tcells, UMAP_RDS)
}

## --- Paleta ----------------------------------------------------------
lineage_cols <- c(CD4 = "#3B82C4", CD8 = "#D9534F", MAIT = "#F0AD4E", gdT = "#5CB85C")

## --- Plots -------------------------------------------------------------
p_cluster <- DimPlot(Tcells, group.by = "Tcell_cluster", label = TRUE, repel = TRUE) +
  ggtitle("Subclusters T (T1-T17)") + theme(legend.position = "none")

p_lineage <- DimPlot(Tcells, group.by = "lineage", cols = lineage_cols, label = TRUE, repel = TRUE) +
  ggtitle("Linaje (CD4 / CD8 / MAIT / gdT)")

p_disease <- DimPlot(Tcells, group.by = "disease_assignment",
                      cols = c(healthy = "#4C9F70", diseased = "#B84C4C")) +
  ggtitle("Sano vs. enfermo")

p_tissue <- DimPlot(Tcells, group.by = "tissue_assignment",
                     cols = c(PBMC = "#8E5FBF", R = "#D98A3D")) +
  ggtitle("Sangre (PBMC) vs. tejido (R = biopsia rectal)")

combo <- p_cluster + p_lineage + p_tissue + p_disease + patchwork::plot_layout(nrow = 1)

ggsave(file.path(OUTPUT_DIR, "Tcell_UMAP_cluster_lineage_tissue_disease.png"),
       combo, width = 24, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Tcell_UMAP_cluster_lineage_tissue_disease.pdf"),
       combo, width = 24, height = 6)

ggsave(file.path(OUTPUT_DIR, "Tcell_UMAP_lineage.png"), p_lineage, width = 7, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Tcell_UMAP_tissue.png"), p_tissue, width = 7, height = 6, dpi = 300)

## --- Tissue x subcluster x lineage table (checks the tissue confounding) ---
tissue_tab <- as.data.frame.matrix(table(Tcells$Tcell_cluster, Tcells$tissue_assignment))
tissue_tab$lineage <- lineage_map[rownames(tissue_tab)]
tissue_tab$pct_tissue_R <- round(100 * tissue_tab$R / (tissue_tab$R + tissue_tab$PBMC), 1)
write.csv(tissue_tab, file.path(OUTPUT_DIR, "Tcell_cluster_tissue_composition.csv"))
print(tissue_tab)

cat("Saved in", OUTPUT_DIR, "\n")
