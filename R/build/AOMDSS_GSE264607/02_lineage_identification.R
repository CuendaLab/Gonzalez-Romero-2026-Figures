###############################################################
## Step 2: identify the T/NK/ILC compartment via canonical marker
## module scores, with an ABSOLUTE floor on the winning score - not
## just "highest among the 7 candidates".
##
## Why a floor: an earlier version of this pipeline used pure argmax
## (whichever lineage scores highest, no matter how weak) and it
## produced a confident-looking "NK" call for a cluster whose real
## FindAllMarkers top genes were Cd74/H2-Aa/H2-Eb1/Flt3 - textbook
## dendritic cells, not NK. Its NK score was only 0.066 (essentially
## noise) but still beat every other candidate by a hair. FLOOR = 0.3
## was picked by checking the actual score distribution: every
## genuine Tcell/NK cluster in this dataset scores 0.75-1.3+, while
## contaminating/ambiguous clusters top out under 0.15. There is a
## clean gap between those two groups at this resolution - re-check
## output/02_lineage_scores.csv if you change the clustering upstream.
###############################################################

library(Seurat)
library(dplyr)

DATA_DIR <- "data"; OUTPUT_DIR <- "output"

combined <- readRDS(file.path(DATA_DIR, "01_combined_integrated.rds"))

lineage_markers <- list(
  Tcell      = c("Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd5", "Trbc1", "Trbc2"),
  NK         = c("Nkg7", "Klrb1c", "Ncr1", "Klrk1", "Gzma"),
  Bcell      = c("Cd19", "Ms4a1", "Cd79a", "Cd79b"),
  Myeloid    = c("Cd14", "Lyz2", "Itgam", "Csf1r", "Adgre1"),
  Epithelial = c("Epcam", "Krt8", "Krt18"),
  Stromal    = c("Col1a1", "Col1a2", "Pecam1", "Dcn"),
  # Mast cells score near-zero on every OTHER lineage here (none of
  # Myeloid's markers are mast-specific), so without their own
  # candidate they fell through to "Tcell" by default in an earlier
  # version of this pipeline.
  Mast       = c("Cpa3", "Mcpt1", "Mcpt2", "Fcer1a")
)
lineage_markers <- lapply(lineage_markers, function(g) intersect(g, rownames(combined)))

combined <- AddModuleScore(combined, features = lineage_markers, name = "lineage_", assay = "RNA")
score_cols <- paste0("lineage_", seq_along(lineage_markers))
names(score_cols) <- names(lineage_markers)

per_cluster_scores <- combined@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), n_cells = n(), .groups = "drop")
colnames(per_cluster_scores) <- c("cluster", names(lineage_markers), "n_cells")

best_lineage <- apply(per_cluster_scores[, names(lineage_markers)], 1, function(row) names(lineage_markers)[which.max(row)])
names(best_lineage) <- as.character(per_cluster_scores$cluster)
per_cluster_scores$best_lineage <- unname(best_lineage)

write.csv(per_cluster_scores, file.path(OUTPUT_DIR, "02_lineage_scores.csv"), row.names = FALSE)
cat("Lineage scores per cluster:\n")
print(as.data.frame(per_cluster_scores))

FLOOR <- 0.3
tnkilc_clusters <- per_cluster_scores$cluster[
  best_lineage[as.character(per_cluster_scores$cluster)] %in% c("Tcell", "NK") &
    pmax(per_cluster_scores$Tcell, per_cluster_scores$NK) > FLOOR
]
cat("\nClusters T/NK/ILC (floor", FLOOR, "):", paste(tnkilc_clusters, collapse = ", "), "\n")

Tcells <- subset(combined, subset = seurat_clusters %in% as.character(tnkilc_clusters))
cat("T/NK/ILC compartment:", ncol(Tcells), "cells\n")

saveRDS(Tcells, file.path(DATA_DIR, "02_tnkilc_raw.rds"))
cat("\nGuardado:", file.path(DATA_DIR, "02_tnkilc_raw.rds"), "\n")
