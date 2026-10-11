###############################################################
## Step 4: annotate T/NK/ILC subclusters.
##
## Deliberately marker-driven (module-score argmax + absolute floors),
## NOT hardcoded by numeric cluster ID. An earlier version of this
## pipeline hardcoded overrides like `"7" = "gdT"` for a specific run,
## and a later change to the upstream clustering silently relabeled a
## DIFFERENT population (Cycling cells, Mki67/Top2a/Birc5 all >80%)
## as gdT, because cluster IDs are renumbered from scratch by any
## change to FindClusters/the embedding. Every rule here keys off real
## gene expression, so it stays correct even if clustering shifts
## slightly on a re-run - but ALWAYS re-check
## output/04_subcluster_scores.csv and output/03_subcluster_markers.csv
## after any upstream change before trusting the output.
##
## Three passes:
##   A) Contamination gate - drop subclusters whose non-immune module
##      score clears an ABSOLUTE floor (not "beats the T/NK score" -
##      lineage-negative ILC2 cells score negative on CoreLymphoid too,
##      so a relative rule wrongly flagged them in an earlier version).
##      CONTAM_FLOOR=0.5 sits between genuine contamination (0.8-1.0 in
##      this dataset: a near-pure epithelial/DC cluster, a near-pure
##      B-cell cluster) and the ~0.3-0.4 ambient background gut tissue
##      carries even in clean, real immune clusters.
##   B) Cluster-level subtype argmax, with its own absolute floor
##      (MIN_CALL_SCORE) - a cluster with nothing but pseudogene/
##      ribosomal/mitochondrial genes as its real FindAllMarkers
##      output still "wins" some candidate by default under pure
##      argmax; below-floor clusters are left "Unassigned" instead.
##   C) Per-cell CD4 vs CD8 rescue, symmetric, for cells in Unassigned
##      clusters (Cd4 alone vs Cd8a+Cd8b1 co-expression - not a single
##      gene at >0, which is noisy at low depth). This is NOT applied
##      to already-confidently-labeled clusters (their own marker
##      profile takes priority over this coarser fallback rule).
##   D) Per-cell Th1/Th17 rescue, applied to BOTH unspecialized CD4
##      buckets (Naive and Activated/Conventional) - not just one.
##      Checked directly: restricting this to only one of the two
##      buckets missed real cells in the other (Ifng/Il17a positivity
##      in "unchecked" buckets looked like background at the
##      cluster-mean level, ~3%, but individual cells there still
##      cleared the per-cell score floor). Treg/Tfh are deliberately
##      excluded - their own lineage marker (Foxp3, Cxcr5/Bcl6) is more
##      specific and should not be overridden by secondary Th1/Th17
##      plasticity signal.
###############################################################

library(Seurat)
library(dplyr)

DATA_DIR <- "data"; OUTPUT_DIR <- "output"

Tcells <- readRDS(file.path(DATA_DIR, "03_tnkilc_rpca.rds"))

###############################################################
## A. CONTAMINATION GATE
###############################################################

## NonImmune was originally one combined module (Epcam/Krt8/Col1a1/
## Pecam1/Cntn1/Cx3cl1). Splitting it by tissue-of-origin fixed a real
## miss: a purely epithelial cluster (Epcam/Krt8/Cldn7/Pigr/Agr2, real
## FindAllMarkers log2FC 1.8-2.0) and a purely neural/glial cluster
## (Cntn1/Foxs1/Cx3cl1, log2FC 6.6-8.9) each only carry HALF of the
## combined module's genes - averaging in the other tissue's genes
## (which read as background in a single-tissue contaminant) diluted
## both clusters below the floor that correctly caught a real Bcell/
## Myeloid-contaminated cluster in the same run (score >1.0). Scoring
## each tissue separately (and lowering the floor accordingly - see
## below) catches both real cases without re-inflating the Bcell/Mast/
## Myeloid floor that was already working.
contam_markers <- list(
  CoreLymphoid  = c("Cd3d", "Cd3e", "Cd3g", "Ptprc", "Il7r", "Thy1", "Cd2", "Nkg7", "Ncr1", "Klrb1c"),
  Bcell         = c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Ebf1", "Bank1"),
  Mast          = c("Cpa3", "Mcpt1", "Mcpt2", "Fcer1a"),
  Myeloid_DC    = c("Cd14", "Lyz2", "Itgam", "Csf1r", "Adgre1", "H2-Aa", "H2-Eb1", "Cd74", "Flt3"),
  Epithelial    = c("Epcam", "Krt8", "Krt18", "Cldn7", "Agr2", "Pigr"),
  Neural        = c("Cntn1", "Foxs1", "Cx3cl1"),
  EndoStromal   = c("Pecam1", "Col1a1", "Col1a2", "Dcn")
)
contam_markers <- lapply(contam_markers, function(g) intersect(g, rownames(Tcells)))

Tcells <- AddModuleScore(Tcells, features = contam_markers, name = "contam_", assay = "RNA")
contam_cols <- paste0("contam_", seq_along(contam_markers)); names(contam_cols) <- names(contam_markers)

per_cluster_contam <- Tcells@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(contam_cols), mean), .groups = "drop")
colnames(per_cluster_contam) <- c("cluster", names(contam_markers))
write.csv(per_cluster_contam, file.path(OUTPUT_DIR, "04_contamination_scores.csv"), row.names = FALSE)

## 0.15, not 0.3: checked against this run's actual scores (see
## output/04_contamination_scores.csv) - a cluster with real, if
## moderate, ambient-epithelial contamination (mt-Rnr1/2, mt-Nd5,
## Epcam/Krt8/Krt19/Lgals4/Lypd8 all in its own top FindAllMarkers hits)
## scored 0.23, below the original 0.3 floor, while EVERY genuinely
## clean cluster in this run scores negative on every contaminant
## category (next-highest clean-cluster score across all 5 categories
## is -0.0009) - 0.15 sits in the wide gap between those two groups and
## catches the real case without touching anything else.
CONTAM_FLOOR <- 0.15
contaminant_cols <- c("Bcell", "Mast", "Myeloid_DC", "Epithelial", "Neural", "EndoStromal")
contaminated <- per_cluster_contam$cluster[
  apply(per_cluster_contam[, contaminant_cols, drop = FALSE], 1, max) > CONTAM_FLOOR
]
if (length(contaminated) > 0) {
  n_dropped <- sum(Tcells$seurat_clusters %in% contaminated)
  cat("Descartando", length(contaminated), "subcluster(s) contaminados (",
      paste(contaminated, collapse = ", "), "):", n_dropped, "cells\n")
  Tcells <- subset(Tcells, subset = !(seurat_clusters %in% contaminated))
} else {
  cat("Sin subclusters contaminados detectados.\n")
}

###############################################################
## B. CLUSTER-LEVEL SUBTYPE ARGMAX (with absolute floor)
###############################################################

subtype_markers <- list(
  CD4_Naive     = c("Sell", "Ccr7", "Tcf7", "Lef1", "S1pr1"),
  CD4_Treg      = c("Foxp3", "Il2ra", "Ctla4", "Tnfrsf18"),
  CD4_Tfh       = c("Cxcr5", "Pdcd1", "Bcl6", "Icos"),
  CD4_Activated = c("Cd4", "Cd40lg", "Icos", "Ctla4", "Cd28"),
  CD8_Effector  = c("Cd8a", "Cd8b1", "Gzmk", "Gzma", "Ccl5", "Nkg7"),
  CD8_IEL       = c("Cd8a", "Itgae", "Itgb7"),
  CD8_TRM       = c("Cd8a", "Itgae", "Itgb7", "Gzmb", "Tbx21"),
  gdT17         = c("Trdc", "Trdv4", "Tcrg-C1", "Il17a"),
  # Klrd1 (CD94)/Cd160/Fcer1g/Xcl1 added after checking a real
  # FindAllMarkers cluster that the narrower NK panel (Ncr1/Klrb1c/
  # Klrb1b/Tyrobp/Klrk1 alone) scored too low to call (0.14, below
  # MIN_CALL_SCORE) despite an unambiguous innate-cytotoxic profile
  # (Klrd1/Cd160/Fcer1g log2FC 2.4-2.6, Xcl1 log2FC 3.0) - these are
  # standard mouse NK/innate-lymphocyte markers, not an ad hoc addition.
  NK            = c("Ncr1", "Klrb1c", "Klrb1b", "Tyrobp", "Klrk1", "Klrd1", "Cd160", "Fcer1g", "Xcl1"),
  ILC2          = c("Il17rb", "Il13", "Il1rl1", "Il4", "Gata3"),
  Cycling       = c("Mki67", "Top2a", "Pclaf", "Stmn1")
)
subtype_markers <- lapply(subtype_markers, function(g) intersect(g, rownames(Tcells)))

Tcells <- AddModuleScore(Tcells, features = subtype_markers, name = "subtype_", assay = "RNA")
sub_cols <- paste0("subtype_", seq_along(subtype_markers)); names(sub_cols) <- names(subtype_markers)

per_cluster_scores <- Tcells@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(sub_cols), mean), .groups = "drop")
colnames(per_cluster_scores) <- c("cluster", names(subtype_markers))
write.csv(per_cluster_scores, file.path(OUTPUT_DIR, "04_subcluster_scores.csv"), row.names = FALSE)
cat("\nPuntuaciones de subtipo por subcluster (REVISAR antes de confiar):\n")
print(as.data.frame(per_cluster_scores))

MIN_CALL_SCORE <- 0.35
best_call <- apply(per_cluster_scores[, names(subtype_markers), drop = FALSE], 1, function(row) {
  if (max(row) > MIN_CALL_SCORE) names(subtype_markers)[which.max(row)] else "Unassigned"
})
names(best_call) <- as.character(per_cluster_scores$cluster)

## ---------------------------------------------------------------
## MANUAL OVERRIDE HOOK: edit after reviewing output/04_subcluster_scores.csv
## and output/03_subcluster_markers.csv against the CURRENT run's cluster
## IDs (they are NOT stable across reruns - see header note). Empty by
## default.
## ---------------------------------------------------------------
OVERRIDE_LABELS <- c(
  # "3" = "CD4_Tfh",
)
for (id in names(OVERRIDE_LABELS)) best_call[id] <- OVERRIDE_LABELS[[id]]

Tcells$ClusterFull <- unname(best_call[as.character(Tcells$seurat_clusters)])

cat("\nDistribucion tras el paso B (antes de los rescates C y D):\n")
print(table(Tcells$ClusterFull))

###############################################################
## C. PER-CELL CD4/CD8 RESCUE (symmetric, for Unassigned cells only)
###############################################################

cd4  <- GetAssayData(Tcells, assay = "RNA", layer = "data")["Cd4", ]
cd8a <- GetAssayData(Tcells, assay = "RNA", layer = "data")["Cd8a", ]
cd8b1<- GetAssayData(Tcells, assay = "RNA", layer = "data")["Cd8b1", ]
is_cd8 <- cd8a > 0 & cd8b1 > 0 & cd4 == 0
is_cd4 <- cd4 > 0 & !(cd8a > 0 & cd8b1 > 0)

is_unassigned <- Tcells$ClusterFull == "Unassigned"
new_full <- Tcells$ClusterFull
new_full[is_unassigned & is_cd8] <- "CD8_Conventional"
new_full[is_unassigned & is_cd4] <- "CD4_Conventional"
Tcells$ClusterFull <- new_full

cat("\nDistribucion tras el rescate CD4/CD8 (paso C):\n")
print(table(Tcells$ClusterFull))

###############################################################
## D. PER-CELL Th1/Th17 RESCUE (both unspecialized CD4 buckets)
###############################################################

cd4_sub_markers <- list(
  CD4_Th1  = c("Tbx21", "Ifng", "Cxcr3", "Stat4"),
  CD4_Th17 = c("Rorc", "Il17a", "Il17f", "Il23r")
)
cd4_sub_markers <- lapply(cd4_sub_markers, function(g) intersect(g, rownames(Tcells)))

Tcells <- AddModuleScore(Tcells, features = cd4_sub_markers, name = "cd4sub_", assay = "RNA")
cd4_cols <- paste0("cd4sub_", seq_along(cd4_sub_markers)); names(cd4_cols) <- names(cd4_sub_markers)

SCORE_FLOOR <- 0.15
is_unspecialized <- Tcells$ClusterFull %in% c("CD4_Naive", "CD4_Activated", "CD4_Conventional")
scores <- Tcells@meta.data[is_unspecialized, cd4_cols, drop = FALSE]
colnames(scores) <- names(cd4_sub_markers)
current_label <- Tcells$ClusterFull[is_unspecialized]

best_sub <- character(length(current_label))
for (i in seq_along(current_label)) {
  row <- as.numeric(scores[i, ])
  best_sub[i] <- if (max(row) > SCORE_FLOOR) names(cd4_sub_markers)[which.max(row)] else current_label[i]
}
new_full <- Tcells$ClusterFull
new_full[is_unspecialized] <- best_sub
Tcells$ClusterFull <- new_full

###############################################################
## MINIMUM POPULATION SIZE - a real safety net, not just a display
## filter. The per-cell CD4/CD8 rescue (step C) requires 2-gene
## co-detection for the CD8 call (Cd8a>0 & Cd8b1>0) but only 1-gene
## detection for the CD4 call (Cd4>0) - checked directly: on one run,
## that asymmetry let 282 real leftover cells become CD4_Conventional
## but only 1 cell (pure dropout luck, out of ~700 mostly low-quality
## "Unassigned" cells) happened to co-detect both CD8 genes and became
## "CD8_Conventional". A single cell is not a population - no
## variance, no possible statistics, and every downstream figure/test
## quietly treated it as if it were one. Any ClusterFull label with
## fewer than MIN_POPULATION_SIZE cells total (across ALL mice/groups,
## not per mouse - see MIN_CELLS_PER_SAMPLE in 05_composition_figure.R
## for the per-mouse equivalent) is folded back into "Unassigned"
## here, at the source, so no downstream script can accidentally
## report a 1-cell "population" as if it were real.
###############################################################

MIN_POPULATION_SIZE <- 10
pop_counts <- table(Tcells$ClusterFull)
too_small <- names(pop_counts)[pop_counts > 0 & pop_counts < MIN_POPULATION_SIZE]
if (length(too_small) > 0) {
  cat("\nFolding", length(too_small), "population(s) with <", MIN_POPULATION_SIZE,
      "total cells back into Unassigned (not a real population, just rescue noise):\n")
  print(pop_counts[too_small])
  Tcells$ClusterFull[Tcells$ClusterFull %in% too_small] <- "Unassigned"
}

###############################################################
## ClusterMidway (TCD4 / TCD8 / other)
###############################################################

Tcells$ClusterMidway <- dplyr::case_when(
  Tcells$ClusterFull %in% c("CD4_Naive", "CD4_Activated", "CD4_Conventional", "CD4_Th1", "CD4_Th17", "CD4_Treg", "CD4_Tfh") ~ "TCD4",
  Tcells$ClusterFull %in% c("CD8_Effector", "CD8_IEL", "CD8_TRM", "CD8_Conventional") ~ "TCD8",
  Tcells$ClusterFull == "Cycling" & is_cd8 ~ "TCD8",
  Tcells$ClusterFull == "Cycling" & is_cd4 ~ "TCD4",
  TRUE ~ "other"
)

Tcells$Grupo <- factor(Tcells$condition, levels = c("Normal", "Colitis", "CAC"))

cat("\n=== Distribucion FINAL ===\n")
print(table(Tcells$ClusterFull))
cat("\nClusterMidway x Grupo:\n")
print(table(Tcells$ClusterMidway, Tcells$Grupo))

cat("\nTCD8 por muestra:\n")
t <- as.data.frame(table(Tcells$sample[Tcells$ClusterMidway == "TCD8"]))
colnames(t) <- c("sample", "n"); print(t[order(-t$n), ], row.names = FALSE)
cat("Samples with >=20 cells TCD8:", sum(t$n >= 20), "of 14\n")

cat("\nTCD4 por muestra:\n")
t2 <- as.data.frame(table(Tcells$sample[Tcells$ClusterMidway == "TCD4"]))
colnames(t2) <- c("sample", "n"); print(t2[order(-t2$n), ], row.names = FALSE)
cat("Samples with >=20 cells TCD4:", sum(t2$n >= 20), "of 14\n")

markers_final <- FindAllMarkers(
  SetIdent(Tcells, value = "ClusterFull"),
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, max.cells.per.ident = 2000
)
write.csv(markers_final, file.path(OUTPUT_DIR, "04_final_population_markers.csv"), row.names = FALSE)

saveRDS(Tcells, file.path(DATA_DIR, "GSE264607_Tcells_final_con_labels.rds"))
cat("\nGuardado:", file.path(DATA_DIR, "GSE264607_Tcells_final_con_labels.rds"), "\n")
cat("Revisa output/04_final_population_markers.csv para verificar cada poblacion antes de generar figuras.\n")
