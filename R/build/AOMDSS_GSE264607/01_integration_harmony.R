###############################################################
## Step 1: merge, normalize, integrate the WHOLE object (all cell
## types) with Harmony (grouped by sample), cluster at a coarse
## resolution - just for broad lineage separation (Epithelial/
## Myeloid/T/NK/B/Stromal/Mast), not final cell-type calls.
##
## Harmony is appropriate HERE (unlike for the T/NK/ILC subset in step
## 03) because per-sample cell counts at this whole-object level are
## reasonably balanced (~2300-11900 cells/sample, ~5x range) - nothing
## like the ~250x range the T/NK/ILC compartment alone has per sample.
## See code/03_subcluster_rpca.R for why Harmony breaks down there and
## RPCA is used instead.
###############################################################

library(Seurat)
library(harmony)
library(dplyr)

DATA_DIR <- "data"; OUTPUT_DIR <- "output"

obj.list <- readRDS(file.path(DATA_DIR, "00_obj_list_qc_doublet_filtered.rds"))

combined <- merge(x = obj.list[[1]], y = obj.list[-1])
combined <- JoinLayers(combined)

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 50, verbose = FALSE)

combined <- RunHarmony(combined, group.by.vars = "sample")
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

combined$condition <- factor(combined$condition, levels = c("Normal", "Colitis", "CAC"))

cat("Total integrated cells:", ncol(combined), "\n")
cat("Numero de clusters:", length(unique(Idents(combined))), "\n")

saveRDS(combined, file.path(DATA_DIR, "01_combined_integrated.rds"))
cat("\nGuardado:", file.path(DATA_DIR, "01_combined_integrated.rds"), "\n")
