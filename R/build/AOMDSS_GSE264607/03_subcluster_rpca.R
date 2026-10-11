###############################################################
## Step 3: fine-subcluster the T/NK/ILC compartment using RPCA
## integration (NOT Harmony).
##
## Why not Harmony here: per-mouse cell counts within just the T/NK/
## ILC compartment range from ~12 to ~1900 cells - a ~250x imbalance
## (vs ~5x for the whole object in step 01). Harmony's per-batch
## correction becomes unstable at that imbalance: checked directly by
## looking at where clearly Cd8a+Cd8b1+ cells (real CD8 T cells, present
## at comparable rates in every one of the 14 mice - confirmed by raw
## gene expression, independent of any clustering) actually landed -
## with fresh per-subset Harmony, the CD8 cluster was reproducibly
## built almost entirely from the 3 biggest mice (only 3-9 of 14 mice
## cleared a 20-cell floor across repeated attempts, several different
## resolutions, and even after regressing out nCount_RNA as a
## covariate - none of which fixed it). The same big/small split also
## silently fractured Treg, gammadelta T, and ILC2 into two clusters
## each (one per batch-of-origin, not two real subtypes - confirmed by
## checking that each pair partitioned samples 100/0 between the 3
## biggest mice and the other 11).
##
## RPCA integration (more conservative than Harmony, does not try to
## force-merge batches with weak anchor support) resolved this:
## 9-11 of 14 mice clear the 20-cell floor for CD8, and Treg/gdT/ILC2
## each stay as ONE cluster. k.weight is lowered from Seurat's default
## (100) because the smallest sample has ~30-40 cells in this
## compartment - the default errors out ("k.weight larger than
## smallest object").
##
## NPCS is capped at 20 (not Seurat's usual 30-50) for the same reason:
## RPCAIntegration computes its own per-sample PCA internally, and
## RunPCA requires npcs strictly less than the sample's own cell count
## - the smallest sample in this run has only 30 cells (see
## output/02_lineage_scores.csv / per-sample counts), so npcs=30
## (the original value here) fails outright ("max(nu, nv) must be
## strictly less than min(nrow(A), ncol(A))"). 20 leaves a 10-cell
## margin below the smallest sample. NPCS=15 was tried first and ran
## without error, but produced clusters coarse enough to merge
## CD8_IEL/CD8_TRM into CD8_Effector and CD4_Tfh into CD4_Treg/
## CD4_Activated (checked against output/03_subcluster_markers.csv -
## none of those clusters' own top markers showed the Itgae/Itgb7 or
## Cxcr5/Bcl6 signal on their own) - 20 recovers that resolution.
##
## k.filter=NA, k.score=10: Seurat's anchor-finding defaults
## (k.filter=200, k.score=30) assume every dataset has hundreds of
## cells. With the smallest sample at 30 cells, the neighbor search
## those defaults rely on runs out of neighbors partway through and
## crashes with a low-level matrix-assignment error ("numero de
## elements to replace is not a multiple of the replacement length)
## reemplazo"). k.filter=NA skips that filtering pass entirely
## (Seurat's own documented recommendation when the smallest dataset
## is small); k.score=10 keeps anchor scoring's kNN below the
## smallest sample's cell count with margin to spare.
###############################################################

library(Seurat)
options(future.globals.maxSize = 4 * 1024^3)

DATA_DIR <- "data"; OUTPUT_DIR <- "output"
NPCS <- 20

Tcells <- readRDS(file.path(DATA_DIR, "02_tnkilc_raw.rds"))
DefaultAssay(Tcells) <- "RNA"

Tcells[["RNA"]] <- split(Tcells[["RNA"]], f = Tcells$sample)
Tcells <- NormalizeData(Tcells, verbose = FALSE)
Tcells <- FindVariableFeatures(Tcells, nfeatures = 2000, verbose = FALSE)
Tcells <- ScaleData(Tcells, verbose = FALSE)
Tcells <- RunPCA(Tcells, npcs = NPCS, verbose = FALSE)

Tcells <- IntegrateLayers(
  object = Tcells, method = RPCAIntegration, orig.reduction = "pca",
  new.reduction = "integrated.rpca", k.weight = 20, dims = 1:NPCS,
  k.filter = NA, k.score = 10, verbose = FALSE
)
Tcells[["RNA"]] <- JoinLayers(Tcells[["RNA"]])

Tcells <- FindNeighbors(Tcells, reduction = "integrated.rpca", dims = 1:NPCS)
Tcells <- FindClusters(Tcells, resolution = 0.8)
Tcells <- RunUMAP(Tcells, reduction = "integrated.rpca", dims = 1:NPCS)

cat("Subclusters:", length(unique(Idents(Tcells))), "\n")

markers <- FindAllMarkers(Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4, max.cells.per.ident = 3000)
write.csv(markers, file.path(OUTPUT_DIR, "03_subcluster_markers.csv"), row.names = FALSE)

saveRDS(Tcells, file.path(DATA_DIR, "03_tnkilc_rpca.rds"))
cat("\nGuardado:", file.path(DATA_DIR, "03_tnkilc_rpca.rds"), "\n")
cat("Revisa output/03_subcluster_markers.csv antes de anotar (paso 04).\n")
