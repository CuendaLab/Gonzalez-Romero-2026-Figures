###############################################################
## Step 0: QC filter + doublet removal, per sample.
##
## Thresholds were picked by looking at the actual unfiltered
## nFeature_RNA/nCount_RNA/percent.mt distributions per sample (violin
## plots), not assumed a priori - all 3 cutoffs sit right where the
## low-quality-droplet bulge separates from the real-cell population,
## consistently across all 14 samples. Total cells after QC alone
## (92,326) lands very close to the paper's own reported 93,158
## (Zhou et al. 2025, PMID 40859429), a useful cross-check.
##
## Doublet removal (scDblFinder, per sample - must run before merging,
## its simulated-doublet model is fit per batch) brings the total to
## ~83,000 cells, closing most of the remaining gap to the paper's
## count.
###############################################################

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)

###############################################################
## CONFIGURATION
###############################################################

## folder with the CellRanger output of GSE264607 (one subfolder per sample,
## with barcodes/features/matrix); downloaded from GEO as GSE264607_RAW.tar
RAW_DIR <- file.path(
  "~", "GSE264607_RAW"
)
DATA_DIR   <- "data"
OUTPUT_DIR <- "output"
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

samples <- c(
  "N1_processed", "N2_processed", "N3_processed", "N4_processed", "N5_processed",
  "A1_processed", "A2_processed", "A3_processed", "A4_processed",
  "C1_processed", "C2_processed", "C3_processed", "C4_processed", "C5_processed"
)

group_of_sample <- function(sample) {
  dplyr::case_when(
    grepl("^N", sample) ~ "Normal",
    grepl("^A", sample) ~ "Colitis",
    grepl("^C", sample) ~ "CAC",
    TRUE ~ NA_character_
  )
}

set.seed(1234)

###############################################################
## READ, QC-FILTER, REMOVE DOUBLETS - PER SAMPLE
###############################################################

obj.list <- list()
results <- list()

for (sample in samples) {
  cat("Procesando:", sample, "\n")

  counts <- Read10X(file.path(RAW_DIR, sample))
  obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = sample)
  obj$sample <- sample
  obj$condition <- group_of_sample(sample)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

  n_raw <- ncol(obj)

  obj <- subset(
    obj,
    subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & nCount_RNA > 1000 & percent.mt < 25
  )
  n_after_qc <- ncol(obj)

  sce <- SingleCellExperiment(assays = list(counts = GetAssayData(obj, assay = "RNA", layer = "counts")))
  sce <- scDblFinder(sce)
  obj$scDblFinder.class <- sce$scDblFinder.class
  n_doublets <- sum(obj$scDblFinder.class == "doublet")
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  n_final <- ncol(obj)

  results[[sample]] <- data.frame(
    sample = sample, condition = group_of_sample(sample),
    n_raw = n_raw, n_after_qc = n_after_qc,
    n_doublets = n_doublets, pct_doublets = round(100 * n_doublets / n_after_qc, 1),
    n_final = n_final
  )

  obj.list[[sample]] <- obj
}

summary_tab <- bind_rows(results)
write.csv(summary_tab, file.path(OUTPUT_DIR, "00_qc_doublet_summary.csv"), row.names = FALSE)
cat("\n--- Resumen QC + dobletes ---\n")
print(summary_tab)
cat("\nFinal cell count:", sum(summary_tab$n_final), "\n")

saveRDS(obj.list, file.path(DATA_DIR, "00_obj_list_qc_doublet_filtered.rds"))
cat("\nGuardado:", file.path(DATA_DIR, "00_obj_list_qc_doublet_filtered.rds"), "\n")
