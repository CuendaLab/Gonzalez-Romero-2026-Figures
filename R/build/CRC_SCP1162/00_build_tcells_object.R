###############################################################
## Build Tcells_final_con_UMAP_y_labels.rds
##
## Subsets the T/NK/ILC compartment (26 subpopulations, ClusterFull
## cTNI01-26) out of the full CRC atlas and attaches the portal UMAP
## coordinates, cluster labels, and Normal/Tumor group - producing the
## object consumed by fig1_tcell_composition.R.
##
## Source data: Broad Single Cell Portal study SCP1162 (Pelka et al.,
## Cell 2021 - colorectal cancer atlas). See data/README.md for the
## three required input files and where to obtain them.
##
## IMPORTANT (normalization): scp1162_raw_seurat.rds already contains
## the portal's normalized matrix (decimal values). This script does
## NOT normalize - it only subsets and copies values as-is. Do not add
## NormalizeData()/log1p() here or in the figure script.
###############################################################

library(Seurat)
library(dplyr)

###############################################################
## CONFIGURATION
###############################################################

RAW_RDS     <- file.path("data", "raw", "scp1162_raw_seurat.rds")
PORTAL_TSV  <- file.path("data", "raw", "crc10x_tSNE_cl_TNKILC.tsv")
META_TSV    <- file.path("data", "raw", "metatable_v3_fix_v3.tsv")
OUTPUT_RDS  <- file.path("data", "Tcells_final_con_UMAP_y_labels.rds")

for (f in c(RAW_RDS, PORTAL_TSV, META_TSV)) {
  if (!file.exists(f)) {
    stop(
      "Required input not found: ", f,
      "\nSee data/README.md for where to obtain the SCP1162 files."
    )
  }
}

dir.create(dirname(OUTPUT_RDS), showWarnings = FALSE, recursive = TRUE)

###############################################################
## 1. LOAD FULL scRNA-seq OBJECT
###############################################################

obj <- readRDS(RAW_RDS)

###############################################################
## 2. LOAD PORTAL T/NK/ILC CLUSTER TABLE
###############################################################

TN <- read.table(PORTAL_TSV, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

TN <- subset(TN, NAME != "TYPE")
TN <- TN[!duplicated(TN$NAME), ]

###############################################################
## 3. SUBSET T/NK/ILC CELLS
###############################################################

Tcells <- subset(obj, cells = TN$NAME)

###############################################################
## 4. ATTACH PORTAL UMAP COORDINATES
###############################################################

TN$X <- as.numeric(gsub(",", ".", TN$X))
TN$Y <- as.numeric(gsub(",", ".", TN$Y))

idx <- match(colnames(Tcells), TN$NAME)

emb <- cbind(
  UMAP_1 = TN$X[idx],
  UMAP_2 = TN$Y[idx]
)
rownames(emb) <- colnames(Tcells)

Tcells[["umap"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "UMAP_",
  assay = DefaultAssay(Tcells)
)

###############################################################
## 5. ATTACH ClusterFull / ClusterMidway
###############################################################

if ("ClusterFull" %in% colnames(TN)) {
  Tcells$ClusterFull <- TN$ClusterFull[idx]
}
if ("ClusterMidway" %in% colnames(TN)) {
  Tcells$ClusterMidway <- TN$ClusterMidway[idx]
}

###############################################################
## 6. HUMAN-READABLE LABELS FOR THE 26 SUBPOPULATIONS
##    (ClusterFull -> CT_label)
###############################################################

label_map <- c(
  "cTNI01 (CD4+ IL7R+)"          = "CD4 IL7R+",
  "cTNI02 (CD4+ IL7R+SELL+)"     = "CD4 IL7R+ SELL+",
  "cTNI03 (CD4+ IL7R+HSP+)"      = "CD4 IL7R+ HSP+",
  "cTNI04 (CD4+ IL7R+CCL5+)"     = "CD4 IL7R+ CCL5+",
  "cTNI05 (CD4+ IL17+)"          = "CD4 Th17",
  "cTNI06 (CD4+ TFH)"            = "CD4 TFH",
  "cTNI07 (CD4+ CXCL13+)"        = "CD4 CXCL13+",
  "cTNI08 (CD4+ Treg)"           = "Treg",
  "cTNI09 (CD4+ Treg prolif)"    = "Treg prolif",

  "cTNI10 (CD8+ IL7R+)"          = "CD8 IL7R+",
  "cTNI11 (CD8+GZMK+)"           = "CD8 GZMK+",
  "cTNI12 (CD8+ IL7R+)"          = "CD8 IL7R+ (mem)",
  "cTNI13 (CD8+ T IL17+)"        = "CD8 IL17+",
  "cTNI14 (CD8+ CXCL13+)"        = "CD8 CXCL13+",
  "cTNI15 (CD8+ CXCL13+ HSP+)"   = "CD8 CXCL13+ HSP+",
  "cTNI16 (CD8+ CXCL13+ prolif)" = "CD8 CXCL13+ prolif",

  "cTNI17 (gd-like T)"           = "γδ T",
  "cTNI18 (gd-like T PDCD1+)"    = "γδ T PDCD1+",
  "cTNI19 (gd-like T prolif)"    = "γδ T prolif",

  "cTNI20 (PLZF+ T)"             = "PLZF+ T",
  "cTNI21 (PLZF+ T prolif)"      = "PLZF+ T prolif",

  "cTNI22 (cTNI22)"              = "cTNI22",
  "cTNI23 (NK CD16A+)"           = "NK CD16A+",
  "cTNI24 (NK GZMK+)"            = "NK GZMK+",
  "cTNI25 (NK XCL1+)"            = "NK XCL1+",
  "cTNI26 (ILC3)"                = "ILC3"
)

if ("ClusterFull" %in% colnames(Tcells[[]])) {
  Tcells$CT_label <- ifelse(
    Tcells$ClusterFull %in% names(label_map),
    label_map[Tcells$ClusterFull],
    Tcells$ClusterFull
  )
}

###############################################################
## 7. NORMAL / TUMOR GROUP FROM SAMPLE METADATA
###############################################################

meta <- read.delim(META_TSV, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(meta) <- meta$NAME

grupo_raw <- meta$SpecimenType[match(colnames(Tcells), rownames(meta))]

Grupo <- rep(NA_character_, length(grupo_raw))
Grupo[grupo_raw %in% c("N", "Normal", "Heal", "Healthy")] <- "Normal"
Grupo[grupo_raw %in% c("T", "Tumor")] <- "Tumor"
Tcells$Grupo <- factor(Grupo, levels = c("Normal", "Tumor"))

###############################################################
## 8. SANITY CHECKS
###############################################################

cat("T/NK/ILC cells:", ncol(Tcells), "\n")
cat("Subpopulations (ClusterFull):", length(unique(na.omit(Tcells$ClusterFull))), "\n")
print(table(Tcells$Grupo, useNA = "ifany"))

###############################################################
## 9. SAVE
###############################################################

saveRDS(Tcells, OUTPUT_RDS)
cat("Saved:", OUTPUT_RDS, "\n")
