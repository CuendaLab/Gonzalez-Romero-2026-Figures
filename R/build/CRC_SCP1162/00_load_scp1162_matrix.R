###############################################################
## Step 0: build the full-atlas Seurat object from the portal files
##
## The Broad Single Cell Portal distributes study SCP1162 as an expression
## matrix in Matrix Market format plus its gene and barcode tables. This
## script only reads those three files and saves them as a Seurat object,
## which is what step 01 then subsets.
##
## Download from https://singlecell.broadinstitute.org/single_cell/study/SCP1162
## into DATA_DIR:
##   matrix.mtx.gz, matrix.genes.tsv, matrix.barcodes.tsv
###############################################################

library(Seurat)

DATA_DIR <- "~/SCP1162"      # folder with the three files downloaded from the portal

mat <- ReadMtx(
  mtx      = file.path(DATA_DIR, "matrix.mtx.gz"),
  features = file.path(DATA_DIR, "matrix.genes.tsv"),
  cells    = file.path(DATA_DIR, "matrix.barcodes.tsv")
)

obj <- CreateSeuratObject(mat, project = "SCP1162")

saveRDS(obj, file = file.path(DATA_DIR, "scp1162_raw_seurat.rds"))
cat("Saved:", file.path(DATA_DIR, "scp1162_raw_seurat.rds"), "\n")
