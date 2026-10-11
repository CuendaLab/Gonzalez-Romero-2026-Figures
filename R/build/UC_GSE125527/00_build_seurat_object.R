#############################################################
## GSE125527 - Smillie et al. 2019 (Cell)
## scRNA-seq of CD45+ immune cells in ulcerative colitis (UC)
## Building the Seurat object
#############################################################

## --- 0. Paquetes ------------------------------------------------
pkgs <- c("Matrix", "data.table", "Seurat", "R.utils")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst)

library(Matrix)
library(data.table)
library(Seurat)

## --- 1. Download the files (only if not already present) --------
options(timeout = 1200)  # the sparse file is ~188 MB; the default 60 s timeout is not enough

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125527/suppl/"
data_dir <- "GSE125527_data"
dir.create(data_dir, showWarnings = FALSE)

## Returns the path to the file, whether it is a .csv (already unpacked,
## e.g. if the browser unpacked it on download) or a .csv.gz
resolve_file <- function(fname_gz) {
  fname_csv <- sub("\\.gz$", "", fname_gz)
  path_csv  <- file.path(data_dir, fname_csv)
  path_gz   <- file.path(data_dir, fname_gz)

  if (file.exists(path_csv)) return(path_csv)
  if (!file.exists(path_gz)) {
    for (attempt in 1:3) {
      message("Descargando ", fname_gz, " (intento ", attempt, ")")
      ok <- tryCatch({
        download.file(paste0(base_url, fname_gz), path_gz, mode = "wb")
        TRUE
      }, error = function(e) {
        message("  fallo: ", conditionMessage(e))
        FALSE
      })
      if (ok) break
      if (file.exists(path_gz)) file.remove(path_gz)
    }
    if (!file.exists(path_gz)) stop("No se pudo descargar ", fname_gz, " tras 3 intentos.")
  }
  path_gz
}

files <- c(
  "GSE125527_UMI_cell_table_sparse.csv.gz",
  "GSE125527_gene_id_rownames.csv.gz",
  "GSE125527_cell_id_colnames.csv.gz",
  "GSE125527_cell_metadata.csv.gz",
  "GSE125527_Bcell_cluster.csv.gz",
  "GSE125527_Tcell_cluster.csv.gz",
  "GSE125527_oldPatientId-newPatientId.csv.gz"
)

paths <- setNames(lapply(files, resolve_file), files)

## --- 2. Load genes and cells (they define the matrix dimnames) --
genes <- fread(paths[["GSE125527_gene_id_rownames.csv.gz"]], header = FALSE)[[1]]
cells <- fread(paths[["GSE125527_cell_id_colnames.csv.gz"]], header = FALSE)[[1]]

## --- 3. Load the sparse matrix (triplet format: row,col,value) --
triplet <- fread(paths[["GSE125527_UMI_cell_table_sparse.csv.gz"]])
# columns: row (gene index, 1-based), col (cell index, 1-based), value (UMI count)

mat <- sparseMatrix(
  i    = triplet$row,
  j    = triplet$col,
  x    = triplet$value,
  dims = c(length(genes), length(cells)),
  dimnames = list(genes, cells)
)

## --- 4. Per-cell metadata ---------------------------------------
## The CSV was written with R row names: it has 6 columns but the
## header lists only 5 (the row-index name is missing).
meta <- fread(paths[["GSE125527_cell_metadata.csv.gz"]], header = TRUE)
setnames(meta, c("row_idx", "cell_id", "patient_assignment",
                  "tissue_assignment", "disease_assignment", "celltype"))

## Align by barcode (more robust than assuming the same order in both files)
idx <- match(cells, meta$cell_id)
if (anyNA(idx)) stop(sum(is.na(idx)), " cells from cell_id_colnames.csv are missing en cell_metadata.csv")

meta_df <- as.data.frame(meta[idx, -"row_idx"])
rownames(meta_df) <- meta_df$cell_id

## --- 5. Fine annotation of B and T subclusters (separate columns) --
bcell <- fread(paths[["GSE125527_Bcell_cluster.csv.gz"]])
tcell <- fread(paths[["GSE125527_Tcell_cluster.csv.gz"]])

meta_df$Bcell_cluster <- bcell$Cluster_id[match(rownames(meta_df), bcell$cell_id)]
meta_df$Tcell_cluster <- tcell$Cluster_id[match(rownames(meta_df), tcell$cell_id)]

## --- 6. Create the Seurat object ("raw" checkpoint: no QC, no normalisation) -
seu <- CreateSeuratObject(
  counts    = mat,
  meta.data = meta_df,
  project   = "GSE125527_colitis",
  min.cells = 3,
  min.features = 200
)

saveRDS(seu, file.path(data_dir, "GSE125527_seurat_raw.rds"))
cat("Saved:", file.path(data_dir, "GSE125527_seurat_raw.rds"), "\n")
cat("dims:", paste(dim(seu), collapse = " x "), "\n")
print(table(seu$celltype))

## QC (nFeature/percent.mt filtering), normalisation, PCA/UMAP/clustering
## and the CD4/CD8 subsets are done in the analysis scripts (01_...),
## as in the CRC project: every figure starts from this raw checkpoint
## and applies its own processing, instead of baking a single
## one-size-fits-all normalisation into this script.
