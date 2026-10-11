# Building the Seurat objects

These scripts are **not** called by `run_all.R`: they are the processing pipeline of
each dataset, run once, and they produce the object that the figure scripts then use.
Each pipeline is run **from its own dataset folder** (the scripts use relative paths
such as `data/` and `output/`), and afterwards `R/00_paths.R` is pointed at the
resulting objects.

The input data are not included here: they are public and large.

---

## Human CRC — SCP1162 / GSE178341 (Pelka et al., *Cell* 2021)

`CRC_SCP1162/00_build_tcells_object.R`

Input (Broad Single Cell Portal, study SCP1162), placed in `data/raw/`:

```
scp1162_raw_seurat.rds       full-atlas Seurat object (matrix ALREADY normalised by the portal)
crc10x_tSNE_cl_TNKILC.tsv    T/NK/ILC portal clusters and UMAP coordinates
metatable_v3_fix_v3.tsv      sample-level metadata (SpecimenType, donor_id)
```

Output: `data/Tcells_final_con_UMAP_y_labels.rds` → `CRC_RDS`.

It selects the T/NK/ILC compartment (cTNI01–26) and attaches the coordinates, the
cluster labels and the Normal/Tumour group. It does **not** normalise: the portal values
are already normalised, which is why the figure scripts do not normalise either.

## Mouse AOM/DSS — GSE264607 (Zhou et al., 2025)

Run in order:

```
AOMDSS_GSE264607/00_qc_and_doublets.R        per-sample QC + scDblFinder
AOMDSS_GSE264607/01_integration_harmony.R    normalisation, PCA and Harmony integration
AOMDSS_GSE264607/02_lineage_identification.R T/NK/ILC compartment from marker module scores
AOMDSS_GSE264607/03_subcluster_rpca.R        RPCA reintegration and UMAP of that compartment
AOMDSS_GSE264607/04_annotate_populations.R   population annotation with absolute thresholds
```

Input: `GSE264607_RAW.tar` from GEO, extracted into a folder with one subfolder per
sample (CellRanger format). The path is set in `RAW_DIR`, at the top of script `00`.

Output: `data/GSE264607_Tcells_final_con_labels.rds` → `AOM_RDS`.

This is the longest pipeline (several hours; step `00` needs a fair amount of memory).

## Ulcerative colitis — GSE125527 (Boland et al.)

```
UC_GSE125527/00_build_seurat_object.R   downloads from GEO and builds the object
UC_GSE125527/02_tcell_umap.R            T cells, normalisation, PCA and UMAP
UC_GSE125527/03_tcell_tissue_only.R     keeps the tissue cells only
```

Script `00` downloads the files from the GEO FTP into `GSE125527_data/` by itself.

Output: `GSE125527_data/GSE125527_Tcells_tissue.rds` → `COLITIS_RDS`.

## mLN proteomics

No object to build: the Excel file with the per-sample copy numbers is read directly
(`PROT_XLSX` in `R/00_paths.R`, sheet "Complete analysis").
