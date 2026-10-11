###############################################################
## Paths and output folders
##
## The ONLY place to edit when running this code on another machine:
## the four data paths below. Everything else is relative to the root
## of this repository (scripts are run from there).
###############################################################

## --- input data (see README for where each one is downloaded from) ----------
CRC_DIR     <- "/Users/ester/Desktop/SCP1162 CRC Tcell Analysis _claude code"   # Pelka et al. 2021, SCP1162
AOM_DIR     <- "/Users/ester/Desktop/AOM-DSS_GSE264607_FINAL"                   # GSE264607, AOM/DSS
COLITIS_DIR <- "/Users/ester/Desktop/Colitis GSE1255227_claude code"            # GSE125527, ulcerative colitis
PROT_XLSX   <- "/Users/ester/Desktop/proteómica Diego /Diegos adj p value analysis_LS_DGR.xlsx"  # mLN proteomics

## Seurat objects and gene signatures used by the scripts
CRC_RDS     <- file.path(CRC_DIR,     "data", "Tcells_final_con_UMAP_y_labels.rds")
AOM_RDS     <- file.path(AOM_DIR,     "data", "GSE264607_Tcells_final_con_labels.rds")
COLITIS_RDS <- file.path(COLITIS_DIR, "GSE125527_data", "GSE125527_Tcells_tissue.rds")
SIG_HUMAN   <- file.path("gene_signatures", "cd4_gene_level_prognosis.csv")        # 412 genes, 9 processes
SIG_MOUSE   <- file.path("gene_signatures", "cd4_gene_level_prognosis_mouse.csv")  # 409 mouse orthologs

## --- output folders, one per figure -----------------------------------------
DIR_F8 <- file.path("figures", "Figure_8")
DIR_S7 <- file.path("figures", "Suppl_Figure_7")
DIR_S8 <- file.path("figures", "Suppl_Figure_8")
DIR_S9 <- file.path("figures", "Suppl_Figure_9")
DIR_TAB <- "tables"      # supporting tables (per-patient scores, statistics, genes)

for (d in c(DIR_F8, DIR_S7, DIR_S8, DIR_S9, DIR_TAB))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
