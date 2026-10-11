###############################################################
## Dataset names as they appear in the figures (single source)
##
## The mLN proteomics carries a different name in each block because the
## comparison is different (colitis vs day 0 and cancer vs day 0).
## 
###############################################################

DS_PRO_COL <- "Colitis mLN proteomic\n(mouse)"   # mLN proteomics, Colitis vs day 0
DS_UC      <- "UC colon scRNA\n(human)"          # GSE125527 (Boland), UC vs healthy
DS_PRO_CAC <- "CAC mLN proteomic\n(mouse)"       # mLN proteomics, Cancer vs day 0
DS_AOM     <- "CAC colon scRNA\n(mouse)"         # AOM-DSS GSE264607, Healthy -> CAC
DS_CRC     <- "CRC colon scRNA\n(human)"         # Pelka et al. (SCP1162), Tumour vs Normal

## column order: Colitis block [1-2] | Cancer block [3-5]
DATASETS <- c(DS_PRO_COL, DS_UC, DS_PRO_CAC, DS_AOM, DS_CRC)

## three-line version, for narrow columns (single-lineage variants)
NARROW_LABELS <- setNames(c("Colitis mLN\nproteomic\n(mouse)", "UC colon\nscRNA\n(human)",
                            "CAC mLN\nproteomic\n(mouse)", "CAC colon\nscRNA\n(mouse)",
                            "CRC colon\nscRNA\n(human)"), DATASETS)

## description of each comparison for the figure footnote; key = first line of the name
CAPTION_LINES <- setNames(
  c("Colitis mLN proteomic, Colitis vs day 0 (copy numbers/cell, Welch t, n = 3 vs 3)",
    "UC colon scRNA, UC vs healthy (Wilcoxon)",
    "CAC mLN proteomic, Cancer vs day 0 (copy numbers/cell, Welch t, n = 3 vs 3)",
    "CAC colon scRNA, Cancer vs Healthy (Wilcoxon)",
    "CRC colon scRNA, Cancer vs Healthy (paired Wilcoxon, n = 36)"),
  sub("\n.*", "", DATASETS))
