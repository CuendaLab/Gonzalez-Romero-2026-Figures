###############################################################
## Step required before any CD4/CD8 analysis of GSE125527:
## unlike the CRC object (which already carries ClusterMidway = TCD4/TCD8
## assigned by the original atlas), the T1-T17 subclusters of
## GSE125527_Tcell_cluster.csv do not carry that label in the
## supplementary file. The paper that deposited this GEO series does
## annotate them explicitly:
##
##   Boland et al. 2020, Science Immunology
##   "Heterogeneity and clonal relationships of adaptive immune cells
##    in ulcerative colitis revealed by single-cell analyses"
##   Fig. 4A, tabla "Phenotypic annotations of T lymphocyte clusters"
##   (read directly from the manuscript PDF, not derived here):
##
##     Naive/Memory : T3, T4, T6, T9, T11, T12, T17
##     gd           : T8, T13, T16
##     MAIT         : T15
##     Tfh          : T5
##     Treg         : T7
##     CD8+ TRM     : T1, T2, T10, T14
##
## (OJO: NO es Smillie et al. 2019 Cell -- ese es un dataset distinto
## de colitis, GSE125527 es de este paper.)
##
## El paper NO separa el grupo "Naive/Memory" en CD4 vs CD8 -- eso lo
## hacemos aca por marcadores (CD8A/CD8B %positivas), solo dentro de
## ese subgrupo. Para TRM/Treg/Tfh/gd/MAIT usamos la anotacion oficial
## del paper, no el umbral por marcadores (el umbral por si solo
## clasificaba mal T1/T2 como CD4 y T8 como CD8 -- ver mas abajo).
###############################################################

library(Seurat)
library(dplyr)

INPUT_RDS <- file.path("GSE125527_data", "GSE125527_seurat_raw.rds")
OUTPUT_DIR <- file.path("output", "cd4_cd8_classification")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

seu <- readRDS(INPUT_RDS)

## Solo celulas anotadas como "T" (incluye T convencionales + NKT-like;
## las NK puras ya estan aparte en celltype == "NK")
Tcells <- subset(seu, subset = celltype == "T")
cat("Celulas T:", ncol(Tcells), "en", length(unique(Tcells$Tcell_cluster)), "subclusters\n")

## El objeto raw trae `counts` = UMI crudos (no como en el proyecto CRC,
## donde `counts` ya eran valores normalizados del portal SCP1162).
## Aca SI corresponde el pipeline estandar de Seurat.
Tcells <- NormalizeData(Tcells)

markers <- c("CD3D", "CD3E",                    # pan-T (control de calidad del subset)
             "CD4", "CD8A", "CD8B",             # linaje convencional
             "FOXP3", "IL2RA",                  # Treg (se deja dentro de CD4)
             "IL7R", "GZMK",                     # memoria / citotoxico, contexto
             "TRAV1-2", "KLRB1", "ZBTB16", "DPP4",  # MAIT
             "TRDC", "TRGC1", "TRGC2")          # gamma-delta

markers <- intersect(markers, rownames(Tcells))

## --- % de celulas positivas y media por subcluster -----------------
## CD4 se captura mal en 10x 3' (transcrito corto/poco expresado): su
## media suele ser baja incluso en celulas CD4+ de verdad. Por eso NO
## clasificamos por media de CD4, sino por %positivas de CD8A/CD8B
## (marcador robusto): CD8 si el subcluster es claramente CD8A/B+;
## el resto queda como CD4 "por defecto" (CD8-negativo). Los Tregs
## (T7, FOXP3+) se dejan dentro de CD4 -- no se separan como grupo
## aparte para este primer analisis.
Idents(Tcells) <- Tcells$Tcell_cluster
Idents(Tcells) <- factor(Idents(Tcells), levels = paste0("T", 1:17))

dp <- DotPlot(Tcells, features = markers)
dp_data <- dp$data  # columnas: avg.exp, pct.exp, features.plot, id, avg.exp.scaled

pct <- reshape(dp_data[, c("id", "features.plot", "pct.exp")],
               idvar = "id", timevar = "features.plot", direction = "wide")
names(pct) <- sub("^pct\\.exp\\.", "", names(pct))
names(pct)[1] <- "cluster"
pct$cluster <- as.character(pct$cluster)  # evitar bug de indexar por codigo de factor
pct <- pct[order(as.integer(sub("^T", "", pct$cluster))), ]

avgexp <- reshape(dp_data[, c("id", "features.plot", "avg.exp")],
                   idvar = "id", timevar = "features.plot", direction = "wide")
names(avgexp) <- sub("^avg\\.exp\\.", "", names(avgexp))
names(avgexp)[1] <- "cluster"
avgexp$cluster <- as.character(avgexp$cluster)
avgexp <- avgexp[match(pct$cluster, avgexp$cluster), ]

## --- Anotacion oficial del paper (Boland et al. 2020, Fig. 4A) ------
paper_annotation <- c(
  T1 = "CD8_TRM", T2 = "CD8_TRM", T10 = "CD8_TRM", T14 = "CD8_TRM",
  T7 = "Treg", T5 = "Tfh",
  T8 = "gd", T13 = "gd", T16 = "gd",
  T15 = "MAIT",
  T3 = "Naive_Memory", T4 = "Naive_Memory", T6 = "Naive_Memory",
  T9 = "Naive_Memory", T11 = "Naive_Memory", T12 = "Naive_Memory",
  T17 = "Naive_Memory"
)

## CD4/CD8A-CD8B por marcadores, solo para partir el grupo Naive/Memory
## (el paper no lo hace); para el resto se usa la anotacion oficial.
marker_call <- ifelse(pmax(pct$CD8A, pct$CD8B) > 20, "CD8", "CD4")
names(marker_call) <- pct$cluster

avg <- pct
avg$paper_annotation <- paper_annotation[avg$cluster]
avg$call <- dplyr::case_when(
  avg$paper_annotation == "CD8_TRM"     ~ "CD8",
  avg$paper_annotation == "Treg"        ~ "CD4",   # Treg dentro de CD4, a pedido
  avg$paper_annotation == "Tfh"         ~ "CD4",
  avg$paper_annotation == "gd"          ~ "gdT",
  avg$paper_annotation == "MAIT"        ~ "MAIT",
  avg$paper_annotation == "Naive_Memory" ~ marker_call[avg$cluster]
)

cat("\n--- % of positive cells per T subcluster (marker > 0) ---\n")
print_df <- avg[, c("cluster", "paper_annotation", "call", markers)]
print_df[markers] <- round(print_df[markers], 1)
print(print_df, row.names = FALSE)

write.csv(avg, file.path(OUTPUT_DIR, "Tcell_cluster_marker_pct_positive.csv"), row.names = FALSE)
write.csv(avgexp, file.path(OUTPUT_DIR, "Tcell_cluster_marker_avgexp.csv"), row.names = FALSE)

## --- DotPlot: the standard way to see this at a glance ---
p <- dp +
  ggplot2::coord_flip() +
  ggplot2::ggtitle("CD4/CD8/MAIT/gd markers per T subcluster (GSE125527)")

ggplot2::ggsave(file.path(OUTPUT_DIR, "Tcell_cluster_CD4_CD8_dotplot.png"),
                 p, width = 7, height = 6, dpi = 300)
ggplot2::ggsave(file.path(OUTPUT_DIR, "Tcell_cluster_CD4_CD8_dotplot.pdf"),
                 p, width = 7, height = 6)

cat("\nSaved:\n -", file.path(OUTPUT_DIR, "Tcell_cluster_marker_means.csv"), "\n -",
    file.path(OUTPUT_DIR, "Tcell_cluster_CD4_CD8_dotplot.png"), "\n")

## --- Composicion de cada subcluster por condicion (sano/enfermo) ---
comp <- table(Tcells$Tcell_cluster, Tcells$disease_assignment)
cat("\n--- n celulas por subcluster x condicion ---\n")
print(comp)
write.csv(as.data.frame.matrix(comp),
          file.path(OUTPUT_DIR, "Tcell_cluster_by_condition.csv"))
