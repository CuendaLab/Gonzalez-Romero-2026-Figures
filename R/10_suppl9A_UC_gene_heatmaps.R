###############################################################
## Gene-by-gene heatmaps of the 9 processes in CD4 T cells of ulcerative colitis (GSE125527)
##
## Figure  : Suppl. Fig. 9A
## Input   : COLITIS_RDS, SIG_HUMAN
## Output  : figures/Suppl_Figure_9/SupplFigure9_A_GSE125527_CD4_gene_heatmaps/<proceso>_CD4.{png,pdf}
###############################################################

source("R/00_paths.R")

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

###############################################################
## CONFIGURATION
###############################################################

INPUT_RDS  <- COLITIS_RDS
GENE_CSV   <- SIG_HUMAN
OUTPUT_DIR <- file.path(DIR_S9, "SupplFigure9_A_GSE125527_CD4_gene_heatmaps")

MIN_CELLS <- 15

CELL_H_MM <- 4
CELL_W_MM <- 6

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (f in c(INPUT_RDS, GENE_CSV)) {
  if (!file.exists(f)) stop("Required input not found: ", f)
}

###############################################################
## LOAD DATA
###############################################################

Tissue <- readRDS(INPUT_RDS)

Tissue$lineage <- dplyr::case_when(
  Tissue$Tcell_cluster %in% c("T1", "T2", "T10", "T14") ~ "CD8",
  Tissue$Tcell_cluster %in% c("T5", "T7")                ~ "CD4",
  TRUE ~ NA_character_
)
Tissue <- subset(Tissue, subset = !is.na(lineage))
Tissue$lineage <- factor(Tissue$lineage, levels = c("CD4", "CD8"))

cat("Pacientes:", length(unique(Tissue$patient_assignment)),
    "(", sum(!duplicated(Tissue$patient_assignment[Tissue$disease_assignment == "healthy"])), "sanos,",
    sum(!duplicated(Tissue$patient_assignment[Tissue$disease_assignment == "diseased"])), "UC )\n")

###############################################################
## GENE -> SIGNATURE MAP
###############################################################

gene_level <- read.csv(GENE_CSV, stringsAsFactors = FALSE)

signames <- c(
  "Glycolysis", "Transporters", "Methyltransferase", "Cytokine",
  "TCA", "Trafficking", "OXPHOS", "Tcell_activation", "Translation"
)

sig_labels <- c(
  Glycolysis = "Glycolysis", Transporters = "Transporters",
  Methyltransferase = "Methyltransferase", Cytokine = "Cytokine Regulation",
  TCA = "TCA cycle", Trafficking = "Trafficking", OXPHOS = "OXPHOS",
  Tcell_activation = "T cell activation", Translation = "Translation"
)

sig2genes <- setNames(vector("list", length(signames)), signames)
for (i in seq_len(nrow(gene_level))) {
  labs <- trimws(strsplit(gene_level$signatures[i], ";")[[1]])
  for (s in labs) {
    if (nzchar(s) && s %in% signames) {
      sig2genes[[s]] <- c(sig2genes[[s]], gene_level$gene[i])
    }
  }
}
sig2genes <- lapply(sig2genes, function(g) intersect(g, rownames(Tissue)))
cat("Genes de la firma encontrados en el panel de colitis (de", nrow(gene_level), "totales):\n")
print(sapply(sig2genes, length))

###############################################################
## CD4 / CD8 SUBSETS
###############################################################

CD4 <- subset(Tissue, subset = lineage == "CD4")
CD8 <- subset(Tissue, subset = lineage == "CD8")

###############################################################
## PER-PATIENT, Z-SCORED, CLIPPED EXPRESSION MATRIX
###############################################################

get_matrix <- function(obj, genes) {
  if (length(genes) == 0) return(NULL)

  obj$group_patient <- paste(obj$disease_assignment, obj$patient_assignment, sep = "-")

  n_per_group <- table(obj$group_patient)
  keep <- names(n_per_group)[n_per_group >= MIN_CELLS]
  dropped <- names(n_per_group)[n_per_group < MIN_CELLS]
  if (length(dropped) > 0) {
    cat("  (excluded, fewer than", MIN_CELLS, "cells:", paste(dropped, collapse = ", "), ")\n")
  }
  obj <- subset(obj, subset = group_patient %in% keep)

  avg <- AverageExpression(
    obj, group.by = "group_patient", features = genes, assays = "RNA", layer = "data"
  )$RNA
  if (is.null(avg)) return(NULL)

  mat_raw <- as.matrix(avg)
  genes_present <- genes[genes %in% rownames(mat_raw)]
  if (length(genes_present) == 0) return(NULL)
  mat_raw <- mat_raw[genes_present, , drop = FALSE]

  cols <- colnames(mat_raw)
  mat_raw <- mat_raw[, c(cols[grepl("^healthy", cols)], cols[grepl("^diseased", cols)]), drop = FALSE]
  is_healthy <- grepl("^healthy", colnames(mat_raw))

  row_mean <- rowMeans(mat_raw)
  row_sd   <- apply(mat_raw, 1, sd)

  z <- (mat_raw - row_mean) / row_sd
  z[is.na(z)] <- 0
  z[z > 2]  <- 2
  z[z < -2] <- -2

  H_raw <- rowMeans(mat_raw[, is_healthy, drop = FALSE])
  D_raw <- rowMeans(mat_raw[, !is_healthy, drop = FALSE])
  summary_z <- cbind(H = (H_raw - row_mean) / row_sd, D = (D_raw - row_mean) / row_sd)
  summary_z[is.na(summary_z)] <- 0
  summary_z[summary_z > 2]  <- 2
  summary_z[summary_z < -2] <- -2

  ## --- SORT BY TREND: from highest to lowest Z in the disease column
  ## (mas subido en Diseased arriba, mas bajado abajo) -- unica
  ## diferencia real respecto a 04a.
  ord <- order(-summary_z[, "D"])
  z <- z[ord, , drop = FALSE]
  summary_z <- summary_z[ord, , drop = FALSE]

  list(patient = z, summary = summary_z)
}

###############################################################
## HEATMAP: Healthy/Disease summary strip + one column per patient
###############################################################

col_fun <- colorRamp2(seq(-2, 2, length.out = 100),
                       colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100))

make_ht <- function(res, titulo) {
  if (is.null(res)) return(NULL)
  mat <- res$patient
  summary_mat <- res$summary

  grupos <- factor(
    ifelse(grepl("^healthy", colnames(mat)), "Healthy", "Diseased"),
    levels = c("Healthy", "Diseased")
  )

  ht_summary <- Heatmap(
    summary_mat, name = "HD_summary", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    rect_gp = gpar(col = NA),
    border = "black", border_gp = gpar(lwd = 2),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    width = unit(CELL_W_MM * 2, "mm"),
    height = unit(CELL_H_MM * nrow(mat), "mm")
  )

  ht_main <- Heatmap(
    mat, name = "Z-score", col = col_fun,
    column_split = grupos, column_gap = unit(3, "mm"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = NA),
    border = "black", border_gp = gpar(lwd = 2),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 45,
    width = unit(CELL_W_MM * ncol(mat), "mm"),
    height = unit(CELL_H_MM * nrow(mat), "mm"),
    column_title = titulo
  )

  ht_summary + ht_main
}

###############################################################
## GENERATE ONE HEATMAP PER PROCESS x CELL TYPE
###############################################################

for (name in names(sig2genes)) {

  genes <- sig2genes[[name]]
  label <- sig_labels[name]

  cat("\n==", label, "==\n")
  mat_cd4 <- get_matrix(CD4, genes)
  mat_cd8 <- get_matrix(CD8, genes)

  if (is.null(mat_cd4) && is.null(mat_cd8)) {
    cat("Skipping", name, "(no genes found)\n")
    next
  }

  max_genes <- max(nrow(mat_cd4$patient) %||% 0, nrow(mat_cd8$patient) %||% 0)
  max_cols  <- max(ncol(mat_cd4$patient) %||% 0, ncol(mat_cd8$patient) %||% 0)
  pdf_height <- max(6, (max_genes * CELL_H_MM) / 25.4 + 2)
  pdf_width  <- max(10, max_cols * 0.25) + 1

  for (ct in c("CD4", "CD8")) {
    mat <- if (ct == "CD4") mat_cd4 else mat_cd8
    if (is.null(mat)) next
    ht <- make_ht(mat, paste(ct, label, "por paciente, ordenado por tendencia (Sano vs Enfermo)"))

    pdf(file.path(OUTPUT_DIR, paste0(name, "_", ct, ".pdf")), width = pdf_width, height = pdf_height)
    draw(ht, column_title = label, column_title_gp = gpar(fontsize = 16, fontface = "bold"))
    decorate_heatmap_body("HD_summary", {
      grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 0.8, col = "black"))
    })
    dev.off()

    png(file.path(OUTPUT_DIR, paste0(name, "_", ct, ".png")),
        width = pdf_width, height = pdf_height, units = "in", res = 300)
    draw(ht, column_title = label, column_title_gp = gpar(fontsize = 16, fontface = "bold"))
    decorate_heatmap_body("HD_summary", {
      grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 0.8, col = "black"))
    })
    dev.off()
  }
}

cat("\nDone. Heatmaps saved to", OUTPUT_DIR, "\n")
print(list.files(OUTPUT_DIR))
