###############################################################
## Per-patient module scores in CRC (CD4 and CD8)
##
## Figure  : Fig. 8C (CD4) and Suppl. Fig. 8B (CD8)
## Input   : CRC_RDS, SIG_HUMAN
## Output  : Figure8_C_CRC_CD4_per_patient_scores / SupplFigure8_B_CRC_CD8_per_patient_scores (+ tables)
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(Seurat); library(Matrix); library(dplyr); library(tidyr); library(ggplot2) })

INPUT_RDS <- CRC_RDS
GENE_CSV  <- SIG_HUMAN
OUT_DIR <- DIR_TAB

PANEL_E <- MODULES
COL_CD4 <- "#376199"   # same CD4/CD8 colours as Fig. 8C
COL_CD8 <- "#B37CB8"

## ---- module scores per cell and per patient ----
Tcells <- readRDS(INPUT_RDS)
Tcells_cd <- subset(Tcells, subset = ClusterMidway %in% c("TCD4", "TCD8"))
rm(Tcells); invisible(gc())

signames <- c("Glycolysis", "Ribosome", "Transporters", "Methyltransferase", "Cytokine",
              "TCA", "Trafficking", "OXPHOS", "Tcell_activation", "Translation")
gene_level <- read.csv(GENE_CSV, stringsAsFactors = FALSE)
sig2genes <- setNames(vector("list", length(signames)), signames)
for (i in seq_len(nrow(gene_level))) {
  for (s in trimws(strsplit(gene_level$signatures[i], ";")[[1]])) {
    if (nzchar(s) && s %in% signames) sig2genes[[s]] <- c(sig2genes[[s]], gene_level$gene[i])
  }
}
sig2genes <- lapply(sig2genes, function(g) intersect(g, rownames(Tcells_cd)))

all_genes <- unique(unlist(sig2genes))
expr <- as.matrix(GetAssayData(Tcells_cd, assay = "RNA", layer = "counts")[all_genes, , drop = FALSE])
Z <- (expr - rowMeans(expr)) / (apply(expr, 1, sd) + 1e-9)

S <- as.data.frame(sapply(signames, function(s) colMeans(Z[sig2genes[[s]], , drop = FALSE])))
S$donor <- Tcells_cd$orig.ident; S$mid <- Tcells_cd$ClusterMidway; S$spec <- Tcells_cd$Grupo

pb <- S |> group_by(donor, mid, spec) |>
  summarise(across(all_of(signames), \(x) mean(x, na.rm = TRUE)), n_cells = n(), .groups = "drop")

## ---- per-patient table (all modules, CD4 and CD8) ----
pb_long <- pb |>
  mutate(lineage = sub("^T", "", mid)) |>
  pivot_longer(all_of(signames), names_to = "sig", values_to = "score") |>
  select(donor, lineage, spec, n_cells, sig, score)
write.csv(pb_long, file.path(DIR_TAB, "CRC_module_scores_per_patient.csv"), row.names = FALSE)

## ---- paired test (only donors with both specimens), BH across the 9 modules ----
paired <- pb_long |> filter(sig %in% MODULES) |>
  select(-n_cells) |>
  pivot_wider(names_from = spec, values_from = score) |>
  filter(!is.na(Normal), !is.na(Tumor))

stats <- paired |> group_by(lineage, sig) |>
  summarise(n = n(), delta_mean = mean(Tumor - Normal),
            p = suppressWarnings(wilcox.test(Tumor, Normal, paired = TRUE)$p.value), .groups = "drop") |>
  group_by(lineage) |> mutate(fdr = p.adjust(p, "BH")) |> ungroup()
write.csv(stats, file.path(DIR_TAB, "CRC_paired_stats_9modules.csv"), row.names = FALSE)

## ---- figure: one panel per process, CD4 and CD8 separately ----
fmt_fdr <- function(f) ifelse(f < 0.001, sprintf("FDR = %.1e", f), sprintf("FDR = %.3f", f))
LIN_COLS <- c(CD4 = COL_CD4, CD8 = COL_CD8)

for (lin in c("CD4", "CD8")) {
  dE <- paired |> filter(lineage == lin, sig %in% PANEL_E) |>
    pivot_longer(c(Normal, Tumor), names_to = "spec", values_to = "score") |>
    ## Healthy/Cancer labels in the figure; the metadata keeps Normal/Tumor
    mutate(spec = factor(spec, levels = c("Normal", "Tumor"), labels = c("Healthy", "Cancer")),
           panel = factor(MOD_LABELS[sig], levels = MOD_LABELS[PANEL_E]))
  lab <- stats |> filter(lineage == lin, sig %in% PANEL_E) |>
    mutate(panel = factor(MOD_LABELS[sig], levels = MOD_LABELS[PANEL_E]),
           label = paste0(fmt_fdr(fdr), "\nn = ", n, " pairs"))
  mn <- dE |> group_by(panel, spec) |> summarise(score = mean(score), .groups = "drop")

  ## figure builder; fonts and labels are passed in so the same code serves
  ## both lineages
  make_pE <- function(lab_df, base_pt, strip_pt, x_pt, fdr_pt, wrap_n = 30) {
    ggplot(dE, aes(spec, score)) +
      geom_line(aes(group = donor), colour = "grey75", linewidth = 0.3) +
      geom_point(aes(fill = spec), shape = 21, size = 1.8, colour = "white", stroke = 0.3) +
      geom_line(data = mn, aes(group = 1), colour = "black", linewidth = 0.8) +
      geom_point(data = mn, shape = 21, size = 3, fill = "white", colour = "black", stroke = 0.8) +
      geom_text(data = lab_df, aes(x = 1.5, y = Inf, label = label), vjust = 1.2, size = fdr_pt,
                colour = "grey20", lineheight = 0.95) +
      facet_wrap(~ panel, nrow = 3, scales = "free_y", labeller = as_labeller(function(x) lab_wrap(x, wrap_n))) +
      scale_fill_manual(values = c(Healthy = "grey60", Cancer = LIN_COLS[[lin]]), guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
      labs(x = NULL, y = "Module score (mean Z per patient)",
           title = sprintf("%s T cells, CRC (Pelka): per-patient module scores", lin),
           caption = "Lines join the same patient; black = mean. Paired Wilcoxon, BH across the 9 modules.") +
      theme_classic(base_size = base_pt) +
      theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = strip_pt),
            axis.text.x = element_text(colour = "grey15", size = x_pt),
            plot.title.position = "plot",
            plot.caption = element_text(size = 7, colour = "grey35", hjust = 0),
            plot.caption.position = "plot")
  }
  ## version used in the figures: large labels, FDR only (no "n = 36 pairs")
  lab_big <- lab |> mutate(label = fmt_fdr(fdr))
  pE_big  <- make_pE(lab_big, 12, 15, 15, 4, wrap_n = 20)
  out_dir <- if (lin == "CD4") DIR_F8 else DIR_S8
  stub    <- if (lin == "CD4") "Figure8_C_CRC_CD4_per_patient_scores"
             else              "SupplFigure8_B_CRC_CD8_per_patient_scores"
  ggsave(file.path(out_dir, paste0(stub, ".pdf")), pE_big, width = 8.5, height = 9.5, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(stub, ".png")), pE_big, width = 8.5, height = 9.5, dpi = 200, bg = "white")

  cat("\n", lin, "\n"); print(lab |> select(panel, n, delta_mean, p, fdr))
}
cat("OK\n")
