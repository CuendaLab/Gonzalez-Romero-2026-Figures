###############################################################
## Cross-dataset summary table (change and FDR per module, dataset and lineage)
##
## Figure  : input for Suppl. Fig. 9B
## Input   : input_tables/*.csv and tables/Table_proteomics_mLN_module_delta_vs_Basal.csv
## Output  : tables/panelD_table.csv
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")
source("R/00_labels_datasets.R")

suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

OUT_DIR <- DIR_TAB
CRC_CSV <- file.path("input_tables", "CRC_Pelka_module_delta.csv")      # per-module change, CRC
AOM_CSV <- file.path("input_tables", "AOMDSS_module_delta.csv")          # per-module change, AOM/DSS
UC_CSV  <- file.path("input_tables", "GSE125527_module_delta.csv")       # per-module change, ulcerative colitis
PRO_CSV <- file.path(DIR_TAB, "Table_proteomics_mLN_module_delta_vs_Basal.csv")  # script 11a

COL_UP <- "#C23B22"; COL_DOWN <- "#2C7BB6"

std <- function(d, block, dataset, delta, p) {
  d |> transmute(block = block, dataset = dataset, lineage = comp, sig, delta = {{ delta }}, p = {{ p }}) |>
    filter(sig %in% MODULES)
}

## Proteomics: "absolute" version; the copy numbers are already normalised to
## the histone content per cell.
build <- function(prot_version = "absolute") {
  aom <- read.csv(AOM_CSV); pro <- read.csv(PRO_CSV) |> filter(version == prot_version)
  bind_rows(
    ## the AOM/DSS Healthy -> Colitis comparison is not shown
    std(filter(pro, comparison == "Colitis vs Basal"),   "Colitis", DS_PRO_COL, delta_mean, p),
    std(read.csv(UC_CSV),                                "Colitis", DS_UC,      delta_mean, p),
    std(filter(pro, comparison == "Cancer vs Basal"),    "Cancer",  DS_PRO_CAC, delta_mean, p),
    std(filter(aom, comparison == "Healthy -> Cancer"),  "Cancer",  DS_AOM,     Delta, Pvalue),
    std(read.csv(CRC_CSV),                               "Cancer",  DS_CRC,     delta_mean, p)) |>
    group_by(block, dataset, lineage) |>
    mutate(fdr = p.adjust(p, "BH"), scaled = delta / max(abs(delta))) |> ungroup() |>
    mutate(sig_lab = factor(MOD_LABELS[match(sig, MODULES)], levels = rev(MOD_LABELS)),
           dataset = factor(dataset, levels = DATASETS),
           block = factor(block, levels = c("Colitis", "Cancer")),
           filled = fdr < 0.05,
           sig_status = factor(ifelse(filled, "FDR < 0.05", "FDR ≥ 0.05"), levels = c("FDR < 0.05", "FDR ≥ 0.05")))
}

plot_block <- function(d, title, show_y, border_legend = TRUE) {
  ggplot(d, aes(lineage, sig_lab)) +
    ## very light grey outline on the non-significant dots so that those with change ~ 0 are visible
    ## two layers with fixed stroke; only one block draws the border legend
    ## (the Colitis block has no significant dots and its legend would differ)
    geom_point(data = filter(d, !filled), aes(fill = scaled, colour = sig_status), shape = 21, size = 5.5, stroke = 0.4) +
    geom_point(data = filter(d, filled),  aes(fill = scaled, colour = sig_status), shape = 21, size = 5.5, stroke = 1.1) +
    facet_grid(. ~ dataset, switch = "x") +
    scale_y_discrete(labels = lab_wrap) +
    scale_fill_gradient2(low = COL_DOWN, mid = "#F7F7F7", high = COL_UP, midpoint = 0,
                         limits = c(-1, 1), name = "Δ relative to\nlargest |Δ|\nin its column") +
    scale_colour_manual(values = c("FDR < 0.05" = "black", "FDR ≥ 0.05" = "grey80"), name = NULL, drop = FALSE) +
    guides(fill = guide_colourbar(order = 1, barheight = unit(25, "mm")),
           colour = if (border_legend) guide_legend(order = 2, override.aes = list(size = 5, fill = "white", stroke = c(1.1, 0.4))) else "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = 8, lineheight = 0.9),
          panel.spacing.x = unit(3, "mm"),
          axis.text.y = if (show_y) element_text(colour = "grey15") else element_blank(),
          axis.text.x = element_text(colour = "grey15", face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}

d <- build("absolute")
write.csv(d, file.path(OUT_DIR, "panelD_table.csv"), row.names = FALSE)

cap <- paste0(
  "Colitis: ", paste(CAPTION_LINES[sub("\n.*", "", c(DS_PRO_COL, DS_UC))], collapse = "; "), ".\n",
  "Cancer: ", paste(CAPTION_LINES[sub("\n.*", "", c(DS_PRO_CAC, DS_AOM, DS_CRC))], collapse = "; "), ".\n",
  "BH across the 9 modules per column. Colour: Δ / max|Δ| within each column; intensities are NOT comparable between columns. Black circle: FDR < 0.05.")

p <- (plot_block(filter(d, block == "Colitis"), "Colitis", TRUE, border_legend = FALSE) |
      plot_block(filter(d, block == "Cancer"),  "Cancer",  FALSE)) +
  plot_layout(guides = "collect", widths = c(2, 3)) +
  plot_annotation(title = "Disease vs control, per module", caption = cap,
                  theme = theme(plot.title = element_text(size = 13),
                                plot.caption = element_text(size = 7, colour = "grey35", hjust = 0, lineheight = 1.1))) &
  theme(legend.justification = "top", legend.box.just = "left")

print(d |> group_by(block, dataset, lineage) |>
        summarise(n_up = sum(delta > 0), n_sig = sum(filled), min_fdr = min(fdr), .groups = "drop") |>
        mutate(dataset = gsub("\n", " ", dataset)), n = Inf)
