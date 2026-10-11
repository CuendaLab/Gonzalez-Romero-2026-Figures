###############################################################
## Cross-dataset dot matrix, CD4 and CD8 in a single figure
##
## Figure  : Suppl. Fig. 9B
## Input   : tables/panelD_table.csv (script 11b)
## Output  : figures/Suppl_Figure_9/SupplFigure9_B_crossdataset_dotmatrix.{png,pdf}
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")
source("R/00_labels_datasets.R")

suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

OUT_DIR <- DIR_S9
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_UP <- "#C23B22"; COL_DOWN <- "#2C7BB6"
SIG_LEVELS <- c("FDR < 0.05", "FDR ≥ 0.05")
SPACING <- 0.36   # inches between columns
PANEL_H <- 3.2    # height of the dot area, inches (9 rows)

d_all <- read.csv(file.path(DIR_TAB, "panelD_table.csv")) |>
  mutate(sig_lab = factor(sig_lab, levels = rev(MOD_LABELS)),
         dataset = factor(dataset, levels = DATASETS),
         block = factor(block, levels = c("Colitis", "Cancer")),
         sig_status = factor(sig_status, levels = SIG_LEVELS))

make_caption <- function(d, width_in, show_sig = TRUE) {
  keys <- d |> distinct(block, dataset) |> arrange(block, dataset) |>
    mutate(key = sub("\n.*", "", dataset))     # the names are already unique within each block
  lines <- sapply(levels(droplevels(keys$block)), function(b)
    paste0(b, ": ", paste(CAPTION_LINES[keys$key[keys$block == b]], collapse = "; "), "."))
  lines <- c(lines, paste0("BH across the 9 modules per column. Colour: Δ / max|Δ| within each column; intensities are NOT comparable between columns.",
                           if (show_sig) " Black circle: FDR < 0.05." else ""))
  paste(unlist(lapply(lines, strwrap, width = floor(width_in * 17))), collapse = "\n")   # ~17 characters per inch at 7 pt
}

plot_block <- function(d, title, show_y, single = FALSE, show_sig = TRUE) {
  ## the border legend comes from an invisible helper layer that always has both
  ## categories, so it is identical in every block (patchwork collects them into one)
  ## even if a block has no significant dots
  key <- d[c(1, 1), ] |> mutate(sig_status = factor(SIG_LEVELS, levels = SIG_LEVELS))
  ## single lineage: datasets become columns of one panel (names rotated),
  ## so that the column spacing matches the CD4+CD8 layout
  xvar <- if (single) "dataset" else "lineage"
  ggplot(d, aes(.data[[xvar]], sig_lab)) +
    (if (show_sig) geom_point(data = key, aes(colour = sig_status), alpha = 0, size = 5.5)) +
    geom_point(aes(fill = scaled), colour = "grey80", shape = 21, size = 5.5, stroke = 0.4) +
    (if (show_sig) geom_point(data = filter(d, filled), aes(fill = scaled), colour = "black", shape = 21,
                              size = 5.5, stroke = 1.1)) +
    (if (single) scale_x_discrete(labels = function(x) sub("\n", " ", x))
     else facet_grid(. ~ dataset, switch = "x", labeller = as_labeller(NARROW_LABELS))) +
    scale_y_discrete(limits = rev(MOD_LABELS), labels = lab_wrap) +   # fixed module order (the helper layer would otherwise change it)
    scale_fill_gradient2(low = COL_DOWN, mid = "#F7F7F7", high = COL_UP, midpoint = 0,
                         limits = c(-1, 1), name = "Δ relative to\nlargest |Δ|\nin its column") +
    scale_colour_manual(values = c("FDR < 0.05" = "black", "FDR ≥ 0.05" = "grey80"),
                        limits = SIG_LEVELS, drop = FALSE, name = NULL) +
    guides(fill = guide_colourbar(order = 1, barheight = unit(25, "mm")),
           colour = if (show_sig) guide_legend(order = 2, override.aes = list(size = 5, shape = 21, fill = "white", alpha = 1, stroke = c(1.1, 0.4))) else "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(size = 7, lineheight = 0.9),
          panel.spacing.x = unit(0.5, "mm"),
          axis.text.y = if (show_y) element_text(colour = "grey15") else element_blank(),
          axis.text.x = if (single) element_text(colour = "grey15", size = 7, angle = 90, hjust = 1, vjust = 0.5)
                        else element_text(colour = "grey15", face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}

make_variant <- function(lineages, blocks, name, title, show_sig = TRUE) {
  d <- d_all |> filter(lineage %in% lineages, block %in% blocks)
  blk <- levels(droplevels(d$block))
  n_ds <- sapply(blk, function(b) n_distinct(d$dataset[d$block == b]))
  single <- length(lineages) == 1
  panels <- lapply(seq_along(blk), function(i)
    plot_block(filter(d, block == blk[i]), blk[i], show_y = i == 1, single = single, show_sig = show_sig))
  n_cols <- sum(n_ds) * length(lineages)
  ## absolute sizes of the dot area so that column and row spacing is the same
  ## in every layout:
  ## ggplot leaves 0.6 units of margin on each side, so a panel with k columns
  ## takes (k - 1 + 1.2) units; 1 unit = SPACING inches
  w_blk <- if (single) (n_ds + 0.2) * SPACING
           else n_ds * 2.2 * SPACING + (n_ds - 1) * 0.02
  w <- sum(w_blk) + 2.7
  cap <- make_caption(d, w, show_sig)
  h <- PANEL_H + 0.9 + (if (single) 1.5 else 0.9) + 0.15 * (lengths(regmatches(cap, gregexpr("\n", cap))) + 1)
  p <- wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect", widths = unit(w_blk, "in"), heights = unit(PANEL_H, "in")) +
    plot_annotation(title = title, caption = cap,
                    theme = theme(plot.title = element_text(size = 12),
                                  plot.caption = element_text(size = 7, colour = "grey35", hjust = 0, lineheight = 1.1))) &
    theme(legend.justification = "top", legend.box.just = "left",
          plot.margin = margin(4, 2, 4, 2))                       # smaller gap between blocks
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  cat(sprintf("%-22s columnas = %2d  ancho = %.1f in\n", name, n_cols, w))
}

both <- c("Colitis", "Cancer")

## ---- CD4 and CD8 side by side in a single figure ----
## [ CD4: Colitis | Cancer ]   [ CD8: Colitis | Cancer ], same column spacing
make_side_by_side <- function(name, title, show_sig = TRUE) {
  d <- d_all
  blk <- levels(droplevels(d$block))
  n_ds <- sapply(blk, function(b) n_distinct(d$dataset[d$block == b]))
  w_blk <- (n_ds + 0.2) * SPACING
  panels <- list()
  for (lin in c("CD4", "CD8")) for (i in seq_along(blk))
    panels[[length(panels) + 1]] <- plot_block(filter(d, lineage == lin, block == blk[i]), blk[i],
                                               show_y = lin == "CD4" && i == 1, single = TRUE, show_sig = show_sig)
  spanner <- function(lab) ggplot() + annotate("text", x = 0.5, y = 0.3, label = lab, fontface = "bold", size = 4.6) +
    annotate("segment", x = 0.02, xend = 0.98, y = 0, yend = 0, linewidth = 0.5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.05, 0.8), expand = c(0, 0)) +
    theme_void()
  gap <- 0.35                                                    # gap between the CD4 and CD8 blocks
  design <- c(area(1, 1, 1, 2), area(1, 4, 1, 5),
              area(2, 1), area(2, 2), area(2, 4), area(2, 5))
  p <- wrap_plots(c(list(spanner("CD4"), spanner("CD8")), panels), design = design) +
    plot_layout(guides = "collect", widths = unit(c(w_blk, gap, w_blk), "in"),
                heights = unit(c(0.35, PANEL_H), "in")) +
    plot_annotation(title = title, caption = make_caption(d, 2 * sum(w_blk) + gap + 2.7, show_sig),
                    theme = theme(plot.title = element_text(size = 12),
                                  plot.caption = element_text(size = 7, colour = "grey35", hjust = 0, lineheight = 1.1))) &
    theme(legend.justification = "top", legend.box.just = "left", plot.margin = margin(4, 2, 4, 2))
  w <- 2 * sum(w_blk) + gap + 2.9
  cap <- make_caption(d, w, show_sig)
  h <- PANEL_H + 0.35 + 0.9 + 1.5 + 0.15 * (lengths(regmatches(cap, gregexpr("\n", cap))) + 1)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  cat(sprintf("%-22s ancho = %.1f in\n", name, w))
}
## without significance marks
make_side_by_side("SupplFigure9_B_crossdataset_dotmatrix", "Disease vs control, per module", show_sig = FALSE)
