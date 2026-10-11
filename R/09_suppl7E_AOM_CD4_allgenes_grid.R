###############################################################
## All-gene heatmaps of the 9 processes, AOM/DSS CD4 (3x3 grid)
##
## Figure  : Suppl. Fig. 7E
## Input   : AOM_RDS, SIG_MOUSE
## Output  : SupplFigure7_E_AOMDSS_CD4_allgenes_grid.{png,pdf}
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(Seurat); library(Matrix); library(dplyr)
                   library(ComplexHeatmap); library(circlize); library(grid) })
OUT_DIR <- DIR_S7

CELL_W_MAX_MM <- 12       # with 5+5 mice, cap the cell width
CELL_H_MM <- 6            # target cell height
BODY_MAX_MM <- c(170, 290, 235)   # maximum body height per row (each row has its own)
## panel order, by rows
SQUARE_ROWS <- 1                  # rows drawn with square cells (height = width)
ROWS <- list(c("Glycolysis", "TCA", "Transporters"),
             c("OXPHOS", "Methyltransferase", "Cytokine"),
             c("Trafficking", "Tcell_activation", "Translation"))
PANEL_W_IN <- 8.6; OVERHEAD_IN <- 3.2   # titulo + Healthy/Cancer + margenes
GROUP_PT <- 34; TITLE_PT <- 60; HC_PT <- 22
GROUP_GAP_MM <- 2; STRIP_GAP_MM <- 4; GROUPS_GAP_MM <- 7
HC_W_MM <- 10
BODY_W_MM <- 240          # body width, the same in all 9 panels
col_fun <- colorRamp2(seq(-2, 2, length.out = 100), colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100))

o <- readRDS(AOM_RDS)
o <- subset(o, subset = Grupo %in% c("Normal", "CAC"))          # colitis mice are not used here

gl <- read.csv(SIG_MOUSE, stringsAsFactors = FALSE)
s2g <- setNames(vector("list", length(MODULES)), MODULES)
for (i in seq_len(nrow(gl))) for (sg in trimws(strsplit(gl$signatures[i], ";")[[1]]))
  if (sg %in% MODULES) s2g[[sg]] <- c(s2g[[sg]], gl$gene[i])
s2g <- lapply(s2g, function(g) intersect(unique(g), rownames(o)))

all_matrix <- function(obj, genes) {
  obj$gd <- paste(obj$Grupo, obj$sample, sep = "-")
  avg <- as.matrix(AverageExpression(obj, group.by = "gd", features = genes, assays = "RNA", layer = "data")$RNA)
  colnames(avg) <- gsub("_", "-", colnames(avg))
  avg <- avg[intersect(genes, rownames(avg)), , drop = FALSE]
  cols <- colnames(avg)
  avg <- avg[, c(cols[startsWith(cols, "Normal-")], cols[startsWith(cols, "CAC-")]), drop = FALSE]
  is_ctrl <- startsWith(colnames(avg), "Normal-")
  ## drop genes not detected in this lineage (all zeros -> undefined Z, blank row)
  avg <- avg[rowSums(avg) > 0 & apply(avg, 1, sd) > 0, , drop = FALSE]
  mu <- rowMeans(avg); sdv <- apply(avg, 1, sd)
  clip <- function(m) { m[is.na(m)] <- 0; pmin(pmax(m, -2), 2) }
  summ <- clip(cbind(H = (rowMeans(avg[, is_ctrl, drop = FALSE]) - mu) / sdv,
                     C = (rowMeans(avg[, !is_ctrl, drop = FALSE]) - mu) / sdv))
  trend <- summ[, "C"] - summ[, "H"]; ord <- order(trend, decreasing = TRUE)
  list(patient = clip((avg - mu) / sdv)[ord, , drop = FALSE], summary = summ[ord, , drop = FALSE], is_ctrl = is_ctrl)
}

## ---- one panel inside the current viewport ------------------------------
## lab_case: label of the disease group
panel <- function(res, title, idx, body_h, lab_case = "Cancer") {
  grupos <- factor(ifelse(res$is_ctrl, "Healthy", lab_case), levels = c("Healthy", lab_case))
  sm <- res$summary; colnames(sm) <- c("H", substr(lab_case, 1, 1)); pm <- res$patient
  n <- nrow(pm)
  ch <- body_h / n                          # same height within the row
  gene_pt <- max(3, min(13, ch * 2.4))
  cw <- min(BODY_W_MM / ncol(pm), CELL_W_MAX_MM)
  nm <- paste0("Z", idx)
  ht_s <- Heatmap(sm, name = paste0("HC", idx), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                  show_heatmap_legend = FALSE, rect_gp = gpar(col = NA), border = "black", border_gp = gpar(lwd = 2),
                  show_row_names = FALSE, column_names_side = "top", column_names_rot = 0,
                  column_names_centered = TRUE, column_names_gp = gpar(fontsize = HC_PT),
                  width = unit(HC_W_MM * 2, "mm"), height = unit(body_h, "mm"))
  ht_m <- Heatmap(pm, name = nm, col = col_fun, column_split = grupos, column_gap = unit(GROUPS_GAP_MM, "mm"),
                  cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = NA),
                  border = "black", border_gp = gpar(lwd = 2), show_column_names = FALSE,
                  show_heatmap_legend = FALSE, row_names_gp = gpar(fontsize = gene_pt, fontface = "italic"),
                  column_title = rep("", nlevels(grupos)),
                  width = unit(cw * ncol(pm), "mm"), height = unit(body_h, "mm"))
  draw(ht_s + ht_m, newpage = FALSE, padding = unit(c(3, 2, 34, 2), "mm"), ht_gap = unit(STRIP_GAP_MM, "mm"))
  xl <- xr <- yt <- NULL
  decorate_heatmap_body(paste0("HC", idx), { l <- deviceLoc(unit(0, "npc"), unit(1, "npc")); xl <<- l$x; yt <<- l$y })
  decorate_heatmap_body(nm, { xr <<- deviceLoc(unit(1, "npc"), unit(0, "npc"))$x }, column_slice = 2)
  decorate_heatmap_body(paste0("HC", idx),
    grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 0.8)))
  for (k in 1:2) decorate_heatmap_body(nm, {
    grid.lines(x = unit(c(0, 1), "npc"), y = unit(1, "npc") + unit(GROUP_GAP_MM, "mm"), gp = gpar(lwd = 1.5))
    grid.text(levels(grupos)[k], y = unit(1, "npc") + unit(GROUP_GAP_MM + 2.6, "mm"), vjust = 0,
              gp = gpar(fontsize = GROUP_PT, fontface = "bold"))
  }, column_slice = k)
  list(xl = xl, xr = xr, yt = yt)
}

title_at <- function(loc, title, w_in) {
  xlc <- convertX(loc$xl, "in", valueOnly = TRUE); xrc <- convertX(loc$xr, "in", valueOnly = TRUE)
  xc <- xlc + (xrc - xlc) / 2
  if (nchar(title) * TITLE_PT * 0.0073 > (xrc - xlc)) {
    wds <- strsplit(title, " ")[[1]]
    cut <- which.min(sapply(seq_len(length(wds) - 1), function(k)
      abs(nchar(paste(wds[1:k], collapse = " ")) - nchar(paste(wds[-(1:k)], collapse = " ")))))
    title <- paste(paste(wds[1:cut], collapse = " "), paste(wds[-(1:cut)], collapse = " "), sep = "\n")
  }
  g <- textGrob(title, gp = gpar(fontsize = TITLE_PT, fontface = "bold", lineheight = 0.95))
  tw <- convertWidth(grobWidth(g), "in", valueOnly = TRUE)
  th <- convertHeight(grobHeight(g), "in", valueOnly = TRUE)   # 1 o 2 lineas
  y <- unit(convertY(loc$yt, "in", valueOnly = TRUE), "in") + unit(26, "mm") + unit(th / 2 - 0.1, "in")
  grid.draw(editGrob(g, x = unit(xc, "in"), y = y))
  grid.lines(x = unit(xc + c(-tw / 2 - 0.1, tw / 2 + 0.1), "in"),
             y = y - unit(th / 2, "in") - unit(2, "mm"), gp = gpar(lwd = 2.5))
}

for (lin in "CD4") {
  obj <- subset(o, subset = ClusterMidway == paste0("T", lin))
  res <- lapply(MODULES, function(m) all_matrix(obj, s2g[[m]])); names(res) <- MODULES

  ng   <- vapply(res, function(x) nrow(x$patient), 1L)
  rows <- ROWS
  ## body height per panel: shared within the row, except in SQUARE_ROWS, where
  ## each panel uses n_genes x cell width so that the cells come out square
  cwm <- min(BODY_W_MM / ncol(res[[1]]$patient), CELL_W_MAX_MM)
  bh  <- lapply(1:3, function(i) if (i %in% SQUARE_ROWS) ng[rows[[i]]] * cwm
                else setNames(rep(min(BODY_MAX_MM[i], CELL_H_MM * max(ng[rows[[i]]])), 3), rows[[i]]))
  row_h <- sapply(bh, max) / 25.4 + OVERHEAD_IN

  lay <- grid.layout(3, 4, widths = unit(c(rep(PANEL_W_IN, 3), 2.6), "in"),
                     heights = unit(row_h, "in"))
  W <- 3 * PANEL_W_IN + 2.6; H <- sum(row_h)
  for (lab_case in "Cancer") for (dev in c("pdf", "png")) {
    f <- file.path(OUT_DIR, sprintf("SupplFigure7_E_AOMDSS_%s_allgenes_grid.%s", lin, dev))
    if (dev == "pdf") pdf(f, width = W, height = H) else png(f, width = W, height = H, units = "in", res = 300)
    grid.newpage()
    k <- 0
    for (r in 1:3) for (cn in 1:3) {
      k <- k + 1; m <- rows[[r]][cn]
      upViewport(0); pushViewport(viewport(layout = lay))
      pushViewport(viewport(layout.pos.row = r, layout.pos.col = cn))
      loc <- panel(res[[m]], MOD_LABELS[[m]], k, bh[[r]][[m]], lab_case)
      upViewport(0)
      title_at(loc, MOD_LABELS[[m]], PANEL_W_IN)
    }
    upViewport(0); pushViewport(viewport(layout = lay))
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
    draw(Legend(col_fun = col_fun, title = "Z-score", at = -2:2,
                title_gp = gpar(fontsize = 30, fontface = "bold"), labels_gp = gpar(fontsize = 26),
                legend_height = unit(85, "mm"), grid_width = unit(13, "mm")))
    dev.off()
  }
  cat(lin, "OK - filas:", paste(sapply(rows, paste, collapse = "/"), collapse = " | "), "\n")
}
