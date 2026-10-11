###############################################################
## All-gene heatmaps of the 9 processes, CRC CD8 (4-column layout)
##
## Figure  : Suppl. Fig. 8C
## Input   : CRC_RDS, SIG_HUMAN
## Output  : SupplFigure8_C_CRC_CD8_allgenes_4columns.{png,pdf}
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(Seurat); library(Matrix); library(dplyr)
                   library(ComplexHeatmap); library(circlize); library(grid) })
OUT_DIR <- DIR_S8
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sig_map <- function(csv, genes_in_obj) {
  gl <- read.csv(csv, stringsAsFactors = FALSE)
  s2g <- setNames(vector("list", length(MODULES)), MODULES)
  for (i in seq_len(nrow(gl))) for (sg in trimws(strsplit(gl$signatures[i], ";")[[1]]))
    if (sg %in% MODULES) s2g[[sg]] <- c(s2g[[sg]], gl$gene[i])
  lapply(s2g, function(g) intersect(unique(g), genes_in_obj))
}
FONT <- "Arial"
## spacing (mm) between the heatmap and the labels
GROUP_GAP_MM  <- 6     # heatmap -> Healthy/Cancer rule
GROUP_TXT_MM  <- 3.5   # rule -> Healthy/Cancer text
TITLE_GAP_MM  <- 17    # group text -> process title
STRIP_GAP_MM  <- 5     # H/C strip -> heatmap
GROUPS_GAP_MM <- 7     # Healthy -> Cancer
HC_W_MM <- 7
col_fun <- colorRamp2(seq(-2, 2, length.out = 100), colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100))
## keep the Z-score legend away from the gene names (default is 2 mm)
ht_opt$HEATMAP_LEGEND_PADDING <- unit(12, "mm")

###############################################################
## matrices
###############################################################
all_matrix <- function(obj, genes, donor, case) {
  obj$gd <- paste(obj$Grupo, obj[[donor, drop = TRUE]], sep = "-")
  avg <- as.matrix(AverageExpression(obj, group.by = "gd", features = genes, assays = "RNA", layer = "data")$RNA)
  colnames(avg) <- gsub("_", "-", colnames(avg))
  avg <- avg[intersect(genes, rownames(avg)), , drop = FALSE]
  cols <- colnames(avg)
  avg <- avg[, c(cols[startsWith(cols, "Normal-")], cols[startsWith(cols, paste0(case, "-"))]), drop = FALSE]
  is_ctrl <- startsWith(colnames(avg), "Normal-")
  avg <- avg[rowSums(avg) > 0 & apply(avg, 1, sd) > 0, , drop = FALSE]   # fuera genes no detectados
  mu <- rowMeans(avg); sdv <- apply(avg, 1, sd)
  clip <- function(m) { m[is.na(m)] <- 0; pmin(pmax(m, -2), 2) }
  summ <- clip(cbind(H = (rowMeans(avg[, is_ctrl, drop = FALSE]) - mu) / sdv,
                     C = (rowMeans(avg[, !is_ctrl, drop = FALSE]) - mu) / sdv))
  trend <- summ[, "C"] - summ[, "H"]; ord <- order(trend, decreasing = TRUE)
  list(patient = clip((avg - mu) / sdv)[ord, , drop = FALSE], summary = summ[ord, , drop = FALSE], is_ctrl = is_ctrl,
       genes = data.frame(rank = seq_along(ord), gene = rownames(avg)[ord], trend_C_minus_H = trend[ord]))
}

###############################################################
## drawing of one panel (used both alone and inside the grid)
###############################################################
make_hts <- function(res, idx, body_h, cw, gene_pt, hc_pt, show_leg) {
  sm <- res$summary; colnames(sm) <- c("H", "C"); pm <- res$patient
  grupos <- factor(ifelse(res$is_ctrl, "Healthy", "Cancer"), levels = c("Healthy", "Cancer"))
  ht_s <- Heatmap(sm, name = paste0("HC", idx), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
                  show_heatmap_legend = FALSE, rect_gp = gpar(col = NA), border = "black", border_gp = gpar(lwd = 2),
                  show_row_names = FALSE, show_column_names = FALSE,   # "HC" is drawn separately, at the same size as Healthy/Cancer
                  width = unit(HC_W_MM * 2, "mm"), height = unit(body_h, "mm"))
  ht_m <- Heatmap(pm, name = paste0("Z", idx), col = col_fun, column_split = grupos,
                  column_gap = unit(GROUPS_GAP_MM, "mm"), cluster_rows = FALSE, cluster_columns = FALSE,
                  rect_gp = gpar(col = NA), border = "black", border_gp = gpar(lwd = 2),
                  show_column_names = FALSE, show_heatmap_legend = show_leg,
                  row_names_gp = gpar(fontsize = gene_pt, fontface = "italic", fontfamily = FONT),
                  column_title = rep("", nlevels(grupos)),
                  width = unit(cw * ncol(pm), "mm"), height = unit(body_h, "mm"),
                  heatmap_legend_param = list(title = "Z-score", at = -2:2,
                                              title_gp = gpar(fontsize = 15, fontface = "bold", fontfamily = FONT),
                                              labels_gp = gpar(fontsize = 14, fontfamily = FONT),
                                              legend_height = unit(35, "mm"), grid_width = unit(6, "mm")))
  list(ht = ht_s + ht_m, grupos = grupos, nm_s = paste0("HC", idx), nm_m = paste0("Z", idx))
}

## mid_mm: if given, the text (HC / Healthy / Cancer) is centred at that height above
## the heatmap instead of sitting on its rule, so that it falls halfway between the
## group rule (below) and the title underline (above)
decorate_panel <- function(x, group_pt, gap = GROUP_GAP_MM, txt = GROUP_TXT_MM, mid_mm = NULL) {
  ty  <- if (is.null(mid_mm)) unit(gap + txt, "mm") else unit(mid_mm, "mm")
  tvj <- if (is.null(mid_mm)) 0 else 0.5
  ## summary strip: white divider and the "HC" label above it, at the same size
  ## and height as Healthy/Cancer (no underline)
  decorate_heatmap_body(x$nm_s, {
    grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 0.8))
    grid.text("HC", y = unit(1, "npc") + ty, vjust = tvj,
              gp = gpar(fontsize = group_pt, fontfamily = FONT))
  })
  for (k in 1:2) decorate_heatmap_body(x$nm_m, {
    grid.lines(x = unit(c(0, 1), "npc"), y = unit(1, "npc") + unit(gap, "mm"), gp = gpar(lwd = 1.2))
    grid.text(levels(x$grupos)[k], y = unit(1, "npc") + ty, vjust = tvj,
              gp = gpar(fontsize = group_pt, fontfamily = FONT))
  }, column_slice = k)
  xl <- xr <- yt <- NULL
  decorate_heatmap_body(x$nm_s, { l <- deviceLoc(unit(0, "npc"), unit(1, "npc")); xl <<- l$x; yt <<- l$y })
  decorate_heatmap_body(x$nm_m, { xr <<- deviceLoc(unit(1, "npc"), unit(0, "npc"))$x }, column_slice = 2)
  list(xl = xl, xr = xr, yt = yt)
}

draw_title <- function(loc, title, title_pt, wrap_at = 0.0073, above_mm = GROUP_GAP_MM + GROUP_TXT_MM + TITLE_GAP_MM, under_mm = 4) {
  xlc <- convertX(loc$xl, "in", valueOnly = TRUE); xrc <- convertX(loc$xr, "in", valueOnly = TRUE)
  xc <- xlc + (xrc - xlc) / 2
  if (nchar(title) * title_pt * wrap_at > (xrc - xlc)) {
    w <- strsplit(title, " ")[[1]]
    cut <- which.min(sapply(seq_len(length(w) - 1), function(k)
      abs(nchar(paste(w[1:k], collapse = " ")) - nchar(paste(w[-(1:k)], collapse = " ")))))
    title <- paste(paste(w[1:cut], collapse = " "), paste(w[-(1:cut)], collapse = " "), sep = "\n")
  }
  g <- textGrob(title, gp = gpar(fontsize = title_pt, fontface = "bold", lineheight = 0.95, fontfamily = FONT))
  tw <- convertWidth(grobWidth(g), "in", valueOnly = TRUE)
  th <- convertHeight(grobHeight(g), "in", valueOnly = TRUE)
  y <- unit(convertY(loc$yt, "in", valueOnly = TRUE), "in") + unit(above_mm, "mm") + unit(th / 2, "in")
  grid.draw(editGrob(g, x = unit(xc, "in"), y = y))
  grid.lines(x = unit(xc + c(-tw / 2 - 0.12, tw / 2 + 0.12), "in"),
             y = y - unit(th / 2, "in") - unit(under_mm, "mm"), gp = gpar(lwd = 2))
  if (grepl("\n", title)) 11 else 0
}

open_dev <- function(f, w, h, dev) {
  if (dev == "pdf") cairo_pdf(f, width = w, height = h, family = FONT)
  else png(f, width = w, height = h, units = "in", res = 300, family = FONT)
}


###############################################################
## combinada en 4 columnas
###############################################################
COLS4 <- list(c("Glycolysis", "TCA", "Transporters"),
              c("OXPHOS", "Trafficking"),
              c("Methyltransferase", "Tcell_activation"),
              c("Cytokine", "Translation"))
## panels are given a similar height regardless of how many genes they have, so
## the cell height and the gene font adapt: height = n x CELL_H4_MM, capped
## (the largest processes are compressed)
CELL_H4_MM <- 8         # cell height in the small processes
MAXH4_MM <- 220         # maximum body height of a panel
BODY_W4_MM <- 125       # body width (AOM, 5+5 mice)
BODY_W4_CRC_MM <- 160   # CRC, 36 patients per group, wider
PANEL_W4_IN <- 8.0      # column width (leaves room for the gene names)
PANEL_W4_CRC_IN <- 9.2
OVERHEAD4_IN <- 4.2     # header of each panel (title + group labels)
GAP4_IN <- 0.5          # vertical gap between panels of the same column
G_GAP <- 10; G_TXT <- 15; T_GAP <- 24; UND <- 11
GROUP_PT4 <- 40; TITLE_PT4 <- 46; GENE_PT4 <- 7

combined4 <- function(res, stub, body_w = BODY_W4_MM, panel_w = PANEL_W4_IN) {
  ng <- vapply(res, function(x) nrow(x$patient), 1L)
  bh <- lapply(COLS4, function(cc) pmin(MAXH4_MM, ng[cc] * CELL_H4_MM))          # body height, mm
  ph <- lapply(bh, function(v) v / 25.4 + OVERHEAD4_IN)                          # panel height, in
  col_h <- sapply(ph, function(v) sum(v) + GAP4_IN * (length(v) - 1))
  W <- 4 * panel_w + 2.2; H <- max(col_h) + 0.4
  for (dev in c("pdf", "png")) {
    open_dev(paste0(stub, ".", dev), W, H, dev); grid.newpage()
    idx <- 0
    for (j in seq_along(COLS4)) {
      ytop <- 0.2
      for (i in seq_along(COLS4[[j]])) {
        m <- COLS4[[j]][i]; idx <- idx + 1
        bhm <- bh[[j]][i]
        x <- make_hts(res[[m]], idx, bhm, body_w / ncol(res[[m]]$patient),
                      max(4, min(13, (bhm / ng[m]) * 2.4)), GROUP_PT4, show_leg = FALSE)
        upViewport(0)
        pushViewport(viewport(x = unit((j - 1) * panel_w, "in"), y = unit(1, "npc") - unit(ytop, "in"),
                              width = unit(panel_w, "in"), height = unit(ph[[j]][i], "in"),
                              just = c("left", "top")))
        draw(x$ht, newpage = FALSE,
             padding = unit(c(2, 2, 10 + G_GAP + G_TXT + T_GAP, 2), "mm"), ht_gap = unit(STRIP_GAP_MM, "mm"))
        ## text centred between the group rule (G_GAP) and the title underline
        loc <- decorate_panel(x, GROUP_PT4, G_GAP, G_TXT,
                              mid_mm = (2 * G_GAP + G_TXT + T_GAP - UND) / 2)
        upViewport(0)
        draw_title(loc, MOD_LABELS[[m]], TITLE_PT4, above_mm = G_GAP + G_TXT + T_GAP, under_mm = UND)
        ytop <- ytop + ph[[j]][i] + GAP4_IN
      }
    }
    upViewport(0)
    pushViewport(viewport(x = unit(4 * panel_w + 0.3, "in"), y = unit(0.8, "npc"),
                          width = unit(1.8, "in"), height = unit(0.3, "npc"), just = c("left", "centre")))
    draw(Legend(col_fun = col_fun, title = "Z-score", at = -2:2,
                title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = FONT),
                labels_gp = gpar(fontsize = 22, fontfamily = FONT),
                legend_height = unit(70, "mm"), grid_width = unit(11, "mm")))
    dev.off()
  }
  cat(basename(stub), "OK\n")
}

###############################################################
## CRC Pelka
###############################################################
o <- readRDS(CRC_RDS)
tb <- table(o$orig.ident, o$Grupo)
o <- subset(o, subset = orig.ident %in% rownames(tb)[tb[, "Normal"] > 0 & tb[, "Tumor"] > 0])
o <- SetAssayData(o, assay = "RNA", layer = "data",
                  new.data = log1p(GetAssayData(o, assay = "RNA", layer = "counts")))
s2g <- sig_map(SIG_HUMAN, rownames(o))
for (lin in "CD8") {
  obj <- subset(o, subset = ClusterMidway == paste0("T", lin))
  res <- lapply(MODULES, function(m) all_matrix(obj, s2g[[m]], "orig.ident", "Tumor")); names(res) <- MODULES
  combined4(res, file.path(OUT_DIR, sprintf("SupplFigure8_C_CRC_%s_allgenes_4columns", lin)),
            body_w = BODY_W4_CRC_MM, panel_w = PANEL_W4_CRC_IN)
}
