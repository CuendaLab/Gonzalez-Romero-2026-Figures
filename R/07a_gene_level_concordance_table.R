###############################################################
## Direction of change of every gene of the 9 processes in the three datasets
##
## Figure  : input for Fig. 8F
## Input   : AOM_RDS, CRC_RDS, PROT_XLSX, SIG_MOUSE
## Output  : tables/Concordance_gene_level_table.csv (+ summaries)
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(Seurat); library(Matrix); library(dplyr); library(tidyr); library(readxl) })

OUT_DIR <- DIR_TAB

## ---- mouse signature with its human ortholog and modules (Ribosome excluded) ----
sig <- read.csv(SIG_MOUSE,
                stringsAsFactors = FALSE) |>
  mutate(modules = sapply(strsplit(signatures, ";"), function(x) paste(intersect(trimws(x), MODULES), collapse = ";"))) |>
  filter(modules != "") |>
  select(mouse = gene, human = human_ortholog, modules)

## ---- utilidades ----
group_means <- function(M, key) {               # M: genes x cells (sparse); group means
  f <- factor(key); ind <- sparse.model.matrix(~ 0 + f); colnames(ind) <- levels(f)
  as.matrix(M %*% ind) / rep(as.numeric(table(f)), each = nrow(M))
}

## ---------------- AOM-DSS ----------------
obj <- readRDS(AOM_RDS)
obj <- subset(obj, subset = ClusterMidway %in% c("TCD4", "TCD8") & Grupo %in% c("Normal", "CAC"))
key <- paste(obj$sample, sub("^T", "", obj$ClusterMidway), obj$Grupo, sep = "___")
keep <- names(which(table(key) >= 20))
cells <- key %in% keep
gm <- group_means(GetAssayData(obj, assay = "RNA", layer = "data")[, cells], key[cells])
rm(obj); invisible(gc())
info <- do.call(rbind, strsplit(colnames(gm), "___")); colnames(info) <- c("sample", "lineage", "group")
aom <- bind_rows(lapply(c("CD4", "CD8"), function(lin) {
  a <- gm[, info[, "lineage"] == lin & info[, "group"] == "Normal", drop = FALSE]
  b <- gm[, info[, "lineage"] == lin & info[, "group"] == "CAC", drop = FALSE]
  ok <- rowSums(a) + rowSums(b) > 0
  data.frame(gene = rownames(gm)[ok], lineage = lin,
             delta = rowMeans(b[ok, , drop = FALSE]) - rowMeans(a[ok, , drop = FALSE]),
             p = sapply(which(ok), function(i) suppressWarnings(wilcox.test(b[i, ], a[i, ])$p.value)))
}))

## ---------------- CRC Pelka ----------------
obj <- readRDS(CRC_RDS)
obj <- subset(obj, subset = ClusterMidway %in% c("TCD4", "TCD8"))
key <- paste(obj$orig.ident, sub("^T", "", obj$ClusterMidway), obj$Grupo, sep = "___")
gm <- group_means(GetAssayData(obj, assay = "RNA", layer = "counts"), key)
rm(obj); invisible(gc())
info <- do.call(rbind, strsplit(colnames(gm), "___")); colnames(info) <- c("donor", "lineage", "spec")
crc <- bind_rows(lapply(c("CD4", "CD8"), function(lin) {
  dn <- intersect(info[info[, "lineage"] == lin & info[, "spec"] == "Normal", "donor"],
                  info[info[, "lineage"] == lin & info[, "spec"] == "Tumor", "donor"])
  N <- gm[, paste(dn, lin, "Normal", sep = "___"), drop = FALSE]
  T <- gm[, paste(dn, lin, "Tumor",  sep = "___"), drop = FALSE]
  D <- T - N
  ok <- rowSums(N) + rowSums(T) > 0
  data.frame(gene = rownames(gm)[ok], lineage = lin, delta = rowMeans(D[ok, , drop = FALSE]),
             p = apply(D[ok, , drop = FALSE], 1, function(d) suppressWarnings(wilcox.test(d)$p.value)))
}))

## ---------------- mLN proteomics ----------------
raw <- suppressMessages(read_excel(PROT_XLSX,
                                   sheet = "Complete analysis", col_names = FALSE))
hdr <- as.character(unlist(raw[2, ])); cn_cols <- c(15:23, 26:34)
body <- raw[-(1:2), ]
CN <- suppressWarnings(apply(as.matrix(body[, cn_cols]), 2, as.numeric)); colnames(CN) <- hdr[cn_cols]
CN[!is.finite(CN) | CN <= 0] <- NA
pg_genes <- as.character(body[[2]])
lin_s <- sub(".*_(CD[48])_.*", "\\1", colnames(CN)); grp_s <- sub(".*_CD[48]_(\\w+)_WT$", "\\1", colnames(CN))
ok_q <- Reduce(`&`, lapply(c("CD4", "CD8"), function(l) Reduce(`&`, lapply(c("CONTROL", "CANCER"), function(g)
  rowSums(!is.na(CN[, lin_s == l & grp_s == g, drop = FALSE])) >= 2))))
prot <- bind_rows(lapply(c("CD4", "CD8"), function(lin) {
  A <- CN[, lin_s == lin & grp_s == "CONTROL"]; B <- CN[, lin_s == lin & grp_s == "CANCER"]
  ix <- which(ok_q & !is.na(pg_genes))
  data.frame(pg = pg_genes[ix], lineage = lin,
             delta = log2(rowMeans(B[ix, ], na.rm = TRUE) / rowMeans(A[ix, ], na.rm = TRUE)),
             p = sapply(ix, function(i) tryCatch(t.test(log10(B[i, ]), log10(A[i, ]))$p.value, error = function(e) NA_real_)))
})) |>
  mutate(gene = strsplit(pg, ";")) |> unnest(gene) |>
  group_by(lineage, gene) |> slice(1) |> ungroup()   # one symbol -> one group of proteins

## ---------------- concordancia ----------------
concord <- function(x, y) {        # x, y: Δ emparejados
  ux <- x > 0; uy <- y > 0
  p1 <- mean(ux); p2 <- mean(uy)
  tab <- table(factor(ux, c(TRUE, FALSE)), factor(uy, c(TRUE, FALSE)))
  data.frame(n = length(x), pct_same_direction = 100 * mean(ux == uy),
             pct_expected_by_chance = 100 * (p1 * p2 + (1 - p1) * (1 - p2)),
             pct_up_1 = 100 * p1, pct_up_2 = 100 * p2,
             odds_ratio = unname(fisher.test(tab)$estimate), fisher_p = fisher.test(tab)$p.value,
             spearman_rho = cor(x, y, method = "spearman"))
}

## columns: delta_1/p_1 = first platform of the pair, delta_2/p_2 = second
pairs <- list()
## AOM scRNA vs proteomics (mouse, same gene symbol)
pairs$aom_prot <- inner_join(aom |> rename(delta_1 = delta, p_1 = p),
                             prot |> select(gene, lineage, delta_2 = delta, p_2 = p),
                             by = c("gene", "lineage")) |>
  inner_join(sig |> select(gene = mouse, human, modules), by = "gene") |>
  mutate(comparison = "AOM-DSS scRNA (mouse) vs mLN proteomics (mouse)")
## mouse -> human orthologs of the signature (only genes of the 9 processes)
orth <- sig |> select(mouse, human, modules)
## proteomics vs CRC (mouse -> human)
pairs$prot_crc <- prot |> select(gene, lineage, delta_1 = delta, p_1 = p) |>
  inner_join(orth, by = c("gene" = "mouse")) |>
  inner_join(crc |> rename(delta_2 = delta, p_2 = p), by = c("human" = "gene", "lineage")) |>
  mutate(comparison = "mLN proteomics (mouse) vs CRC Pelka scRNA (human)")
## AOM scRNA vs CRC scRNA (mouse -> human)
pairs$aom_crc <- aom |> rename(delta_1 = delta, p_1 = p) |>
  inner_join(orth, by = c("gene" = "mouse")) |>
  inner_join(crc |> rename(delta_2 = delta, p_2 = p), by = c("human" = "gene", "lineage")) |>
  mutate(comparison = "AOM-DSS scRNA (mouse) vs CRC Pelka scRNA (human)")

genes_tab <- bind_rows(pairs) |> select(comparison, lineage, gene, human, modules, delta_1, p_1, delta_2, p_2)
write.csv(genes_tab, file.path(OUT_DIR, "Concordance_gene_level_table.csv"), row.names = FALSE)

summ <- function(d, universe) d |> group_by(comparison, lineage) |>
  group_modify(~ concord(.x$delta_1, .x$delta_2)) |> ungroup() |> mutate(universe = universe, .before = 1)
res <- bind_rows(
  summ(genes_tab, "9-process genes"),
  summ(genes_tab |> filter(p_1 < 0.05, p_2 < 0.05), "9-process genes, p<0.05 in both"))
write.csv(res, file.path(OUT_DIR, "Concordance_summary.csv"), row.names = FALSE)

by_mod <- genes_tab |>
  mutate(module = strsplit(modules, ";")) |> unnest(module) |>
  group_by(comparison, lineage, module) |>
  summarise(n = n(), pct_same_direction = 100 * mean((delta_1 > 0) == (delta_2 > 0)),
            pct_up_1 = 100 * mean(delta_1 > 0), pct_up_2 = 100 * mean(delta_2 > 0),
            spearman_rho = cor(delta_1, delta_2, method = "spearman"), .groups = "drop") |>
  mutate(module = factor(MOD_LABELS[module], levels = MOD_LABELS)) |> arrange(comparison, lineage, module)
write.csv(by_mod, file.path(OUT_DIR, "Concordance_by_module.csv"), row.names = FALSE)

options(width = 200)
print(as.data.frame(res), digits = 3)
cat("\n"); print(as.data.frame(by_mod), digits = 3)
