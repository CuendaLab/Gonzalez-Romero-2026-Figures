###############################################################
## Module scores of the 9 processes in the mLN proteomics (CD4 and CD8)
##
## Figure  : input for Suppl. Fig. 9B
## Input   : PROT_XLSX, SIG_MOUSE
## Output  : tables/Table_proteomics_mLN_module_delta_vs_Basal.csv
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(readxl); library(dplyr) })

SIG_CSV   <- SIG_MOUSE
OUT_DIR   <- DIR_TAB
dir.create(OUT_DIR, showWarnings = FALSE)

MODULES <- c("Glycolysis", "OXPHOS", "Transporters", "Trafficking", "TCA",
             "Tcell_activation", "Cytokine", "Methyltransferase", "Translation")

## ---- read copy numbers (row 2 = header; rows 3+ = proteins) ----
raw <- suppressMessages(read_excel(PROT_XLSX, sheet = "Complete analysis", col_names = FALSE))
hdr <- as.character(unlist(raw[2, ]))
cn_cols <- c(15:23, 26:34)                       # bloque COPY NUMBERS CD4 y CD8
stopifnot(all(grepl("_BR[123]_CD[48]_(CONTROL|COLITIS|CANCER)_WT$", hdr[cn_cols])))

body  <- raw[-(1:2), ]
genes <- as.character(body[[2]])
CN <- suppressWarnings(apply(as.matrix(body[, cn_cols]), 2, as.numeric))
colnames(CN) <- hdr[cn_cols]
CN[!is.finite(CN) | CN <= 0] <- NA

meta <- data.frame(sample = colnames(CN)) |>
  mutate(lineage = sub(".*_(CD[48])_.*", "\\1", sample),
         group   = sub(".*_CD[48]_(\\w+)_WT$", "\\1", sample))

## ---- filter: quantified in >= 2/3 replicates of every group x lineage ----
keep_grp <- with(meta, group %in% c("CONTROL", "CANCER"))
ok <- sapply(split(which(keep_grp), paste(meta$lineage, meta$group)[keep_grp]),
             function(ix) rowSums(!is.na(CN[, ix, drop = FALSE])) >= 2)
keep <- rowSums(ok) == ncol(ok) & !is.na(genes)

## ---- Z-score per protein across the 18 samples ----
L <- log10(CN[keep, ])
Z <- (L - rowMeans(L, na.rm = TRUE)) / (apply(L, 1, sd, na.rm = TRUE) + 1e-9)
gsym <- strsplit(genes[keep], ";")                # protein groups with several gene symbols

## ---- signature -> proteins ----
sig <- read.csv(SIG_CSV, stringsAsFactors = FALSE)
cov <- list()
scores <- sapply(MODULES, function(m) {
  mg  <- sig$gene[sapply(strsplit(sig$signatures, ";"), function(x) m %in% trimws(x))]
  hit <- which(sapply(gsym, function(g) any(g %in% mg)))
  cov[[m]] <<- c(n_sig = length(mg), n_found = length(hit))
  colMeans(Z[hit, , drop = FALSE], na.rm = TRUE)
})
coverage <- do.call(rbind, cov)
cat("Module coverage (proteins found / genes in the signature):\n"); print(coverage)

## Copy numbers per cell: the mean Z of the WHOLE proteome shifts between
## grupos (CD4 sube ~1 SD Basal->Cancer, CD8 baja ~0.7 SD; contenido proteico
## groups. The "relative" version = module score minus that global mean per
## sample (module relative to the rest of the proteome).
background <- colMeans(Z, na.rm = TRUE)
scores_rel <- scores - background

score_long <- rbind(cbind(meta, version = "absolute", background = background, scores),
                    cbind(meta, version = "relative", background = background, scores_rel))
write.csv(score_long, file.path(OUT_DIR, "proteomics_mLN_module_scores_per_sample.csv"), row.names = FALSE)

## ---- Colitis and Cancer vs baseline (day 0), Welch t, BH across the 9 modules ----
res <- do.call(rbind, lapply(c("COLITIS", "CANCER"), function(case) do.call(rbind, lapply(c("absolute", "relative"), function(v)
  do.call(rbind, lapply(c("CD4", "CD8"), function(lin) {
    d <- score_long[score_long$lineage == lin & score_long$version == v, ]
    r <- do.call(rbind, lapply(MODULES, function(m) {
      a <- d[[m]][d$group == "CONTROL"]; b <- d[[m]][d$group == case]
      data.frame(comparison = paste(tools::toTitleCase(tolower(case)), "vs Basal"), version = v, sig = m, comp = lin,
                 n_control = length(a), n_case = length(b), delta_mean = mean(b) - mean(a), p = t.test(b, a)$p.value)
    }))
    r$fdr <- p.adjust(r$p, "BH"); r
  }))))))
res <- cbind(res, coverage[res$sig, ])
write.csv(res, file.path(OUT_DIR, "Table_proteomics_mLN_module_delta_vs_Basal.csv"), row.names = FALSE)
rownames(res) <- NULL
print(res[res$version == "absolute", c("comparison", "comp", "sig", "delta_mean", "p", "fdr")], digits = 3)
