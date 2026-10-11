###############################################################
## T/NK/ILC UMAP and composition change in CRC (Healthy vs Cancer)
##
## Figure  : Fig. 8A (UMAP and legend) and Fig. 8B (composition change)
## Input   : CRC_RDS
## Output  : Figure8_AB_CRC_UMAP_composition.{png,pdf} (+ composition tables)
###############################################################

source("R/00_paths.R")

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(cowplot)

###############################################################
## CONFIGURATION
###############################################################

INPUT_RDS   <- CRC_RDS
OUTPUT_DIR  <- DIR_F8
FIGURE_NAME <- "Figure8_AB_CRC_UMAP_composition"

FDR_THRESHOLD <- 0.05


###############################################################
## LOAD OBJECT
###############################################################

if (!file.exists(INPUT_RDS)) {
  stop(
    "Input file not found: ", INPUT_RDS,
    "\nPlace the Seurat object in the 'data/' folder or update INPUT_RDS above."
  )
}

Tcells <- readRDS(INPUT_RDS)

required_cols <- c("ClusterFull", "Grupo", "orig.ident")
missing_cols <- setdiff(required_cols, colnames(Tcells@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_cols, collapse = ", "))
}

###############################################################
## MAJOR GROUP ASSIGNMENT
###############################################################

ClusterFull <- Tcells$ClusterFull

Tcells$MajorGroup <- dplyr::case_when(

  ClusterFull %in% c(
    "cTNI01 (CD4+ IL7R+)",
    "cTNI02 (CD4+ IL7R+SELL+)",
    "cTNI03 (CD4+ IL7R+HSP+)",
    "cTNI04 (CD4+ IL7R+CCL5+)"
  ) ~ "CD4 Conventional",

  ClusterFull %in% c("cTNI05 (CD4+ IL17+)") ~ "CD4 Th17",

  ClusterFull %in% c("cTNI06 (CD4+ TFH)") ~ "CD4 TFH",

  ClusterFull %in% c("cTNI07 (CD4+ CXCL13+)") ~ "CD4 CXCL13+",

  ClusterFull %in% c(
    "cTNI08 (CD4+ Treg)",
    "cTNI09 (CD4+ Treg prolif)"
  ) ~ "CD4 Treg",

  ClusterFull %in% c(
    "cTNI10 (CD8+ IL7R+)",
    "cTNI11 (CD8+GZMK+)",
    "cTNI12 (CD8+ IL7R+)"
  ) ~ "CD8 Conventional",

  ClusterFull %in% c("cTNI13 (CD8+ T IL17+)") ~ "CD8 Tc17",

  ClusterFull %in% c(
    "cTNI14 (CD8+ CXCL13+)",
    "cTNI15 (CD8+ CXCL13+ HSP+)",
    "cTNI16 (CD8+ CXCL13+ prolif)"
  ) ~ "CD8 Exhausted",

  ClusterFull %in% c(
    "cTNI17 (gd-like T)",
    "cTNI18 (gd-like T PDCD1+)",
    "cTNI19 (gd-like T prolif)"
  ) ~ "γδ T",

  ClusterFull %in% c(
    "cTNI20 (PLZF+ T)",
    "cTNI21 (PLZF+ T prolif)"
  ) ~ "Innate-like T",

  ClusterFull %in% c(
    "cTNI22 (cTNI22)",
    "cTNI23 (NK CD16A+)",
    "cTNI24 (NK GZMK+)",
    "cTNI25 (NK XCL1+)"
  ) ~ "NK",

  ClusterFull %in% c("cTNI26 (ILC3)") ~ "ILC3"
)

if (any(is.na(Tcells$MajorGroup))) {
  unassigned <- unique(Tcells$ClusterFull[is.na(Tcells$MajorGroup)])
  stop(
    "Unassigned ClusterFull value(s): ",
    paste(unassigned, collapse = ", ")
  )
}

Tcells$MajorGroup <- factor(
  Tcells$MajorGroup,
  levels = c(
    "CD4 Conventional", "CD4 Th17", "CD4 TFH", "CD4 CXCL13+", "CD4 Treg",
    "CD8 Conventional", "CD8 Tc17", "CD8 Exhausted",
    "γδ T", "Innate-like T", "NK", "ILC3"
  )
)

###############################################################
## PAIRED PATIENTS (present in both Normal and Tumor)
###############################################################

tab <- table(Tcells$orig.ident, Tcells$Grupo)
paired_ids <- rownames(tab)[apply(tab > 0, 1, all)]

cat("\nPaired patients:", length(paired_ids), "\n")

###############################################################
## MAJOR GROUP COLOR PALETTE
##
## 12 vivid, distinct hues (tab20-family categorical palette, standard in
## genomics/UMAP visualization) instead of within-lineage color shades -
## the previous CD4-reds / CD8-blues scheme put same-hue subsets too close
## together to tell apart at a glance. Note: with 12 simultaneous
## categories, no palette can guarantee full colorblind-safe separation
## for every pair (a hard limit of color space, not of this choice) -
## every panel using this palette also carries a text legend/axis label,
## so identity never depends on color alone.
###############################################################

major_cols <- c(
  "CD4 Conventional" = "#D62728",
  "CD4 Th17"         = "#FF7F0E",
  "CD4 TFH"          = "#B8860B",
  "CD4 CXCL13+"      = "#B0389C",
  "CD4 Treg"         = "#1FB895",  # jade/teal - NOT purple/blue (that means CD8 elsewhere in this paper); kept well apart from gamma-delta T's green below

  "CD8 Conventional" = "#1F77B4",
  "CD8 Tc17"         = "#00B4FB",
  "CD8 Exhausted"    = "#683B92",

  "γδ T"          = "#095F00",  # dark grass green - distinct from CD4 Treg's jade above

  "Innate-like T"    = "#BCBD22",

  "NK"               = "#8C564B",

  "ILC3"             = "#7F7F7F"
)

###############################################################
## PANEL A: Normal UMAP
###############################################################

normal_obj <- subset(Tcells, subset = Grupo == "Normal")

pA <- DimPlot(
  normal_obj,
  reduction = "umap",
  group.by = "MajorGroup",
  cols = major_cols,
  pt.size = 0.15
) +
  ggtitle(paste0("Healthy (n=", ncol(normal_obj), " cells)")) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

###############################################################
## PANEL B: Tumor UMAP
###############################################################

tumor_obj <- subset(Tcells, subset = Grupo == "Tumor")

pB <- DimPlot(
  tumor_obj,
  reduction = "umap",
  group.by = "MajorGroup",
  cols = major_cols,
  pt.size = 0.15
) +
  ggtitle(paste0("Cancer (n=", ncol(tumor_obj), " cells)")) +
  coord_fixed() +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

###############################################################
## SHARED LEGEND
###############################################################

legend <- cowplot::get_legend(
  DimPlot(
    Tcells,
    reduction = "umap",
    group.by = "MajorGroup",
    cols = major_cols,
    pt.size = 0.1
  ) +
    guides(color = guide_legend(override.aes = list(size = 5.5))) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 20)
    )
)

###############################################################
## PER-PATIENT COMPOSITION (% of cells per MajorGroup)
###############################################################

comp <- Tcells@meta.data %>%
  filter(orig.ident %in% paired_ids) %>%
  group_by(orig.ident, Grupo, MajorGroup) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(orig.ident, Grupo) %>%
  mutate(Percent = 100 * N / sum(N)) %>%
  ungroup()

###############################################################
## PAIRED TUMOR vs NORMAL COMPARISON PER SUBSET
##
## paired_delta() pivots one MajorGroup's per-patient percentages
## to wide (Normal/Tumor columns) and returns both the paired
## Wilcoxon test result and the per-patient delta - used once per
## subset instead of the two separate loops the original script ran.
###############################################################

paired_delta <- function(group_name, comp_df) {

  wide <- comp_df %>%
    filter(MajorGroup == group_name) %>%
    select(orig.ident, Grupo, Percent) %>%
    pivot_wider(names_from = Grupo, values_from = Percent, values_fill = 0)

  if (!all(c("Normal", "Tumor") %in% colnames(wide))) {
    return(NULL)
  }

  wide$Delta <- wide$Tumor - wide$Normal
  wide$MajorGroup <- group_name

  pval <- suppressWarnings(
    wilcox.test(wide$Tumor, wide$Normal, paired = TRUE)$p.value
  )

  list(
    stats = data.frame(
      MajorGroup = group_name,
      Delta = mean(wide$Delta),
      Pvalue = pval
    ),
    per_patient = wide
  )
}

paired_results <- lapply(unique(comp$MajorGroup), paired_delta, comp_df = comp)
paired_results <- Filter(Negate(is.null), paired_results)

res      <- bind_rows(lapply(paired_results, `[[`, "stats"))
delta_df <- bind_rows(lapply(paired_results, `[[`, "per_patient"))

###############################################################
## MULTIPLE-TESTING CORRECTION (FDR) AND SIGNIFICANT SUBSETS
###############################################################

res$FDR <- p.adjust(res$Pvalue, method = "BH")

res_sig <- res %>%
  filter(FDR < FDR_THRESHOLD) %>%
  arrange(Delta)

###############################################################
## SIGNIFICANCE STARS (FDR-based)
###############################################################

fdr_stars <- function(fdr) {
  dplyr::case_when(
    is.na(fdr)   ~ "",
    fdr < 0.001  ~ "***",
    fdr < 0.01   ~ "**",
    fdr < 0.05   ~ "*",
    TRUE         ~ "ns"
  )
}

###############################################################
## MEAN +/- SEM PER SUBSET (all subsets shown in Panel C)
###############################################################

summary_df <- delta_df %>%
  group_by(MajorGroup) %>%
  summarise(
    Mean = mean(Delta),
    SEM = sd(Delta) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(res %>% select(MajorGroup, FDR), by = "MajorGroup") %>%
  arrange(Mean) %>%
  mutate(
    Stars      = fdr_stars(FDR),
    # error bar drawn only beyond the bar's tip (away from zero), not back across it
    err_min    = ifelse(Mean >= 0, Mean, Mean - SEM),
    err_max    = ifelse(Mean >= 0, Mean + SEM, Mean),
    star_x     = Mean + sign(Mean) * (SEM + 0.6),
    star_hjust = ifelse(Mean >= 0, 0, 1)
  )

summary_df$Label <- factor(summary_df$MajorGroup, levels = summary_df$MajorGroup)

###############################################################
## PANEL C: compositional shift (Tumor - Normal), all subsets
###############################################################

pC <- ggplot(summary_df, aes(x = Mean, y = Label)) +

  geom_col(
    aes(fill = MajorGroup),
    width = 0.75,
    alpha = 0.85,
    colour = "black",
    linewidth = 0.3
  ) +

  geom_errorbar(
    aes(xmin = err_min, xmax = err_max),
    width = 0.20,
    linewidth = 0.4
  ) +

  geom_text(
    aes(x = star_x, label = Stars, hjust = star_hjust),
    size = 6,
    fontface = "bold",
    vjust = 0.5
  ) +

  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7) +

  scale_fill_manual(values = major_cols) +

  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +

  guides(fill = "none") +

  labs(
    title = expression(Delta ~ composition ~ "(%)"),
    subtitle = paste0(
      "Cancer - Healthy per patient (paired, n=", length(paired_ids),
      "); *** FDR<0.001, ** FDR<0.01, * FDR<0.05, ns = not significant"
    ),
    x = expression(Delta ~ composition ~ "(Cancer - Healthy, %)"),
    y = NULL
  ) +

  theme_classic() +

  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12)
  )

###############################################################
## ASSEMBLE FIGURE
###############################################################

top_panels <- pA + pB

# narrow pC horizontally (padded with blank columns either side) so it
# reads closer to square instead of stretching across the full width
pC_row <- plot_grid(NULL, pC, NULL, ncol = 3, rel_widths = c(0.22, 0.56, 0.22))

final_plot <- plot_grid(
  top_panels,
  legend,
  pC_row,
  ncol = 1,
  rel_heights = c(1.0, 0.35, 1.20)
)

print(final_plot)

###############################################################
## SAVE FIGURE
###############################################################

ggsave(
  filename = file.path(OUTPUT_DIR, paste0(FIGURE_NAME, ".pdf")),
  plot = final_plot,
  width = 14,
  height = 12
)

ggsave(
  filename = file.path(OUTPUT_DIR, paste0(FIGURE_NAME, ".png")),
  plot = final_plot,
  width = 14,
  height = 12,
  dpi = 300
)

###############################################################
## SAVE SUMMARY TABLES
###############################################################

write.csv(
  res,
  file.path(DIR_TAB, "CRC_composition_all_populations.csv"),
  row.names = FALSE
)

write.csv(
  res_sig,
  file.path(DIR_TAB, "CRC_composition_significant_populations.csv"),
  row.names = FALSE
)

###############################################################
## SESSION INFO (for reproducibility / methods reporting)
###############################################################


###############################################################
## SUMMARY
###############################################################

cat("\nPaired patients:", length(paired_ids), "\n")
cat("\nSignificant subsets (FDR <", FDR_THRESHOLD, "):", nrow(res_sig), "\n")
print(res_sig[order(res_sig$Delta, decreasing = TRUE), ])

cat("\nTotal subsets analyzed:", nrow(res), "\n")
cat(
  "\nFigure saved to:\n",
  file.path(OUTPUT_DIR, paste0(FIGURE_NAME, ".pdf")), "\n",
  file.path(OUTPUT_DIR, paste0(FIGURE_NAME, ".png")), "\n"
)
