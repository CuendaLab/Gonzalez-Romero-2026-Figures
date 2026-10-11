###############################################################
## UMAP and composition change in mouse colon, AOM/DSS model
##
## Figure  : Suppl. Fig. 7C-D
## Input   : AOM_RDS
## Output  : SupplFigure7_CD_AOMDSS_UMAP_composition.{png,pdf} (+ legend and tables)
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(Seurat); library(Matrix); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

OUT_DIR  <- DIR_S7
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MIN_CELLS_PER_SAMPLE <- 50     # mice with fewer cells are left out of the composition
MIN_CELLS_PSEUDOBULK <- 20     # minimum cells per mouse and lineage for a module score
GROUP_LAB <- c(Normal = "Healthy", CAC = "Cancer")
wrap_cap <- function(x, width) paste(unlist(lapply(strsplit(x, "\n")[[1]], strwrap, width = width)), collapse = "\n")

## ---- mouse populations: name and colour matched to the CRC populations ----
POP <- tibble::tribble(
  ~key,               ~label,              ~col,       ~crc_match,
  "CD4_Naive",        "CD4 Naive",         "#F4A3A3",  "within CD4 Conventional (IL7R+)",
  "CD4_Activated",    "CD4 Activated",     "#8E1B1B",  "within CD4 Conventional (IL7R+)",
  "CD4_Conventional", "CD4 Conventional",  "#D62728",  "CD4 Conventional",
  "CD4_Th1",          "CD4 Th1",           "#E377C2",  "none",
  "CD4_Th17",         "CD4 Th17",          "#FF7F0E",  "CD4 Th17",
  "CD4_Tfh",          "CD4 TFH",           "#B8860B",  "CD4 TFH",
  "CD4_Treg",         "CD4 Treg",          "#1FB895",  "CD4 Treg",
  "CD8_Effector",     "CD8 Effector",      "#1F77B4",  "CD8 Conventional",
  "CD8_IEL",          "CD8 IEL",           "#00B4FB",  "CD8 Tc17 colour (no direct match)",
  "CD8_TRM",          "CD8 TRM",           "#683B92",  "CD8 Exhausted colour (no direct match)",
  "CD8_Conventional", "CD8 Conventional",  "#A6CEE3",  "CD8 Conventional (light)",
  "gdT17",            "γδ17 T",            "#095F00",  "γδ T",
  "NK",               "NK",                "#8C564B",  "NK",
  "ILC2",             "ILC2",              "#7F7F7F",  "ILC3 (ILC)",
  ## cell-cycle-defined cluster (Mki67/Top2a), 88% Cd3e+, both CD4 and CD8:
  ## shown as "Proliferating T"
  "Cycling",          "Proliferating T",   "#C7C7C7",  "none")
write.csv(POP, file.path(DIR_TAB, "AOMDSS_population_colours_CRC_matched.csv"), row.names = FALSE)

## ---- objeto ----
Tall <- readRDS(AOM_RDS)   # full object
## drop "Unassigned" cells and mice with < 50 cells
Tcells <- subset(Tall, subset = ClusterFull != "Unassigned")
cps <- table(Tcells$sample)
Tcells <- subset(Tcells, subset = sample %in% names(cps)[cps >= MIN_CELLS_PER_SAMPLE])
present <- intersect(POP$key, unique(as.character(Tcells$ClusterFull)))
stopifnot(all(unique(as.character(Tcells$ClusterFull)) %in% POP$key))
pop <- filter(POP, key %in% present)
Tcells$Pop <- factor(pop$label[match(as.character(Tcells$ClusterFull), pop$key)], levels = pop$label)
cols <- setNames(pop$col, pop$label)

## =============== PANEL A ===============
umap_one <- function(g, legend = FALSE) {
  o <- subset(Tcells, subset = Grupo == g)
  DimPlot(o, reduction = "umap", group.by = "Pop", cols = cols, pt.size = 0.15) +
    labs(title = sprintf("%s (n = %s cells)", GROUP_LAB[[g]], format(ncol(o), big.mark = ","))) +
    guides(colour = guide_legend(override.aes = list(size = 2.6), ncol = 1)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0),
          legend.position = if (legend) "right" else "none", legend.text = element_text(size = 8),
          axis.title = element_text(size = 8), axis.text = element_text(size = 7))
}
pA_umap <- (umap_one("Normal") | umap_one("CAC", legend = TRUE))

grupo_of <- Tcells@meta.data |> distinct(sample, Grupo) |> mutate(sample = as.character(sample), Grupo = as.character(Grupo))
comp <- Tcells@meta.data |>
  mutate(sample = as.character(sample)) |>
  count(sample, Pop, name = "N", .drop = FALSE) |>          # one row per mouse x population (0 when there are no cells)
  group_by(sample) |> mutate(Percent = 100 * N / sum(N)) |> ungroup() |>
  left_join(grupo_of, by = "sample") |>
  filter(Grupo %in% c("Normal", "CAC")) |>
  mutate(Group = factor(GROUP_LAB[as.character(Grupo)], levels = c("Healthy", "Cancer")))
n_mice <- comp |> distinct(sample, Group) |> count(Group)

comp_test <- comp |> group_by(Pop) |>
  summarise(p = suppressWarnings(wilcox.test(Percent[Group == "Cancer"], Percent[Group == "Healthy"])$p.value),
            delta = mean(Percent[Group == "Cancer"]) - mean(Percent[Group == "Healthy"]), .groups = "drop") |>
  mutate(fdr = p.adjust(p, "BH"))
write.csv(comp_test, file.path(DIR_TAB, "AOMDSS_composition_Healthy_vs_Cancer_tests.csv"), row.names = FALSE)
## composition change (Cancer - Healthy), horizontal bars sorted by effect
## (no per-mouse dots; colours matched to the CRC populations)
comp_d <- comp |> group_by(Pop) |>
  summarise(delta = mean(Percent[Group == "Cancer"]) - mean(Percent[Group == "Healthy"]),
            sem = sqrt(var(Percent[Group == "Cancer"]) / sum(Group == "Cancer") +
                       var(Percent[Group == "Healthy"]) / sum(Group == "Healthy")), .groups = "drop") |>
  left_join(comp_test |> select(Pop, fdr), by = "Pop") |>
  mutate(stars = case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**", fdr < 0.05 ~ "*", TRUE ~ "ns"),
         Pop = factor(as.character(Pop), levels = as.character(Pop)[order(delta)]),
         err_min = ifelse(delta >= 0, delta, delta - sem), err_max = ifelse(delta >= 0, delta + sem, delta),
         lab_x = ifelse(delta >= 0, err_max, err_min), hj = ifelse(delta >= 0, -0.25, 1.25))
lim_c <- max(abs(c(comp_d$err_min, comp_d$err_max))) * 1.2
pA_bar <- ggplot(comp_d, aes(delta, Pop, fill = Pop)) +
  geom_col(width = 0.75, colour = NA) +
  geom_errorbar(aes(xmin = err_min, xmax = err_max), width = 0.3, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(-lim_c, lim_c)) +
  labs(x = "Δ composition (Cancer − Healthy, %)", y = NULL,
       title = "T/NK/ILC composition change,\nCancer vs Healthy") +
  ## larger text in the composition panel
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(colour = "grey10", size = 12), axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12), plot.title = element_text(size = 12.5, hjust = 0.5, lineheight = 1))

## wider UMAPs and no significance marks on the bars
pA <- (umap_one("Normal") | umap_one("CAC", legend = TRUE) | pA_bar) + plot_layout(widths = c(1.3, 1.3, 1)) +
  plot_annotation(caption = wrap_cap(paste0(
    "Δ = mean % of T/NK/ILC cells per mouse in Cancer minus Healthy ± SEM (Healthy n = ", n_mice$n[n_mice$Group == "Healthy"],
    ", Cancer n = ", n_mice$n[n_mice$Group == "Cancer"], " mice with ≥ ", MIN_CELLS_PER_SAMPLE, " cells). ",
    "Wilcoxon Healthy vs Cancer per population, BH: none significant (min FDR = ", sprintf("%.2f", min(comp_test$fdr)), "). ",
    "Proliferating T = cell-cycle-defined cluster (Mki67+/Top2a+) containing both CD4 and CD8 T cells. ",
    "Population colours matched to the human CRC populations (Fig. 1 CRC palette)."), 210),
    theme = theme(plot.caption = element_text(size = 7, colour = "grey35", hjust = 0)))
ggsave(file.path(OUT_DIR, "SupplFigure7_CD_AOMDSS_UMAP_composition.pdf"), pA, width = 14, height = 3.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "SupplFigure7_CD_AOMDSS_UMAP_composition.png"), pA, width = 14, height = 3.5, dpi = 200, bg = "white")
write.csv(comp_d |> select(Pop, delta, sem, fdr), file.path(DIR_TAB, "AOMDSS_composition_delta.csv"), row.names = FALSE)

## colour legend as a separate file (for assembling the figure)
leg <- ggplot(pop |> mutate(label = factor(label, levels = label)), aes(1, label, colour = label)) + geom_point(size = 3) +
  scale_colour_manual(values = cols, name = NULL) + theme_void(base_size = 9)
ggsave(file.path(OUT_DIR, "SupplFigure7_CD_AOMDSS_legend.pdf"), leg, width = 2, height = 3, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "SupplFigure7_CD_AOMDSS_legend.png"), leg, width = 2, height = 3, dpi = 200, bg = "white")
