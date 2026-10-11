###############################################################
## Change in the 9 modules, CD4 versus CD8 (CRC)
##
## Figure  : Fig. 8D
## Input   : tables/CRC_module_scores_per_patient.csv (script 03)
## Output  : Figure8_D_CRC_CD4_vs_CD8_delta.{png,pdf}
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

OUT_DIR <- DIR_F8
PB_CSV  <- file.path(DIR_TAB, "CRC_module_scores_per_patient.csv")

comp_cols <- c(CD4 = "#376199", CD8 = "#B37CB8")   # same colours as the rest of the figure

## ---- per-donor change (only donors with both specimens) ----
deltas <- read.csv(PB_CSV) |>
  filter(sig %in% MODULES) |>
  select(donor, lineage, spec, sig, score) |>
  pivot_wider(names_from = spec, values_from = score) |>
  filter(!is.na(Normal), !is.na(Tumor)) |>
  mutate(delta = Tumor - Normal)

summ <- deltas |> group_by(lineage, sig) |>
  summarise(n = n(), delta_mean = mean(delta), sem = sd(delta) / sqrt(n()), .groups = "drop")

## ---- CD4 vs CD8 per module (paired by donor) ----
cmp <- deltas |> select(donor, lineage, sig, delta) |>
  pivot_wider(names_from = lineage, values_from = delta) |>
  filter(!is.na(CD4), !is.na(CD8)) |>
  group_by(sig) |>
  summarise(n = n(), p = suppressWarnings(wilcox.test(CD4, CD8, paired = TRUE)$p.value), .groups = "drop") |>
  mutate(fdr = p.adjust(p, "BH"))
write.csv(cmp, file.path(DIR_TAB, "CRC_CD4_vs_CD8_comparison.csv"), row.names = FALSE)

## ---- figure ----
fdr_stars <- function(f) case_when(f < 0.001 ~ "***", f < 0.01 ~ "**", f < 0.05 ~ "*", TRUE ~ "ns")

plot_df <- summ |>
  mutate(Label = factor(MOD_LABELS[sig], levels = rev(MOD_LABELS)),
         err_min = ifelse(delta_mean >= 0, delta_mean, delta_mean - sem),
         err_max = ifelse(delta_mean >= 0, delta_mean + sem, delta_mean),
         star_y  = as.numeric(Label) + ifelse(lineage == "CD4", 0.2, -0.2))

bracket_df <- plot_df |>
  select(Label, sig, lineage, err_max, star_y) |>
  pivot_wider(id_cols = c(Label, sig), names_from = lineage, values_from = c(err_max, star_y)) |>
  left_join(cmp |> select(sig, fdr), by = "sig") |>
  mutate(Stars = fdr_stars(fdr),
         tick_start_CD4 = err_max_CD4 + 0.008,
         tick_start_CD8 = err_max_CD8 + 0.008,
         bracket_x = pmax(tick_start_CD4, tick_start_CD8) + 0.015,
         mid_y = (star_y_CD4 + star_y_CD8) / 2)

n_lab <- paste0("n = ", unique(summ$n), " paired patients")

pC <- ggplot(plot_df, aes(delta_mean, Label, fill = lineage, group = lineage)) +
  geom_col(position = position_dodge(width = 0.8, reverse = TRUE), width = 0.7,
           colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(xmin = err_min, xmax = err_max),
                position = position_dodge(width = 0.8, reverse = TRUE), width = 0.2, linewidth = 0.4) +
  geom_segment(data = bracket_df, inherit.aes = FALSE, linewidth = 0.4,
               aes(x = tick_start_CD4, xend = bracket_x, y = star_y_CD4, yend = star_y_CD4)) +
  geom_segment(data = bracket_df, inherit.aes = FALSE, linewidth = 0.4,
               aes(x = tick_start_CD8, xend = bracket_x, y = star_y_CD8, yend = star_y_CD8)) +
  geom_segment(data = bracket_df, inherit.aes = FALSE, linewidth = 0.4,
               aes(x = bracket_x, xend = bracket_x, y = star_y_CD4, yend = star_y_CD8)) +
  geom_text(data = bracket_df, inherit.aes = FALSE, aes(x = bracket_x + 0.005, y = mid_y, label = Stars),
            size = 3.6, fontface = "bold", hjust = 0, vjust = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  scale_fill_manual(values = comp_cols) +
  scale_y_discrete(labels = lab_wrap) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "CD4 vs CD8 functional program changes (CRC)",
       x = expression(Delta ~ "score (Cancer - Healthy)"), y = NULL, fill = NULL,
       caption = paste0("Bars: mean Δ ± SEM per patient (", n_lab, "). ",
                        "Bracket: CD4 vs CD8, paired Wilcoxon, BH across the 9 modules.\n",
                        "*** FDR < 0.001, ** FDR < 0.01, * FDR < 0.05, ns = not significant.")) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.title.position = "plot",
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_text(colour = "black"),
        legend.position = "right",
        plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0, lineheight = 1.1),
        plot.caption.position = "plot")

ggsave(file.path(OUT_DIR, "Figure8_D_CRC_CD4_vs_CD8_delta.pdf"), pC, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure8_D_CRC_CD4_vs_CD8_delta.png"), pC, width = 6.5, height = 5, dpi = 200, bg = "white")
print(cmp |> mutate(module = MOD_LABELS[sig]) |> select(module, n, p, fdr), digits = 3)
cat("OK\n")
