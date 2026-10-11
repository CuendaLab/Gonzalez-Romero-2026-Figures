###############################################################
## Gene-by-gene map of the direction of change in the three datasets
##
## Figure  : Fig. 8F
## Input   : tables/Concordance_gene_level_table.csv (script 07a)
## Output  : Figure8_F_direction_heatmap.{png,pdf} (+ summary tables)
###############################################################

source("R/00_paths.R")
source("R/00_labels_modules.R")

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

OUT_DIR <- DIR_TAB
genes <- read.csv(file.path(OUT_DIR, "Concordance_gene_level_table.csv"))


COMPS <- c("AOM-DSS scRNA (mouse) vs CRC Pelka scRNA (human)",
           "AOM-DSS scRNA (mouse) vs mLN proteomics (mouse)",
           "mLN proteomics (mouse) vs CRC Pelka scRNA (human)")
## dataset labels
SHORT <- c("CAC colon scRNA vs\nCRC colon scRNA", "CAC colon scRNA vs\nCAC mLN proteomic",
           "CAC mLN proteomic vs\nCRC colon scRNA")
COL_UP <- "#C23B22"; COL_DOWN <- "#2C7BB6"

## ---------------- (2) mapa gen a gen ----------------
WRAP <- 30
## append the counts to the last line of the name if they fit; otherwise move them down whole
lab_num <- function(name, num) {
  l <- strwrap(name, WRAP)
  if (is.null(num)) return(paste(l, collapse = "\n"))
  n <- length(l)
  if (nchar(l[n]) + 1 + nchar(num) <= WRAP) l[n] <- paste(l[n], num) else l <- c(l, num)
  paste(l, collapse = "\n")
}
first_mod <- function(m) { x <- strsplit(m, ";")[[1]]; MODULES[MODULES %in% x][1] }
wide <- bind_rows(
  genes |> filter(comparison == COMPS[1]) |> transmute(lineage, gene, human, modules, ds = "Mouse scRNA", d = delta_1),
  genes |> filter(comparison == COMPS[1]) |> transmute(lineage, gene, human, modules, ds = "CRC scRNA", d = delta_2),
  genes |> filter(comparison == COMPS[3]) |> transmute(lineage, gene, human, modules, ds = "Mouse proteomics", d = delta_1)) |>
  distinct(lineage, gene, ds, .keep_all = TRUE) |>
  pivot_wider(names_from = ds, values_from = d) |>
  filter(!is.na(`Mouse scRNA`), !is.na(`Mouse proteomics`), !is.na(`CRC scRNA`)) |>
  mutate(process = factor(MOD_LABELS[sapply(modules, first_mod)], levels = MOD_LABELS),
         n_up = (`Mouse scRNA` > 0) + (`Mouse proteomics` > 0) + (`CRC scRNA` > 0))

## column order: mLN proteomics, mouse scRNA, human scRNA
DS <- c("Mouse proteomics", "Mouse scRNA", "CRC scRNA")
heat <- function(lin, show_n = FALSE) {
  w <- wide |> filter(lineage == lin) |>
    arrange(process, desc(n_up), desc(`Mouse proteomics` > 0), desc(`Mouse scRNA` > 0), human) |>
    group_by(process) |> mutate(row = row_number()) |> ungroup()
  strip_lab <- w |> group_by(process) |>
    ## process name and counts on the same line unless they do not fit:
    ## then the whole block of counts moves to the next line (it is never split)
    summarise(lab = lab_num(as.character(process[1]),
                            if (show_n) sprintf("(\u2191%d \u2193%d / %d)", sum(n_up == 3), sum(n_up == 0), n())),
              .groups = "drop")
  w <- w |> mutate(proc_lab = factor(strip_lab$lab[match(process, strip_lab$process)], levels = strip_lab$lab))
  long <- w |> pivot_longer(all_of(DS), names_to = "ds", values_to = "d") |>
    mutate(ds = factor(ds, DS, labels = c("CAC mLN proteomic\n(mouse)", "CAC colon scRNA\n(mouse)", "CRC colon scRNA\n(human)")),
           dir = ifelse(d > 0, "Up", "Down"))
  ggplot(long, aes(ds, row)) +
    geom_tile(aes(fill = dir), width = 0.92, colour = NA) +
    facet_grid(proc_lab ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(Up = COL_UP, Down = COL_DOWN), name = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    ## fixed margin per process so the TCA cycle label (6 genes) fits; with the
    ## counts the label takes one extra line and needs more room
    scale_y_reverse(expand = expansion(add = if (show_n) 4 else 2.5)) +
    labs(x = NULL, y = NULL, title = lin) +
    theme_minimal(base_size = 26) +   # larger text
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = if (show_n) 22 else 24, lineheight = 0.95),
          strip.placement = "outside", strip.clip = "off", panel.spacing.y = unit(3, "pt"),
          axis.text.x = element_text(size = 24, colour = "grey15", lineheight = 0.9, angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(face = "bold", hjust = 0.5))
}
## two versions: without counts and with the gene counts of each process
draw2 <- function(show_n, stub) {
  n_tot <- wide |> filter(lineage == "CD4") |> nrow()
  p2 <- (heat("CD4", show_n) | heat("CD8", show_n)) +
    plot_layout(guides = "collect") +
    plot_annotation(title = if (show_n)
                      sprintf("Direction of change per gene (%d 9-process genes detected in all 3 datasets)", n_tot)
                    else "Direction of change per gene (9-process genes detected in all 3 datasets)",
                    caption = paste(c(strwrap("Red = up, blue = down (CAC mLN proteomic, Cancer vs day 0; CAC colon scRNA, Cancer vs Healthy; CRC colon scRNA, Cancer vs Healthy).", 110),
                                      strwrap(if (show_n) "Genes in several processes are shown under the first one; in brackets, genes up in all three datasets, down in all three, and genes shown for that process."
                                              else "Genes in several processes are shown under the first one.", 110)), collapse = "\n"),
                    theme = theme(plot.title = element_text(face = "bold", size = 30),
                                  plot.caption = element_text(size = 21, colour = "grey30", hjust = 0, lineheight = 1.1))) &
    ## larger legend
    theme(legend.position = "right", legend.text = element_text(size = 28),
          legend.key.size = unit(16, "mm"), legend.key.spacing.y = unit(4, "mm"))
  ggsave(paste0(stub, ".pdf"), p2, width = 17, height = 22, device = cairo_pdf)
  ggsave(paste0(stub, ".png"), p2, width = 17, height = 22, dpi = 200, bg = "white")
}
draw2(TRUE, file.path(DIR_F8, "Figure8_F_direction_heatmap"))

all3 <- wide |> group_by(lineage) |>
  summarise(n = n(), up_all3 = sum(n_up == 3), down_all3 = sum(n_up == 0), mixed = sum(n_up %in% 1:2),
            pct_same_all3 = 100 * mean(n_up %in% c(0, 3)), .groups = "drop")
write.csv(all3, file.path(DIR_TAB, "Direction_all3_summary.csv"), row.names = FALSE)
write.csv(wide, file.path(DIR_TAB, "Direction_per_gene_3datasets.csv"), row.names = FALSE)
print(all3)
