# Figure code — González-Romero *et al.*

R code that reproduces **Figure 8** and **Supplementary Figures 7, 8 and 9** of

> González-Romero D., Spinelli L., Fajardo P., Díaz-Mora E., Gómez de Cedrón M.,
> Sanz-Ezquerro J.J., Cantrell D.A., Cuenda A. **Differential proteomic and metabolic
> remodelling of CD4⁺ and CD8⁺ T cells in mesenteric lymph nodes during
> colitis-associated colorectal cancer.** *Cell Death & Disease*, 2026. doi:[DOI]
>
> Preprint: bioRxiv, doi:[10.64898/2025.12.26.696593](https://doi.org/10.64898/2025.12.26.696593)

These figures are a reanalysis of three public single-cell RNA-seq datasets (human
colorectal cancer, human ulcerative colitis and the mouse AOM/DSS model of
colitis-associated cancer) together with the mesenteric lymph node (mLN) proteomics
generated in this study.

Each script writes the file that was placed in the figure; panel letters, cropping and
diagrams were added afterwards in PowerPoint.

---

## Quick start

```bash
git clone https://github.com/CuendaLab/Gonzalez-Romero-2026-Figures.git
cd Gonzalez-Romero-2026-Figures
```

1. Get the data (see [Data](#data)) and build the three Seurat objects with the
   pipelines in [`R/build/`](R/build) — this is done once.
2. Edit the four paths at the top of [`R/00_paths.R`](R/00_paths.R).
3. Run everything:

```r
source("run_all.R")      # or: Rscript run_all.R
```

Figures are written to `figures/` (one folder per figure), supporting tables to
`tables/`, and the R and package versions of the run to `tables/sessionInfo.txt`.

## Repository layout

```
├── run_all.R              runs every script in the right order
├── R/
│   ├── 00_paths.R         the four data paths (the only thing to edit)
│   ├── 00_labels_*.R      process and dataset names used in the figures
│   ├── 01 … 11c_*.R       one script per figure panel
│   └── build/             construction of the Seurat objects (run once)
├── gene_signatures/       human signature (412 genes) and mouse orthologs (409)
├── input_tables/          per-module change of each dataset (input of Suppl. Fig. 9B)
├── figures/               output, one folder per figure
└── tables/                output, supporting tables
```

## Which script makes which panel

| Figure | Panel | Script |
|---|---|---|
| Fig. 8 | A (UMAP + legend), B (composition change) | [`R/01_fig8AB_CRC_UMAP_composition.R`](R/01_fig8AB_CRC_UMAP_composition.R) |
| Fig. 8 | C (per-patient scores, CD4) | [`R/03_fig8C_suppl8B_CRC_per_patient_scores.R`](R/03_fig8C_suppl8B_CRC_per_patient_scores.R) |
| Fig. 8 | D (change, CD4 vs CD8) | [`R/04_fig8D_CRC_CD4_vs_CD8_delta.R`](R/04_fig8D_CRC_CD4_vs_CD8_delta.R) |
| Fig. 8 | E (gene-level heatmaps, CD4) | [`R/05_fig8E_suppl8A_CRC_CD4_allgenes_grid.R`](R/05_fig8E_suppl8A_CRC_CD4_allgenes_grid.R) |
| Fig. 8 | F (direction of change per gene) | [`R/07a_gene_level_concordance_table.R`](R/07a_gene_level_concordance_table.R) → [`R/07b_fig8F_direction_heatmap.R`](R/07b_fig8F_direction_heatmap.R) |
| Suppl. 7 | C (UMAP), D (composition change) | [`R/08_suppl7CD_AOM_UMAP_composition.R`](R/08_suppl7CD_AOM_UMAP_composition.R) |
| Suppl. 7 | E (gene-level heatmaps, CD4) | [`R/09_suppl7E_AOM_CD4_allgenes_grid.R`](R/09_suppl7E_AOM_CD4_allgenes_grid.R) |
| Suppl. 8 | A (gene-level heatmaps, CD4) | same file as Fig. 8E (other processes cropped from it) |
| Suppl. 8 | B (per-patient scores, CD8) | [`R/03_fig8C_suppl8B_CRC_per_patient_scores.R`](R/03_fig8C_suppl8B_CRC_per_patient_scores.R) |
| Suppl. 8 | C (gene-level heatmaps, CD8) | [`R/06_suppl8C_CRC_CD8_allgenes_4columns.R`](R/06_suppl8C_CRC_CD8_allgenes_4columns.R) |
| Suppl. 9 | A (gene-level heatmaps, colitis) | [`R/10_suppl9A_UC_gene_heatmaps.R`](R/10_suppl9A_UC_gene_heatmaps.R) |
| Suppl. 9 | B (cross-dataset summary) | [`R/11a`](R/11a_proteomics_mLN_module_delta.R) → [`R/11b`](R/11b_crossdataset_table.R) → [`R/11c_suppl9B_crossdataset_dotmatrix.R`](R/11c_suppl9B_crossdataset_dotmatrix.R) |

Panels A–B of Suppl. Fig. 7 (model scheme and flow cytometry) and the diagrams of
Figure 8 are not generated here. The colour bar shown separately in some figures is
cropped from the gene-level heatmaps themselves.

Script dependencies (handled by `run_all.R`): `04` needs the table from `03`, `07b` the
one from `07a`, `11b` the one from `11a`, and `11c` the one from `11b`.

## Data

| Dataset | Accession | What is used |
|---|---|---|
| Human colorectal cancer (Pelka *et al.*, *Cell* 2021) | Broad Single Cell Portal [SCP1162](https://singlecell.broadinstitute.org/single_cell/study/SCP1162) / GEO [GSE178341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341) | T cells with the authors' annotation; 36 patients with paired adjacent normal colon and tumour |
| Mouse AOM/DSS | GEO [GSE264607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264607) | colonic T/NK/ILC cells, Healthy and CAC groups |
| Ulcerative colitis (Boland *et al.*) | GEO [GSE125527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125527) | tissue T cells |
| mLN proteomics | generated in this study, see Data availability of the paper | copy numbers per cell (normalised to histone content); CD4 and CD8 × control/colitis/cancer × 3 replicates |

The Seurat objects are not distributed here because of their size. They are built from
the public matrices with the three pipelines in [`R/build/`](R/build), documented in
[`R/build/README.md`](R/build/README.md). `input_tables/` holds the per-module change
tables of each dataset, which feed the cross-dataset summary of Suppl. Fig. 9B.

## Gene signatures

`gene_signatures/cd4_gene_level_prognosis.csv` (human, 412 genes) and
`cd4_gene_level_prognosis_mouse.csv` (409 one-to-one mouse orthologs) assign each gene
to one or more of the nine functional processes, in the fixed order used throughout:
glycolysis (GO:0006096), oxidative phosphorylation (GO:0006119), amino acid transporters
(GO:0006865), exocytosis (GO:0006887), TCA cycle, T cell activation (GO:0042110),
positive regulation of cytokine production (GO:0001819), methyltransferase activity
(GO:0008168) and regulation of protein synthesis machinery.

## Method in brief

- **Module score** — per-gene Z-score across all CD4 and CD8 T cells together, mean of
  the genes of the process per cell, mean of its cells per donor. In the proteomics,
  per-protein Z-score of the log10 copy number and mean per sample.
- **Comparisons** — paired Wilcoxon signed-rank test in colorectal cancer (n = 36
  pairs), unpaired Wilcoxon rank-sum test in AOM/DSS and colitis, Welch's *t*-test in
  the proteomics (n = 3 vs 3). Benjamini–Hochberg correction across the nine processes
  within each dataset and lineage.
- **Gene-level heatmaps** — mean per donor and condition, per-gene Z-score across donors
  capped at ±2, genes ordered by the difference between group means; genes not detected
  in that lineage are excluded.
- **Direction of change (Fig. 8F)** — sign of the change of each gene in each dataset;
  the concordance expected by chance is p₁p₂ + (1−p₁)(1−p₂), from the proportion of
  increased genes on each platform.

## Requirements

R ≥ 4.4 (developed on 4.6.0) with **Seurat** 5.5.0, **ComplexHeatmap** 2.28.0,
**circlize** 0.4.18, **patchwork** 1.3.2, **cowplot**, **readxl**, **Matrix** and the
tidyverse (**dplyr** 1.2.1, **tidyr** 1.3.2, **ggplot2** 4.0.3). The `R/build/`
pipelines also need **harmony** and **scDblFinder**.

```r
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork",
                   "cowplot", "readxl", "circlize"))
BiocManager::install(c("ComplexHeatmap", "scDblFinder"))
```

Arial is embedded through `cairo_pdf`; on a system without Arial, change the `FONT`
constant at the top of the heatmap scripts.

## Citation

If you use this code, please cite the paper above. [Add the Zenodo DOI of the
repository here once it is archived.]

## License

Released under the MIT License — see [LICENSE](LICENSE). You may use, modify and
redistribute this code, including for commercial purposes, as long as the copyright
notice is kept. The code is provided as is, without warranty.

The license covers the code only. The public datasets keep the terms of their own
repositories, and the mLN proteomics is governed by the Data availability statement of
the paper.
