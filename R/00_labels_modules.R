###############################################################
## Names of the 9 functional processes (single source for all scripts)
##
## MODULES keeps the internal order used in every figure; MOD_LABELS holds the
## name that is printed on the figures.
###############################################################

MODULES <- c("Glycolysis", "OXPHOS", "Transporters", "Trafficking", "TCA",
             "Tcell_activation", "Cytokine", "Methyltransferase", "Translation")

MOD_LABELS <- setNames(c("Glycolysis",
                         "Oxidative phosphorylation",
                         "Amino acid transporters",
                         "Exocytosis",
                         "TCA cycle",
                         "T cell activation",
                         "Positive regulation of cytokine production",
                         "Methyltransferase activity",
                         "Regulation of protein synthesis machinery"),
                       MODULES)

## wrap long names onto several lines for axes and strips (width = characters per line)
lab_wrap <- function(x, width = 24) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "")
