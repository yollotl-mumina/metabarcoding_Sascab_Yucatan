# =============================================================================
# install_r_packages.R — instala las dependencias del análisis en R
# Uso:  Rscript scripts/install_r_packages.R
# Recomendado: R >= 4.2 con Bioconductor >= 3.16
# =============================================================================

# --- Bioconductor ----------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_pkgs <- c("phyloseq", "ANCOMBC", "DESeq2", "ComplexHeatmap", "microbiome")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

# --- CRAN ------------------------------------------------------------------
cran_pkgs <- c("tidyverse", "vegan", "ggpubr", "pals", "ggsci", "remotes")
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install))
  install.packages(to_install, repos = "https://cloud.r-project.org")

# --- GitHub ----------------------------------------------------------------
# MicEco (ps_venn) puede instalarse desde GitHub si no está en CRAN:
if (!requireNamespace("MicEco", quietly = TRUE))
  remotes::install_github("Russel88/MicEco")

# --- Registro de versiones -------------------------------------------------
writeLines(capture.output(sessionInfo()), "SESSION_INFO.txt")
cat("Dependencias instaladas. Versiones registradas en SESSION_INFO.txt\n")
