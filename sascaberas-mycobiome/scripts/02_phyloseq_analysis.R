# =============================================================================
# 02_phyloseq_analysis.R
# Análisis ecológico de comunidades fúngicas — Sascaberas de Yucatán
# Reconstrucción del flujo documentado (R / phyloseq)
#
# Entradas: data/exported/feature-table.tsv, data/exported/taxonomy.tsv,
#           metadata/metadata.txt
# Salidas:  figuras en figures/ (abundancia, diversidad alfa, rarefacción, Venn)
# =============================================================================

## ---- 0. Librerías --------------------------------------------------------
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggpubr)
library(MicEco)        # ps_venn()
library(ANCOMBC)
library(DESeq2)
library(ComplexHeatmap)
library(pals)
library(ggsci)
library(microbiome)

set.seed(123)
dir.create("figures", showWarnings = FALSE)

## ---- 1. Importar tabla de ASVs, taxonomía y metadatos --------------------
otu <- read.delim("data/exported/feature-table.tsv", skip = 1,
                  row.names = 1, check.names = FALSE)
otu <- as.matrix(otu)

tax_raw <- read.delim("data/exported/taxonomy.tsv", row.names = 1)
tax <- tax_raw %>%
  rownames_to_column("FeatureID") %>%
  separate(Taxon,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right") %>%
  mutate(across(Kingdom:Species, ~ str_remove(.x, "^\\s*[a-z]__"))) %>%
  column_to_rownames("FeatureID") %>%
  select(Kingdom:Species) %>%
  as.matrix()

meta <- read.delim("metadata/metadata.txt", row.names = 1) %>%
  filter(row_number() != 1)   # descarta la fila de tipos (#Sampleid ...)

ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
               tax_table(tax),
               sample_data(meta))

## ---- 2. Filtrado: conservar solo hongos identificados a nivel de Filo ----
ps <- subset_taxa(ps, Kingdom == "Fungi")
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "unidentified"))
ps   # resumen del objeto phyloseq

## ---- 3. Curvas de rarefacción (Observed y Shannon) -----------------------
rc <- ggrare(ps, step = 200, color = "source", se = FALSE) +
  theme_bw()
ggsave("figures/CurvaObservada.svg", rc, width = 7, height = 4)

## ---- 4. Diversidad alfa (Observed, Chao1, Shannon) -----------------------
alpha <- estimate_richness(ps, measures = c("Observed","Chao1","Shannon"))
alpha <- cbind(alpha, meta[rownames(alpha), ])
write.csv(alpha, "results/alpha_diversity.csv")

p_alpha <- plot_richness(ps, x = "source",
                         measures = c("Observed","Chao1","Shannon")) +
  geom_violin(aes(fill = source), alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw() + scale_fill_npg()
ggsave("figures/DiversidadAlpha.svg", p_alpha, width = 7, height = 3.5)

## ---- 5. Diagramas de Venn (taxa compartidos) -----------------------------
svg("figures/DiagramaVennSascab.svg", width = 5, height = 5)
ps_venn(ps, group = "source",   weight = FALSE)
dev.off()

svg("figures/DiagramaVennColor.svg", width = 5, height = 5)
ps_color <- subset_samples(ps, variety == "color")
ps_venn(ps_color, group = "coloration", weight = FALSE)
dev.off()

## ---- 6. Abundancia relativa por Filo y por Género ------------------------
ps_rel <- transform_sample_counts(ps, function(x) 100 * x / sum(x))

# Filo
p_phy <- plot_bar(ps_rel, fill = "Phylum") +
  facet_wrap(~ variety, scales = "free_x") +
  geom_bar(stat = "identity") + theme_bw() +
  labs(y = "%") + scale_fill_manual(values = as.vector(cols25()))
ggsave("figures/PhylumAbd.svg", p_phy, width = 8, height = 4)

# Género (colapsando taxa < 1 %)
glom <- tax_glom(ps_rel, taxrank = "Genus")
df_g <- psmelt(glom) %>%
  mutate(Genus = ifelse(Abundance < 1, "< 1%", Genus))
p_gen <- ggplot(df_g, aes(Sample, Abundance, fill = Genus)) +
  geom_bar(stat = "identity") + theme_bw() + labs(y = "%") +
  theme(legend.text = element_text(size = 6))
ggsave("figures/AbundanciaGenero.svg", p_gen, width = 10, height = 5)

## ---- 7. Diversidad beta y PERMANOVA --------------------------------------
ps_hell <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
dist_bc <- phyloseq::distance(ps_hell, method = "bray")
ord <- ordinate(ps_hell, method = "PCoA", distance = dist_bc)
plot_ordination(ps_hell, ord, color = "source") + theme_bw()

adonis2(dist_bc ~ source, data = as(sample_data(ps), "data.frame"))

## ---- 8. Abundancia diferencial (ANCOM-BC / DESeq2) -----------------------
# ANCOM-BC entre sustratos
anc <- ancombc2(data = ps, fix_formula = "source", p_adj_method = "BH")
head(anc$res)

sessionInfo()
