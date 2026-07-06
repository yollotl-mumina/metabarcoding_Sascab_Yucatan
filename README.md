# Diversidad de comunidades fúngicas en sascaberas de Yucatán

Caracterización de la micobiota de **sascaberas** (depósitos de *sascab*, material
calcáreo blando del karst yucateco) y rizósferas asociadas de la Península de
Yucatán, mediante **metabarcoding de la región ITS** del ADN ribosomal.
El repositorio contiene el flujo bioinformático (QIIME 2), el análisis ecológico
(R / phyloseq), los datos procesados, las figuras y el informe científico.

**Autora:** Daniela Jacqueline Zúñiga Jiménez (Licenciatura en Agroecología, UADY)
**Asesora:** Dra. Miriam Montserrat Ferrer Ortega (FMVZ, UADY)

---

## Resumen del estudio

Se procesaron **8 muestras ambientales** (sascab y rizósfera; morfotipos de
coloración blanca, naranja, verde, roja y beige). Tras el control de calidad,
recorte, extracción de ITS, inferencia de ASVs (DADA2) y asignación taxonómica
(UNITE), se recuperaron **918 ASVs**, todas del reino **Fungi**.

Resultados principales:

- **Ascomycota** fue el filo dominante en todas las muestras, seguido de
  Basidiomycota y Mortierellomycota.
- Los géneros más abundantes fueron ***Acremonium*** (~21 %) y ***Fusarium*** (~18 %).
- Las muestras de **sascab** presentaron mayor riqueza (≈238 ASVs) pero menor
  equidad (Shannon ≈3.56) que las **rizósferas** (≈200 ASVs; Shannon ≈4.29).
- Baja proporción de taxa compartidos entre grupos → fuerte especificidad de hábitat
  (solo ~5 % de ASVs comunes a los cuatro morfotipos de coloración).

---

## Estructura del repositorio

```
sascaberas-mycobiome/
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml              # entorno conda (FastQC, MultiQC, Cutadapt, ITSxpress)
├── manifest.tsv                  # rutas de las lecturas para QIIME 2 (regenerable)
├── raw_reads/                    # 16 FASTQ.gz (paired-end, 8 muestras)
├── metadata/
│   └── metadata.txt              # mapa de muestras (sustrato, variedad, coloración)
├── data/
│   ├── read_counts.tsv           # conteo de lecturas crudas por muestra
│   ├── qiime2/                   # artefactos QIIME 2
│   │   ├── table.qza
│   │   ├── taxonomy.qza
│   │   └── taxonomy.qzv
│   ├── exported/                 # tablas exportadas para R
│   │   ├── feature-table.tsv     # tabla de ASVs (918 ASVs × 8 muestras)
│   │   └── Generos.txt           # abundancias por género (formato largo)
│   └── DatosPhylo.RData          # objeto phyloseq serializado
├── scripts/
│   ├── make_manifest.sh          # genera manifest.tsv con rutas absolutas
│   ├── install_r_packages.R      # instala dependencias de R y registra versiones
│   ├── 01_qiime2_pipeline.sh     # QIIME 2: QC → Cutadapt → ITSxpress → DADA2 → UNITE
│   └── 02_phyloseq_analysis.R    # R: diversidad alfa/beta, rarefacción, Venn, abundancias
├── figures/                      # figuras generadas (SVG)
└── results/
    └── Informe_Sascaberas.docx   # informe científico completo
```

---

## Datos

| Archivo | Descripción |
|---|---|
| `metadata/metadata.txt` | Identificador de muestra y factores: `source` (sascab / rhizosphere), `variety` (white / color), `coloration`. |
| `data/qiime2/table.qza` | Tabla de frecuencias de ASVs (QIIME 2). |
| `data/qiime2/taxonomy.qza` / `.qzv` | Asignación taxonómica (clasificador UNITE) y su visualización. |
| `data/exported/feature-table.tsv` | Tabla de ASVs exportada (entrada de R). |
| `data/exported/Generos.txt` | Abundancias relativas por género y muestra. |
| `data/DatosPhylo.RData` | Objeto `phyloseq` listo para reanalizar en R. |

### Lecturas crudas

Las lecturas ITS paired-end (16 archivos `FASTQ.gz`, ~30 MB) están en `raw_reads/`.
Todos los archivos pasaron la verificación de integridad `gzip` y los conteos de
R1 y R2 coinciden en todas las muestras (longitud de lectura ≈ 291–292 pb).

| Muestra | Lecturas crudas | Muestra | Lecturas crudas |
|---|---:|---|---:|
| KarstCh1 | 5 794 | KarstCh5 | 58 945 |
| KarstCh2 | 64 593 | KarstCh6 | 6 573 |
| KarstCh3 | 13 084 | KarstCh7 | 7 151 |
| KarstCh4 | 6 266 | KarstCh8 | 5 762 |

> Para una publicación pública se recomienda depositar las lecturas en NCBI SRA / ENA
> y enlazar el número de acceso; para repos grandes puede usarse Git LFS.

---

## Reproducción

### 0. Entorno

```bash
# Herramientas de QC/recorte (FastQC, MultiQC, Cutadapt, ITSxpress)
conda env create -f environment.yml && conda activate sascaberas-qc

# QIIME 2 amplicon-2023.9 (distribución oficial, entorno aparte)
conda env create -n qiime2-amplicon-2023.9 \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2023.9-py38-linux-conda.yml

# Paquetes de R
Rscript scripts/install_r_packages.R   # genera SESSION_INFO.txt con las versiones
```

### 1. Procesamiento (QIIME 2)

```bash
# Requiere: qiime2 amplicon-2023.9, FastQC, MultiQC, Cutadapt
bash scripts/make_manifest.sh      # genera manifest.tsv con rutas absolutas
bash scripts/01_qiime2_pipeline.sh # necesita el clasificador UNITE (unite_classifier.qza)
```

Parámetros clave: Cutadapt Q20 (fwd) / Q15 (rev) y remoción de 9 nt; extracción de
ITS1 con **ITSxpress**; **DADA2** con truncado a 140 pb y `max-ee = 6`; asignación
con **classify-sklearn** entrenado sobre **UNITE**.

### 2. Análisis ecológico (R)

```r
# Requiere R ≥ 4.2 y los paquetes indicados abajo
source("scripts/02_phyloseq_analysis.R")
```

**Paquetes de R:** phyloseq, tidyverse, vegan, ggpubr, MicEco, ANCOMBC, DESeq2,
ComplexHeatmap, pals, ggsci, microbiome.

---

## Referencias principales

- Bolyen, E., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology, 37*, 852–857. https://doi.org/10.1038/s41587-019-0209-9
- Callahan, B. J., et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. https://doi.org/10.1038/nmeth.3869
- McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE, 8*(4), e61217. https://doi.org/10.1371/journal.pone.0061217
- Nilsson, R. H., et al. (2019). The UNITE database for molecular identification of fungi. *Nucleic Acids Research, 47*(D1), D259–D264. https://doi.org/10.1093/nar/gky1022
- Rivers, A. R., et al. (2018). ITSxpress: Software to rapidly trim ITS sequences. *F1000Research, 7*, 1418. https://doi.org/10.12688/f1000research.15704.1
- Man, B., et al. (2018). Diversity of fungal communities in Heshang Cave. *Frontiers in Microbiology, 9*, 1400. https://doi.org/10.3389/fmicb.2018.01400

La lista completa está en el informe (`results/Informe_Sascaberas.docx`).

---

## Cómo citar

> Zúñiga Jiménez, D. J., & Ferrer Ortega, M. M. (2025). *Diversidad de comunidades
> fúngicas en sascaberas de la Península de Yucatán caracterizada mediante
> metabarcoding de la región ITS*. Universidad Autónoma de Yucatán.

## License

This repository combines code and research materials, released under separate licenses:

- **Code** (QIIME 2 and R scripts, in `scripts/`): **MIT License** (see `LICENSE`). Free to reuse, modify, and redistribute with attribution.
- **Data, figures, and report** (`raw_reads/`, `data/`, `metadata/`, `figures/`, `results/`): **Creative Commons Attribution 4.0 International (CC-BY 4.0)**. Use and adaptation are permitted provided appropriate credit is given.

When reusing any part of this work, please cite it as indicated in the *How to cite* section above.
