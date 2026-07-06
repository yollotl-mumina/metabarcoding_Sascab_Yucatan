#!/usr/bin/env bash
# =============================================================================
# 01_qiime2_pipeline.sh
# Procesamiento de amplicones ITS (hongos) — Sascaberas de Yucatán
# Reconstrucción del flujo documentado (QIIME 2 amplicon-2023.9)
#
# Requisitos: qiime2 (amplicon-2023.9), FastQC, MultiQC, Cutadapt
# Entrada:   lecturas paired-end demultiplexadas (FASTQ.gz) en raw_reads/
# Salida:    table.qza, taxonomy.qza, taxonomy.qzv y tablas exportadas (.tsv)
# =============================================================================
set -euo pipefail

RAW=raw_reads                 # FASTQ crudos (no incluidos en el repo)
TRIM=trimmed                  # lecturas recortadas con Cutadapt
OUT=qiime2_out
CLASSIFIER=unite_classifier.qza   # clasificador entrenado con UNITE
THREADS=8
mkdir -p "$TRIM" "$OUT" qc

# --- 1. Control de calidad -------------------------------------------------
fastqc "$RAW"/*.fastq.gz -o qc/ -t "$THREADS"
multiqc qc/ -o qc/

# --- 2. Recorte de adaptadores y bases de baja calidad (Cutadapt) ----------
# Q20 forward / Q15 reverse; se remueven los primeros 9 nucleótidos.
# Los archivos siguen el patrón F_KarstCh<n>_S<...>_1_fastq.gz / _2_fastq.gz
for R1 in "$RAW"/*_1_fastq.gz; do
  R2=${R1/_1_fastq.gz/_2_fastq.gz}
  base=$(basename "$R1" | sed 's/_1_fastq.gz//')
  cutadapt \
    -q 20,15 -u 9 -U 9 \
    --minimum-length 50 \
    -o "$TRIM/${base}_1_fastq.gz" -p "$TRIM/${base}_2_fastq.gz" \
    "$R1" "$R2"
done

# --- 3. Importación a QIIME 2 (archivo manifest) ---------------------------
# Se genera el manifest apuntando a las lecturas ya recortadas (carpeta trimmed/).
bash "$(dirname "$0")/make_manifest.sh" "$TRIM"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path manifest.tsv \
  --output-path "$OUT/demux.qza"

# --- 4. Extracción de la región ITS (ITSxpress) ----------------------------
qiime itsxpress trim-pair-output-unmerged \
  --i-per-sample-sequences "$OUT/demux.qza" \
  --p-region ITS1 \
  --p-taxa F \
  --p-threads "$THREADS" \
  --o-trimmed "$OUT/demux_its.qza"

# --- 5. Denoising e inferencia de ASVs (DADA2) -----------------------------
# Truncado forward/reverse a 140 pb; máximo de 6 errores esperados.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$OUT/demux_its.qza" \
  --p-trunc-len-f 140 --p-trunc-len-r 140 \
  --p-max-ee-f 6 --p-max-ee-r 6 \
  --p-n-threads "$THREADS" \
  --o-table "$OUT/table.qza" \
  --o-representative-sequences "$OUT/rep-seqs.qza" \
  --o-denoising-stats "$OUT/denoising-stats.qza"

# --- 6. Asignación taxonómica (clasificador bayesiano UNITE) ---------------
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$OUT/rep-seqs.qza" \
  --o-classification "$OUT/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "$OUT/taxonomy.qza" \
  --o-visualization "$OUT/taxonomy.qzv"

# --- 7. Exportación para análisis en R -------------------------------------
qiime tools export --input-path "$OUT/table.qza"    --output-path exported/
biom convert -i exported/feature-table.biom -o exported/feature-table.tsv --to-tsv
qiime tools export --input-path "$OUT/taxonomy.qza" --output-path exported/

echo "Pipeline QIIME 2 finalizado. Continúa con scripts/02_phyloseq_analysis.R"
