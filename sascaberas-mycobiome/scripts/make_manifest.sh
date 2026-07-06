#!/usr/bin/env bash
# Genera manifest.tsv con rutas ABSOLUTAS para QIIME 2
# (formato PairedEndFastqManifestPhred33V2). Ejecutar desde la raíz del repo.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
READS_DIR="${1:-$ROOT/raw_reads}"     # opcional: carpeta de lecturas (por defecto raw_reads)
OUT="$ROOT/manifest.tsv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$OUT"
for n in 1 2 3 4 5 6 7 8; do
  r1=$(ls "$READS_DIR"/F_KarstCh${n}_*_1_fastq.gz)
  r2=$(ls "$READS_DIR"/F_KarstCh${n}_*_2_fastq.gz)
  printf "KarstCh%s\t%s\t%s\n" "$n" "$r1" "$r2" >> "$OUT"
done
echo "manifest.tsv generado en $OUT (lecturas: $READS_DIR)"
