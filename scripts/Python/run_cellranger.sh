#!/bin/bash

# Configuración de rutas relativas
SCRIPT_DIR=$(dirname "$0")
FASTQ_DIR="$SCRIPT_DIR/YPUR_SAMPLES_DIRECTORY"
TRANSCRIPTOME_DIR="$SCRIPT_DIR/refdata-gex-GRCh38-2024-A"
OUTPUT_ID="OUTPUT_SAMPLE"

# Verificación de la existencia de directorios
if [ ! -d "$FASTQ_DIR" ]; then
  echo "Error: El directorio de FASTQ '$FASTQ_DIR' no existe."
  exit 1
fi

if [ ! -d "$TRANSCRIPTOME_DIR" ]; then
  echo "Error: El directorio del transcriptoma '$TRANSCRIPTOME_DIR' no existe."
  exit 1
fi

# Ejecución de cellranger count
cellranger count \
  --id="$OUTPUT_ID" \
  --fastqs="$FASTQ_DIR" \
  --sample=SAMPLE_ID \
  --transcriptome="$TRANSCRIPTOME_DIR" \
  --create-bam=true

echo "Procesamiento de cellranger completado."
