#!/usr/bin/env bash
set -euo pipefail

# Directory to store pulled SIFs
SIF_DIR="${HOME}/.singularity/mtg-images/"
mkdir -p "$SIF_DIR"
cd "$SIF_DIR"


declare -A IMAGES=(
  [hariszaf-microbetag-nf-0.1.0]="hariszaf/microbetag-nf:0.1.0"
  [hariszaf-carveme-1.6.6]="hariszaf/carveme:1.6.6"
  [hariszaf-flashweave-0.19.2]="hariszaf/flashweave:0.19.2"
  [hariszaf-phenotrex-0.6.0]="hariszaf/phenotrex:0.6.0"
  [biocontainers-prodigal-v1-2.6.3-4-deb_cv1]="biocontainers/prodigal:v1-2.6.3-4-deb_cv1"
  [staphb-hmmer-3.4]="staphb/hmmer:3.4"
)

for name in "${!IMAGES[@]}"; do
  sif="${name}.img"  
  image="docker://${IMAGES[$name]}"
  if [[ -s "$sif" ]]; then
    echo "✅ $sif already exists, skipping."
  else
    echo "⬇️  Pulling $image as $sif..."
    singularity pull --name "$sif" "$image"
  fi
done
