# singularity pull docker://hariszaf/microbetag-nf:0.1.0 
# singularity pull docker://hariszaf/carveme:1.6.6
# singularity pull docker://hariszaf/flashweave:0.19.2
# singularity pull docker://hariszaf/phenotrex:0.6.0 
# singularity pull docker://biocontainers/prodigal:v1-2.6.3-4-deb_cv1
# singularity pull docker://staphb/hmmer:3.4
#!/usr/bin/env bash
set -euo pipefail

# Directory to store pulled SIFs
SIF_DIR="${HOME}/singularity_images"
mkdir -p "$SIF_DIR"
cd "$SIF_DIR"

declare -A IMAGES=(
  [microbetag-nf]="hariszaf/microbetag-nf:0.1.0"
  [carveme]="hariszaf/carveme:1.6.6"
  [flashweave]="hariszaf/flashweave:0.19.2"
  [phenotrex]="hariszaf/phenotrex:0.6.0"
  [prodigal]="biocontainers/prodigal:v1-2.6.3-4-deb_cv1"
  [hmmer]="staphb/hmmer:3.4"
)

for name in "${!IMAGES[@]}"; do
  sif="${name}.sif"
  image="docker://${IMAGES[$name]}"
  if [[ -s "$sif" ]]; then
    echo "✅ $sif already exists, skipping."
  else
    echo "⬇️  Pulling $image..."
    singularity pull "$sif" "$image"
  fi
done