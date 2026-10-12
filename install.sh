#!/usr/bin/env bash
#
# Install microbetag and its dependencies as conda environments.
#
# Every environment is described in envs/*.yml; this script only creates them and runs the few
# steps conda can't: installing FlashWeave inside Julia and, optionally, downloading the kofam database.
#
#   bash install.sh                 # microbetag environment only
#   bash install.sh --phenotrex     # + the mtg-phenotrex environment
#   bash install.sh --all --kofam   # everything

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENVS_DIR="$SCRIPT_DIR/envs"

KOFAM=false
PHENOTREX=false
MODELSEED=false
DNNGIOR=false
FLASHWEAVE=true

usage() {
  cat <<EOF
Installation script for microbetag.

Usage: bash install.sh [options]

Always creates (or updates) the 'microbetag' conda environment from envs/microbetag.yml.

Options:
  -k, --kofam          Download the kofam database into ext_data/kofam_database (~1.5 GB)
  -p, --phenotrex      Create the mtg-phenotrex environment (trait prediction with phenotrex)
  -m, --modelseed      Create the mtg-modelseed environment (GEM reconstruction with ModelSEEDpy)
  -d, --dnngior        Create the mtg-dnngior environment (gap-filling with DNNGIOR)
  -a, --all            Same as --phenotrex --modelseed --dnngior
      --no-flashweave  Skip installing FlashWeave in Julia (only needed to infer networks)
  -h, --help           Show this help message

Requires conda (or mamba). Join https://matrix.to/#/#microbetagcommunity:matrix.org if you run into trouble.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -k|--kofam)      KOFAM=true ;;
    -p|--phenotrex)  PHENOTREX=true ;;
    -m|--modelseed)  MODELSEED=true ;;
    -d|--dnngior)    DNNGIOR=true ;;
    -a|--all)        PHENOTREX=true; MODELSEED=true; DNNGIOR=true ;;
    --no-flashweave) FLASHWEAVE=false ;;
    -h|--help)       usage; exit 0 ;;
    *)               echo "❌ Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

# --- conda / mamba ---------------------------------------------------------------------------
if command -v mamba >/dev/null 2>&1; then
  CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
  CONDA=conda
else
  echo "❌ conda is not available. Install Miniforge first: https://github.com/conda-forge/miniforge" >&2
  exit 1
fi

env_exists() {
  conda env list | awk '{print $1}' | grep -qx "$1"
}

# Create the environment from envs/<name>.yml, or update it if it already exists
create_env() {
  local name=$1
  local file="$ENVS_DIR/$name.yml"
  if env_exists "$name"; then
    echo "⏳ Updating conda environment '$name' from $file"
    "$CONDA" env update -n "$name" -f "$file"
  else
    echo "⏳ Creating conda environment '$name' from $file"
    "$CONDA" env create -n "$name" -f "$file"
  fi
  echo "✅ Environment '$name' is ready."
}

# --- main environment ------------------------------------------------------------------------
create_env microbetag

ENV_PREFIX=$(conda run -n microbetag python -c "import sys; print(sys.prefix)")

if $FLASHWEAVE; then
  echo "⏳ Installing FlashWeave in Julia (this takes a few minutes the first time)"

  # Keep Julia packages inside the conda environment instead of ~/.julia
  JULIA_DEPOT="$ENV_PREFIX/share/julia-depot"
  conda env config vars set -n microbetag JULIA_DEPOT_PATH="$JULIA_DEPOT" >/dev/null

  # Newer kernels refuse to load libopenlibm with an executable stack
  "$ENV_PREFIX/bin/patchelf" --clear-execstack "$ENV_PREFIX/lib/julia/libopenlibm.so"

  # PYTHON makes PyCall use the environment's Python, which pyjulia needs
  JULIA_DEPOT_PATH="$JULIA_DEPOT" PYTHON="$ENV_PREFIX/bin/python" \
    "$ENV_PREFIX/bin/julia" -e 'using Pkg; Pkg.add(["PyCall", "FlashWeave"]); Pkg.build("PyCall")'

  echo "✅ FlashWeave installed."
fi

# --- optional environments -------------------------------------------------------------------
# microbetag calls these through `conda run -n <name>`, so their names must not change.
$PHENOTREX && create_env mtg-phenotrex
$MODELSEED && create_env mtg-modelseed
$DNNGIOR   && create_env mtg-dnngior

if $MODELSEED && ! command -v rast-create-genome >/dev/null 2>&1; then
  cat <<EOF
🔴 ModelSEEDpy reconstructions from genomes also need RASTtk (rast-create-genome etc.),
   which conda can't install. Install the BV-BRC command-line tools, which needs admin rights:
   https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html
EOF
fi

if $DNNGIOR; then
  echo "🔴 DNNGIOR needs a Gurobi license: https://www.gurobi.com/academia/academic-program-and-licenses/"
fi

# --- kofam database --------------------------------------------------------------------------
if $KOFAM; then
  echo "⏳ Downloading the kofam database"
  (
    cd "$SCRIPT_DIR/ext_data/kofam_database"
    wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    gzip -df ko_list.gz
    tar zxf profiles.tar.gz
  )
  echo "✅ kofam database is in ext_data/kofam_database"
else
  echo "⚪ Skipped the kofam database. It is needed for pathway complementarities with your own genomes:"
  echo "   pass --kofam, or point the configuration YAML to an existing copy."
fi

# --- MetaNetX compounds (normally tracked in git) --------------------------------------------
METANETX="$SCRIPT_DIR/microbetag/mtg_maps_models/MetaNetX/chem_xref.tar.gz"
if [[ ! -f "$METANETX" ]]; then
  echo "⏳ Downloading MetaNetX compounds"
  mkdir -p "$(dirname "$METANETX")"
  wget -q --show-progress -O "$METANETX" "https://zenodo.org/records/15102937/files/chem_xref.tar.gz"
fi

echo ""
echo "🎉 Done. Activate the environment with:  conda activate microbetag"
