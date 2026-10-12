# microbetag: annotating microbial co-occurrence networks
#
# One image for both the stand-alone tool and the microbetag-nf Nextflow pipeline.
# It holds the same conda environments that install.sh creates (envs/*.yml):
#   - microbetag:    the package, Prodigal, HMMER, DIAMOND, FragGeneScan, CarveMe, manta, Julia + FlashWeave
#   - mtg-phenotrex: phenotrex, called by microbetag with `conda run`
# ModelSEEDpy and DNNGIOR are not included; see containers/ for their images.
#
# Build, from the root of the repository:
#   docker build -t hariszaf/microbetag:<version> .
#
# Run:
#   docker run --rm -v $PWD:/data \
#       -v /path/to/kofam_database:/kofam_database \
#       -v /path/to/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
#       hariszaf/microbetag:<version> microbetag --config /data/config.yml

ARG MINIFORGE_VERSION=25.11.0-1
FROM condaforge/miniforge3:${MINIFORGE_VERSION}

LABEL org.opencontainers.image.title="microbetag" \
      org.opencontainers.image.source="https://github.com/hariszaf/microbetag" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      maintainer="Haris Zafeiropoulos <haris.zafeiropoulos@kuleuven.be>"

ENV SRC=/opt/microbetag \
    MTG_ENV=/opt/conda/envs/microbetag

# phenotrex environment first: it changes least often, so it stays cached
# deepnog weights in a fixed place, readable also when the image is run read-only (Singularity)
ENV DEEPNOG_DATA=/opt/deepnog
COPY envs/mtg-phenotrex.yml ${SRC}/envs/
RUN conda env create -f ${SRC}/envs/mtg-phenotrex.yml \
    # deepnog downloads its weights on first use; fetch them at build time instead
 && conda run -n mtg-phenotrex python -c "from deepnog.utils import get_weights_path; \
get_weights_path(database='eggNOG5', level='2', architecture='deepencoding')" \
 && conda clean -afy

# microbetag environment; the package is installed in editable mode from ${SRC}
COPY pyproject.toml README.md LICENSE ${SRC}/
COPY envs/microbetag.yml ${SRC}/envs/
COPY microbetag/ ${SRC}/microbetag/
RUN conda env create -f ${SRC}/envs/microbetag.yml \
    # GNU parallel and jq are used by the Nextflow pipeline's shell steps
 && conda install -y -n microbetag -c conda-forge parallel jq \
 && conda clean -afy

# FlashWeave inside the environment's Julia depot (as install.sh does)
ENV JULIA_DEPOT_PATH=${MTG_ENV}/share/julia-depot \
    PYTHON=${MTG_ENV}/bin/python
RUN ${MTG_ENV}/bin/patchelf --clear-execstack ${MTG_ENV}/lib/julia/libopenlibm.so \
 && ${MTG_ENV}/bin/julia -e 'using Pkg; Pkg.add(["PyCall", "FlashWeave"]); Pkg.build("PyCall"); using FlashWeave'

# Activate the microbetag environment for every command
ENV PATH=${MTG_ENV}/bin:${PATH} \
    CONDA_DEFAULT_ENV=microbetag \
    # Mount points for the kofam database and the Gurobi license (used by CarveMe)
    MICROBETAG_KOFAM_DB=/kofam_database \
    GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

WORKDIR /data
CMD ["bash"]
