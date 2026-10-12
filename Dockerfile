# microbetag: annotating microbial co-occurrence networks
#
# The microbetag package and its commands; used stand-alone and by the microbetag-nf pipeline.
# External tools (Prodigal, kofam, phenotrex, FlashWeave, CarveMe, FAPROTAX, manta) have their
# own images in microbetag-nf.
#
# Build, from the root of the repository:
#   docker build -t hariszaf/microbetag:<version> .
#
# Run, with the configuration YAML and the files it points to in the current folder:
#   docker run --rm -v $PWD:/data hariszaf/microbetag:<version> microbetag --config /data/config.yml

FROM python:3.10-slim

LABEL org.opencontainers.image.title="microbetag" \
      org.opencontainers.image.source="https://github.com/hariszaf/microbetag" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      maintainer="Haris Zafeiropoulos <haris.zafeiropoulos@kuleuven.be>"

# libexpat is needed by python-libsbml
RUN apt-get update \
 && apt-get install -y --no-install-recommends libexpat1 \
 && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE /opt/microbetag/
COPY microbetag/ /opt/microbetag/microbetag/
RUN pip install --no-cache-dir "/opt/microbetag[otf]" \
 && rm -rf /opt/microbetag

WORKDIR /data
CMD ["bash"]
