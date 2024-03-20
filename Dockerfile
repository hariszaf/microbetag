# microbetag: annotating microbial co-occurrence networks
# 
# Aim:   this Docker image will encapsulate all the related  
#        tools, databases and software modules for the microbetag
#        network annotator
# 
# Usage: docker build -t hariszaf/microbetag:<tag> .

FROM microbetag_base:latest

LABEL maintainer = "Haris Zafeiropoulos" 
LABEL contact    = "haris.zafeiropoulos@kuleuven.be"
LABEL build_date = "2022-12-01"
LABEL version    = "v.1.0"


# Copy microbetag utils 
WORKDIR /microbetag
ADD microbetagDB/mappings/kegg_mappings/*  ./microbetagDB/mappings/kegg_mappings/
ADD microbetagDB/mappings/MetaNetX/chem_xref.tsv ./microbetagDB/mappings/MetaNetX/chem_xref.tsv
ADD microbetagDB/scripts/flashweave.jl ./microbetagDB/scripts/flashweave.jl
ADD utils.py ./
ADD microbetag.py  ./
ADD config.py ./
ADD build_cx_annotated_graph.py ./

ENTRYPOINT [ "python3", "microbetag.py", "/data/config.yml" ]



