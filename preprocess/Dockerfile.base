# microbetagPre: preprocess for running the microbetag annotator of microbial co-occurrence networks
# 
# Aim:   this Docker image will encapsulate all the related  
#        tools, databases and software modules for the microbetag
#        network 
# 
# Usage: docker run --rm -it -v ./test/:/media microbetagpre bash

FROM ubuntu:20.04 

LABEL maintainer = "Haris Zafeiropoulos" 
LABEL contact    = "haris.zafeiropoulos@kuleuven.be"
LABEL build_date = "2023-11-01"
LABEL version    = "v.0.0.1"

# This mode allows zero interaction while installing or upgrading the system via apt; 
# it accepts the default answer for all questions.
ENV DEBIAN_FRONTEND noninteractive
WORKDIR /home

# Get general software; bzip is required to install R 
RUN apt-get update &&\
    apt-get install -y software-properties-common &&\
    apt-get update --fix-missing && \
    apt-get install -y wget vim git unzip mlocate libbz2-dev && \
    apt-get clean &&\
    rm -rf /var/lib/apt/lists/*

# Set Python
RUN add-apt-repository ppa:deadsnakes/ppa &&\
    apt-get install -y python3 python3-pip

# Set R
# First I need to get some R dependencies 
RUN apt-get install -y gfortran build-essential fort77 xorg-dev libblas-dev gcc-multilib && \
    apt-get install -y gobjc++ aptitude libreadline-dev

RUN export CC=/usr/bin/gcc && \
    export CXX=/usr/bin/g++ && \
    export FC=/usr/bin/gfortran && \
    export PERL=/usr/bin/perl

RUN apt-get install -y libpcre3-dev \
                       libpcre2-dev \
                       libpcre-ocaml-dev \
                       libghc-regex-pcre-dev

# Install some extra staff and leave out later what is not needed
RUN apt-get install -y liblzma-dev \
                       libcurl4-openssl-dev \
                       libglib2.0-0 \
                       libxext6 \
                       libsm6 \
                       libxrender1 \
                       mercurial \
                       subversion \
                       autoconf \
                       autogen \
                       libtool \
                       zlib1g-dev

# Install R 
WORKDIR /usr/local/lib/
RUN wget https://ftp.cc.uoc.gr/mirrors/CRAN/src/base/R-3/R-3.6.0.tar.gz
RUN tar -xf R-3.6.0.tar.gz
WORKDIR /usr/local/lib/R-3.6.0
RUN ./configure &&\
    make &&\
    make install

# Install DECIPHER
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) \
    install.packages('BiocManager', dependencies=TRUE, repos='http://cran.us.r-project.org'); \
    BiocManager::install('DECIPHER')"


# Install FlashWeave
# As it is written in Julia we need to get that too
WORKDIR /opt
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.1-linux-x86_64.tar.gz &&\
    tar -zxvf julia-1.7.1-linux-x86_64.tar.gz &&\
    echo "export PATH=/opt/julia-1.7.1/bin:$PATH" >> /root/.bashrc 

ENV PATH="/opt/julia-1.7.1/bin:${PATH}"

# Get FlashWeave
RUN /opt/julia-1.7.1/bin/julia -e 'using Pkg;Pkg.add("FlashWeave")'


# Set paths to mount
# WORKDIR /media/
RUN chmod 777 /media/ &&\
    chmod g+s /media/

RUN pip install PyYAML pandas 

WORKDIR /pre_microbetag
COPY prep.py gtdb_16s.RData classify.R flashweave.jl  ./

ENTRYPOINT [ "/usr/bin/python3", "/pre_microbetag/prep.py" ]
