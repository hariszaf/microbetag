# microbetag: annotating microbial co-occurrence networks
# 
# Aim:   this Docker image will encapsulate all the related  
#        tools, databases and software modules for the microbetag
#        network annotator
# 
# Usage: docker build -t hariszaf/microbetag:<tag> .

FROM ubuntu:20.04 

LABEL maintainer = "Haris Zafeiropoulos" 
LABEL contact    = "haris.zafeiropoulos@kuleuven.be"
LABEL build_date = "2022-12-01"
LABEL version    = "v.0.0.1-dev"

# This mode allows zero interaction while installing or upgrading the system via apt; it accepts the default answer for all questions.
ENV DEBIAN_FRONTEND noninteractive
WORKDIR /home

# Get general software; bzip is required to install R 
RUN apt-get update &&\
    apt-get install -y software-properties-common &&\
    apt-get update --fix-missing && \
    apt-get install -y wget && \
    apt-get install -y git && \
    apt-get install -y unzip && \
    apt-get install -y mlocate && \
    apt-get install -y libbz2-dev && \
    apt-get clean &&\
    rm -rf /var/lib/apt/lists/*

# Set Python
# Install py39 from deadsnakes repository
# Install pip from standard ubuntu packages
RUN add-apt-repository ppa:deadsnakes/ppa &&\
    apt-get install -y python3 &&\
    apt-get install -y python3-pip

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

# PhenDB
RUN pip install phenotrex[fasta]

# RUN wget https://zenodo.org/records/10562677/files/phen_classes.zip &&\
#     unzip phen_classes.zip &&\
#     rm phen_classes.zip

# Set a text editor as
RUN apt-get install -y vim

# Set paths to mount
WORKDIR /data/
RUN chmod 777 /data/ &&\
    chmod g+s /data/

# Install pandas, dash, plotly 
RUN pip install pandas

# -----------------------------------
#  ADD WHATEVER BEFORE THE COPIES 
# -----------------------------------

WORKDIR /microbetag/microbetagDB/ref-dbs/phenDB/classes/
ADD microbetagDB/ref-dbs/phenDB/classes/*  ./

RUN sed -i 's/np\.integer/int/g' /usr/local/lib/python3.8/dist-packages/phenotrex/io/flat.py  &&\
    sed -i 's/np\.floating/float/g' /usr/local/lib/python3.8/dist-packages/phenotrex/io/flat.py  &&\
    sed -i 's/np\.int/int/g' /usr/local/lib/python3.8/dist-packages/deepnog/data/dataset.py &&\
    sed -i 's/np\.int/int/g' /usr/local/lib/python3.8/dist-packages/deepnog/learning/training.py

RUN pip install networkx[default] &&\
    pip install scikit-learn==1.2.2

WORKDIR /usr/local
RUN git clone https://github.com/hyattpd/Prodigal.git 
RUN cd Prodigal &&\ 
    make install &&\ 
    make

RUN pip install scikit-bio

RUN apt-get install infernal infernal-doc

# To install DRAM (for KEGG annotation of the bins), we need the following
# pandas [done], networkx[done], scikit-bio[done], prodigal[done], mmseqs2, hmmer and 
# tRNAscan-SE (for ) === [TODO] CHECK IF ALL NEEDED 
RUN wget http://lowelab.ucsc.edu/software/trnascan-se-2.0.12.tar.gz &&\ 
    gunzip trnascan-se-2.0.12.tar.gz  &&\
    tar xf trnascan-se-2.0.12.tar &&\
    cd tRNAscan-SE-2.0/ &&\
    ./configure &&\
    make &&\
    make install

RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
ENV PATH="/usr/local/mmseqs/bin:${PATH}"

RUN wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz ; tar xf hmmer-3.4.tar.gz ; cd hmmer-3.4 ; ./configure ; make ; make install

# # This needs to be performed locally and then mounted; otherwise we go to an image of 15G
# RUN mkdir kofam_database &&\
#     cd kofam_database &&\
#     wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz &&\
#     wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz &&\
#     gzip -d ko_list.gz &&\
#     tar zxvf profiles.tar.gz 

WORKDIR /microbetag
RUN git clone https://github.com/xuechunxu/DiTing.git


# Copy microbetag utils 
WORKDIR /microbetag
ADD microbetagDB/mappings/kegg_mappings/*  ./microbetagDB/mappings/kegg_mappings/
ADD utils.py ./
ADD microbetag.py  ./
ADD config.py ./


ENTRYPOINT [ "python3", "microbetag.py", "/data/config.yml" ]



# ==============================================================  to be added
# # FlashWeave
# # As it is written in Julia we need to get that too
# WORKDIR /opt
# RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.1-linux-x86_64.tar.gz &&\
#     tar -zxvf julia-1.7.1-linux-x86_64.tar.gz &&\
#     echo "export PATH=/opt/julia-1.7.1/bin:$PATH" >> /root/.bashrc 

# ENV PATH="/opt/julia-1.7.1/bin:${PATH}"

# # Get FlashWeave
# RUN /opt/julia-1.7.1/bin/julia -e 'using Pkg;Pkg.add("FlashWeave")'

# # FAPRTOTAX
# WORKDIR /opt
# RUN wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.6/FAPROTAX_1.2.6.zip &&\
#     # https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.4/FAPROTAX_1.2.4.zip &&\
#     unzip FAPROTAX_1.2.6.zip &&\
#     rm FAPROTAX_1.2.6.zip
# # Install Python dependencies for FAPROTAX script
# RUN pip install "numpy<1.24" &&\
#     pip install pytest-shutil &&\
#     pip install biom-format

# # Install cmake and then..
# WORKDIR /usr/lib

# RUN apt update &&\
#     apt purge --auto-remove cmake &&\
#     wget https://github.com/Kitware/CMake/releases/download/v3.21.4/cmake-3.21.4.tar.gz &&\
#     tar -xzvf cmake-3.21.4.tar.gz

# # Get an OpenSSl
# RUN apt-get install -y libssl-dev 

# WORKDIR /usr/lib/cmake-3.21.4
# RUN /bin/bash bootstrap
# RUN make -j$(nproc) 
# RUN make install


# ----------------------------------------------------------------

# WORKDIR /home/software/FAPROTAX_1.2.4/
# RUN sed -i "208s/return (s.lower() is not 'nan') and is_number(s);/return (s.lower() != 'nan' and is_number(s))/g" collapse_table.py

# # BugBase 
# # Dependencies
# RUN Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com")' &&\
#     Rscript -e 'install.packages("RColorBrewer2", repos="https://cran.rstudio.com")' &&\
#     Rscript -e 'install.packages("beeswarm", repos="https://cran.rstudio.com")' &&\
#     Rscript -e 'install.packages("reshape2", repos="https://cran.rstudio.com")' &&\
#     Rscript -e 'install.packages("plyr", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("RJSONIO", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("digest", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("optparse", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("Matrix", repos="https://cran.rstudio.com")' && \
#     Rscript -e 'install.packages("labeling", repos="https://cran.rstudio.com")' &&\
#     Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")'


# ==========================================================================================================
# # Set R
# # First I need to get some R dependencies 
# RUN apt-get install -y gfortran && \
#     apt-get install -y build-essential && \
#     apt-get install -y fort77 && \
#     apt-get install -y xorg-dev && \
#     apt-get install -y libblas-dev &&\
#     apt-get install -y gcc-multilib && \
#     apt-get install -y gobjc++ && \
#     apt-get install -y aptitude && \
#     aptitude install -y libreadline-dev

# RUN export CC=/usr/bin/gcc && \
#     export CXX=/usr/bin/g++ && \
#     export FC=/usr/bin/gfortran && \
#     export PERL=/usr/bin/perl

# RUN apt-get install -y libpcre3-dev \
#     libpcre2-dev \
#     libpcre-ocaml-dev \
#     libghc-regex-pcre-dev

# # Install R 
# WORKDIR /usr/local/lib/
# RUN wget https://ftp.cc.uoc.gr/mirrors/CRAN/src/base/R-3/R-3.6.0.tar.gz
# RUN tar -xf R-3.6.0.tar.gz
# WORKDIR /usr/local/lib/R-3.6.0
# RUN ./configure &&\
#     make &&\
#     make install


# ==========================================================================================================


# ==========================================================================================================

# # Install OpenJDK-11
# RUN apt-get update && \
#     apt-get install -y openjdk-11-jdk && \
#     apt-get install -y ant && \
#     apt-get clean;

# # Fix certificate issues
# RUN apt-get update && \
#     apt-get install ca-certificates-java && \
#     apt-get clean && \
#     update-ca-certificates -f;

# # Setup JAVA_HOME -- useful for docker commandline
# ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
# RUN export JAVA_HOME

# ==========================================================================================================

