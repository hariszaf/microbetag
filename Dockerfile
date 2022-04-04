# microbetag: annotating microbial co-occurrence networks
# 
# Aim:   this Docker image will encapsulate all the related  
#        tools, databases and software modules for the microbetag
#        network annotator
# 
# Usage: docker build -t hariszaf/microbetag:<tag> .

FROM ubuntu:20.04 

LABEL maintainer = "Haris Zafeiropoulos" 
LABEL contact    = "haris-zaf@hcmr.gr"
LABEL build_date = "2022-12-01"
LABEL version    = "v.0.0.1"

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


# Set R

## First I need to get some R dependencies 
RUN apt-get install -y gfortran && \
    apt-get install -y build-essential && \
    apt-get install -y fort77 && \
    apt-get install -y xorg-dev && \
    apt-get install -y libblas-dev &&\
    apt-get install -y gcc-multilib && \
    apt-get install -y gobjc++ && \
    apt-get install -y aptitude && \
    aptitude install -y libreadline-dev

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

## Install R 
WORKDIR /usr/local/lib/
RUN wget https://ftp.cc.uoc.gr/mirrors/CRAN/src/base/R-3/R-3.6.0.tar.gz
RUN tar -xf R-3.6.0.tar.gz
WORKDIR /usr/local/lib/R-3.6.0
RUN ./configure &&\
    make &&\
    make install

# Install BugBase dependencies
RUN Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("RColorBrewer2", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("beeswarm", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("reshape2", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("plyr", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("RJSONIO", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("digest", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("optparse", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("Matrix", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("labeling", repos="https://cran.rstudio.com")'


# Install BugBase
WORKDIR /home/software
RUN git clone https://github.com/knights-lab/BugBase.git

RUN echo "export BUGBASE_PATH=/home/software/BugBase" >> /root/.bashrc && \
    echo "export PATH=$BUGBASE_PATH/bin:$PATH" >> /root/.bashrc

ENV PATH="${BUGBASE_PATH}/bin:${PATH}"

RUN Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")'

WORKDIR /usr/local/lib64/R/library
RUN ln -s $PWD/dplyr /home/software/BugBase/R_lib &&\
    ln -s $PWD/RColorBrewer /home/software/BugBase/R_lib &&\
    ln -s $PWD/beeswarm /home/software/BugBase/R_lib &&\
    ln -s $PWD/reshape2 /home/software/BugBase/R_lib &&\
    ln -s $PWD/plyr /home/software/BugBase/R_lib &&\
    ln -s $PWD/gridExtra /home/software/BugBase/R_lib &&\
    ln -s $PWD/RJSONIO /home/software/BugBase/R_lib &&\
    ln -s $PWD/digest /home/software/BugBase/R_lib &&\
    ln -s $PWD/optparse /home/software/BugBase/R_lib &&\
    ln -s $PWD/Matrix /home/software/BugBase/R_lib &&\
    ln -s $PWD/labeling /home/software/BugBase/R_lib &&\
    ln -s $PWD/ggplot2 /home/software/BugBase/R_lib

# Install FAPRTOTAX
WORKDIR /home/software/
RUN wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.4/FAPROTAX_1.2.4.zip &&\
    unzip FAPROTAX_1.2.4.zip &&\
    rm FAPROTAX_1.2.4.zip

# Install Python dependencies for FAPROTAX script
RUN pip install numpy &&\
    pip install pytest-shutil &&\
    pip install biom-format

# Set a text editor as
RUN apt-get install -y vim

# Set paths to mount
WORKDIR /mnt/
RUN chmod 777 /mnt/ &&\
    chmod g+s /mnt/

# Install pandas, dash, plotly 
RUN pip install pandas &&\
    pip install dash &&\
    pip install plotly

# Install EnDED
WORKDIR /home/software
# The boost library is dependency for that
RUN apt-get install -y libboost-dev

# Get and install EnDED
RUN git clone https://github.com/InaMariaDeutschmann/EnDED.git &&\
    cd EnDED &&\
    make

# Install FlashWeave
# As it is written in Julia we need to get that too
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.1-linux-x86_64.tar.gz &&\
    tar -zxvf julia-1.7.1-linux-x86_64.tar.gz &&\
    echo "export PATH=/home/software/julia-1.7.1/bin:$PATH" >> /root/.bashrc 

ENV PATH="/home/software/julia-1.7.1/bin:${PATH}"

# Get FlashWeave
RUN /home/software/julia-1.7.1/bin/julia -e 'using Pkg;Pkg.add("FlashWeave")'


# Install cwl-runner
RUN git clone https://github.com/common-workflow-language/cwltool.git &&\
    cd cwltool &&\
    pip install .[deps] 

RUN pip install cwlref-runner
#------------     SET THE DASH - CYTO - DOCKER SERVER  --------- #

# Set port 
EXPOSE 8050

# Set environmet
ENV NAME world

RUN pip install dash-cytoscape &&\
    pip install dash-vtk &&\
    pip install ipywidgets


# Get the Silva - NCBI Taxonomy Id dump files
RUN wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/taxonomy/taxmap_ncbi_ssu_parc_138.txt.gz


# Install OpenJDK-11
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get install -y ant && \
    apt-get clean;
    
# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME


# Install Graphtools

WORKDIR /home/software
RUN wget http://msysbiology.com/documents/Graphtools.zip && \
    unzip Graphtools.zip && \
    cd Graphtools && \
    export GRAPHTOOLS_ROOT=/home/software/Graphtools/Graphtools && \
    export CLASSPATH=$CLASSPATH:$GRAPHTOOLS_ROOT/lib/NeAT_javatools.jar &&\
    cd REA &&\
    make &&\
    REA_ROOT=$GRAPHTOOLS_ROOT/REA && \
    cd ../kwalks/src/ && \
    make && \
    KWALKS_ROOT=$GRAPHTOOLS_ROOT/kwalks/bin && \
    echo "export CLASSPATH=$CLASSPATH" >> .bashrc && \
    echo "export REA_ROOT=$REA_ROOT" >> .bashrc


# Install cmake and then..
WORKDIR /usr/lib

RUN apt update &&\
    apt purge --auto-remove cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.21.4/cmake-3.21.4.tar.gz &&\
    tar -xzvf cmake-3.21.4.tar.gz

# Get an OpenSSl
RUN apt-get install -y libssl-dev 

WORKDIR /usr/lib/cmake-3.21.4
RUN /bin/bash bootstrap
RUN make -j$(nproc) 
RUN make install

# the App-SpaM placement tool
WORKDIR /usr/lib

RUN git clone https://github.com/matthiasblanke/App-SpaM &&\
    cd App-SpaM &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    make &&\
    echo "export PATH=/usr/lib/App-SpaM/build/:$PATH" >> /root/.bashrc 

ENV PATH="/usr/lib/App-SpaM/build/:${PATH}"


# -----------------------------------
#  ADD WHATEVER BEFORE THE COPIES 
# -----------------------------------

# Copy microbetag utils 
WORKDIR /home

# Add instead of copy.. why?
COPY utils/ ./utils/
COPY tools/ ./tools/
COPY microbetag.py  ./
COPY ref-dbs/silva ./ref-dbs/silva
COPY ref-dbs/kegg_genomes ./ref-dbs/kegg_genomes

# ENV WORKFLOW otu_table



# COPY test/ ./test/
# CMD ["python3", "app.py"]
# # CMD ["cwl-runner", "--debug", "test.cwl", "test-job.yml"]
# # CMD ["cwl-runner", "--debug", "microbetag.cwl", "microbetag-job.yml"]
# CMD ["sh", "-c", "python3 test.py ${WORKFLOW}"]
# CMD ["python3", "scripts/build_a_graph.py", "/mnt/network_output.edgelist"]
# CMD ["python3", "scripts/pass_networkx_to_dash.py"]






# ----------------------------------------------------------------

# WORKDIR /home/software/FAPROTAX_1.2.4/
# RUN sed -i "208s/return (s.lower() is not 'nan') and is_number(s);/return (s.lower() != 'nan' and is_number(s))/g" collapse_table.py


# # Format 
# RUN pip install pyarrow


# # In case Picrust2 is to be integrated // needs to be discussed.. 

# # First, we need to set conda 

# # Install miniconda
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#      /bin/bash ~/miniconda.sh -b -p /opt/conda 
# RUN echo "export PATH=$CONDA_DIR/bin:$PATH" >> /root/.bashrc 

# # Then, get Picrust2
# WORKDIR /home/software
# RUN wget https://github.com/picrust/picrust2/archive/v2.4.2.tar.gz &&\
#     tar xvzf  v2.4.2.tar.gz &&\
#     cd picrust2-2.4.2/ &&\
#     conda env create -f picrust2-env.yaml

# RUN echo "conda activate picrust2" >> ~/.bashrc
# # SHELL ["/bin/bash", "--login", "-c"]   
# RUN pip install --editable . 
