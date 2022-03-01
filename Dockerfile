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
    apt-get install -y wget \ 
                       git \
                       unzip \
                       mlocate \ 
                       libbz2-dev
                     
# Set Python
RUN add-apt-repository ppa:deadsnakes/ppa &&\
# Install py39 from deadsnakes repository
    apt-get install -y python3 &&\
    # Install pip from standard ubuntu packages
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
WORKDIR /home/external_tools
RUN git clone https://github.com/knights-lab/BugBase.git

RUN echo "export BUGBASE_PATH=/home/external_tools/BugBase" >> /root/.bashrc && \
    echo "export PATH=$PATH:$BUGBASE_PATH/bin" >> /root/.bashrc

RUN Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")'

WORKDIR /usr/local/lib64/R/library
RUN ln -s $PWD/dplyr /home/external_tools/BugBase/R_lib &&\ 
    ln -s $PWD/RColorBrewer /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/beeswarm /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/reshape2 /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/plyr /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/gridExtra /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/RJSONIO /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/digest /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/optparse /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/Matrix /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/labeling /home/external_tools/BugBase/R_lib &&\
    ln -s $PWD/ggplot2 /home/external_tools/BugBase/R_lib

# Install FAPRTOTAX
WORKDIR /home/external_tools/
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
WORKDIR /home/external_tools
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
    echo "export PATH=$PATH:/home/external_tools/julia-1.7.1/bin" >> /root/.bashrc 
# Get FlashWeave
RUN /home/external_tools/julia-1.7.1/bin/julia -e 'using Pkg;Pkg.add("FlashWeave")'


# Install cwl-runner
RUN git clone https://github.com/common-workflow-language/cwltool.git &&\
    cd cwltool &&\
    pip install .[deps] 

RUN pip install cwlref-runner
#------------     SET THE DASH - CYTO - DOCKER SERVER  --------- #

# Copy microbetag app 
WORKDIR /app

# Add instead of copy.. why?
COPY app/ ./

# Set port 
EXPOSE 8050

# Set environmet
ENV NAME world

RUN pip install dash-cytoscape &&\
    pip install dash-vtk &&\
    pip install ipywidgets


RUN ln -s /home/external_tools/FAPROTAX_1.2.4/collapse_table.py /app/tools/faprotax

ENV WORKFLOW otu_table
COPY test/ ./test/
CMD ["python3", "app.py"]
# # CMD ["cwl-runner", "--debug", "test.cwl", "test-job.yml"]
# # CMD ["cwl-runner", "--debug", "microbetag.cwl", "microbetag-job.yml"]
# CMD ["sh", "-c", "python3 test.py ${WORKFLOW}"]
# CMD ["python3", "scripts/build_a_graph.py", "/mnt/network_output.edgelist"]
# CMD ["python3", "scripts/pass_networkx_to_dash.py"]
