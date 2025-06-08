---
title: microbetag on containers
layout: default
parent: Additional tutorials
nav_order: 3
description: "running microbetag on a container"
---



<!-- Thus, we provide a version of `microbetag` as a stand-alone, 
containerized tool 
so that users can annotate a co-occurrence network using their own bins/MAGs.
To do that, you need first to make sure you have either 
<a href="https://docs.docker.com/get-docker/" target="_blank">Docker</a> 
or <a href="https://docs.sylabs.io/guides/3.0/user-guide/installation.html" target="_blank">Singularity</a>
/
<a href="https://apptainer.org/docs/user/latest/quick_start.html" target="_blank">Apptainer</a> 
in the computing system to be used for running `microbetag`.
The last is common in HPC systems and if you are about to use such a system, you should ask your admin for more information. -->


# on a container

  - Docker / Singularity (containerization technology)
  - the `microbetag` image based on the containerization technology you are using 
  (see below for how to get `microbetag` as a [Docker](#using-docker) or a [Singularity](#using-singularityapptainer) image)


## .. using Docker 


  Once you have installed Docker locally, you may run 
  ```bash
  docker pull hariszaf/microbetag:v1.0.2
  ```
  
  to get microbetag locally.
  
  ```{important}
  **Version is essential!**
  
  Please, make sure you are aware of the version you are using. 
  Latest versions may fix reported bugs or have new features. 
  It is important to always be aware of the version you are using and report it when you are about to submit any issues. 
  ```
  
  Then, you need to get a copy of the `kofam` database to allow the annotation of your sequences with KEGG ORTHOLOGY terms. 
  You may get this by running the following chunk of code: 

  ```bash 
  mkdir kofam_database &&\
  cd kofam_database &&\
  wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz &&\
  wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz &&\
  gzip -d ko_list.gz &&\
  tar zxvf profiles.tar.gz 
  ```


  Now, you need to [download][1] the `config.yml` file that accompanies `microbetag`, 
  to set the values to the required and optional arguments of your choice. 

  In this file, each argument has a `required` field that denotes whether 
  it is mandatory to be set or not. 

  One may provide just an abundance table and the corresponding bins/MAGs sequence files.

  ```{hint}
  **FILENAMES**

  The filenames of your bins/MAGs need to have the same name, like those in your abundance table.
  For example, if in the abundance table you have bin101, then the corresponding filename of the bin should be bin101.fa or bin101.fasta etc.
  This will soon be changed so a mapping file can be used instead. 
  Until then though `microbetag` will fail if that is not the case. 
  ```


  In case you do not already have GENREs for your bins/MAGs, `microbetag` supports two ways for the reconstruction of metabolic networks: 
  
  - using `modelseedpy` that required RAST annotation of your bins and are based on the [ModelSEED resource](https://modelseed.org) and identifiers, 
  
  - using `carveme` that can be performed in both DNA and protein sequences, make use of the [BiGG identifiers](http://bigg.ucsd.edu) and required a Gurobi license (see section [GEM reconstruction step](./local.md#gem-reconstruction-step))

  This can be a rather time-consuming step, especially using `modelseedpy`.

  As you may already have gene predictions for your bins/MAGs, or even protein annotations, you may also provide them to `microbetag`, so those steps can be skipped. 
  If you have already built metabolic networks, then in case they are based on either ModelSEED or BiGG identifiers, you may provide them so seed scores and seed complementarities can be 
  computed directly on them.



  ```{hint}
  **Input folder**

  To conclude, your input folder to be mounted will look like this:

  u0156635@gbw-l-l0074:microbetag$ ls 
  config.yml
  my_bins/
  my_abundance_table.tsv


  where in the `my_bins` folder you have:

  u0156635@gbw-l-l0074:microbetag$ ls bins/
  bin_101.fa
  bin_151.fa
  bin_19.fa
  ```

  Once your input folder is ready, you can mount it on your Docker container and run `microbetag`: 

  ```bash
  docker run --rm -it  \
  --volume=./tests/dev_io_microbetag/:/data \
  --volume=./microbetagDB/ref-dbs/kofam_database/:/microbetag/microbetagDB/ref-dbs/kofam_database/ \
  --volume=$PWD/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
  --entrypoint /bin/bash  \
  hariszaf/microbetag:v1.0.2
  ```

  The `--volume` flag allows you to mount a local directory to a specific path in the container.
  It is **essential** that the **right** parts of the volumes are kept as above! 
  For example, when using `carveme`, a gurobi license is required; `microbetag` expects the license unde the `/opt/gurobi` path, so you need to make sure all the right parts of the volumes are as above and that the left parts point to your local paths. 

  ```{important}
  **Remember!** It is strongly suggested all the files and folders you mount to be part of your root path; meaning the directory from which you initiate your Docker container. 

  For example, if you observe the last chunk of code, you will notice that both `kofam_database` and `gurobi.lic` and the input-output folder called `dev_io_microbetag` they are all within my root folder
  `~/github_repos/KU/microbetag` from where I run the `docker run` command.

  ```

  ```{note}
  A Web License Service (WLS) [Gurobi license](https://www.gurobi.com/downloads/) in case you are about to use `carveme`.
  You may find the following [link](https://support.gurobi.com/hc/en-us/community/posts/4406485885841-Installing-Gurobi-on-a-Docker-container-Ubuntu) useful on how to do that.
  ```

  Once you have fired a container, you can now run `microbetag` using the following command:

  ```bash
  root@20510f8400f1:/microbetag# python3 microbetag.py /data/config.yml 
  ```



## .. using Singularity/Apptainer

  These technologies are widely used in High Performance Computing (HPC) systems. 
  In case you are about to use `microbetag` in such a system, you first need to build a Singularity image (`.simg`) based on the Docker one:

  ```bash
  sudo singularity build microbetag_v102.simg docker://hariszaf/microbetag:v1.0.2
  ```

  You will need to have sudo rights to run this command. 
  If you do not have `sudo` rights, you can either ask your admin to do so or run the build command in a similar environment, e.g. your own Linux based laptop and move it to the HPC with a single `scp` command.
  Also, you can ask your admin or check your HPC documentation site how they deal with Docker images and follow their lead. 

  Once a `.simg` image is available, you may run `microbetag` again by mounting the necessary paths:

  ```bash
  singularity exec 
    -B tests/dev_io_microbetag/:/data  
    -B microbetagDB/ref-dbs/kofam_database/:/microbetag/microbetagDB/ref-dbs/kofam_database/
    -B $PWD/gurobi.lic:/opt/gurobi/gurobi.lic:ro  
    microbetag_v101.simg 
    python3 /microbetag/microbetag.py /data/config.yml
  ```

[1]:../_static/download/local/config.yml
[2]:../_static/download/local/thirty_Samples.tsv
[3]:../_static//download/local/edgelist.csv

