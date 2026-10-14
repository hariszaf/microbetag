---
title: Install microbetag
layout: default
nav_order: 3
description: "Installation approaches"
---


# Install `microbetag`

There are two things to install, depending on [how you use `microbetag`](./modes.md):

- **`MGG`**, the Cytoscape app: to run `microbetag` through a server, and to browse annotated networks;
- **`microbetag`**: one repository with the **package**, the **pipeline** and the **server**, to annotate networks
  with your own genomes, or to run your own server.

| You want to .. | Install |
|---|---|
| annotate your network with the GTDB genomes, on our server | [`MGG`](#install-mgg-cytoscape) |
| prepare a large amplicon dataset for our server | [`MGG`](#install-mgg-cytoscape) and [the pipeline](#the-pipeline) |
| annotate your network with your own genomes | [the pipeline](#the-pipeline) (it runs the package), and `MGG` to browse the result |
| annotate it with the outputs of the tools you already have (GEMs, KEGG annotations, ..) | [the package](#the-package) |
| serve the GTDB genomes or your own genomes to `MGG` | [a server](#a-server) |


## Where to get what

  - `MGG`, the Cytoscape app:
    - the app <a href="https://apps.cytoscape.org/apps/mgg" target="_blank">on the Cytoscape App Store</a>
    - <a href="https://github.com/ermismd/MGG" target="_blank">source code</a>

  - `microbetag`: <a href="https://github.com/hariszaf/microbetag" target="_blank">source code</a>, with
    - the **package**, in `microbetag/`; its image <a href="https://hub.docker.com/r/hariszaf/microbetag" target="_blank">`hariszaf/microbetag`</a> on Docker Hub
    - the **pipeline** (*microbetag-nf*), in `pipeline/`; the images of its tools are on Docker Hub too, and
      `pipeline/get_sif.sh` downloads them for Singularity/Apptainer
    - the **server**, in `server/`, with the schema of its database

  - the data of our server's database, the annotations of the GTDB representative genomes:
    <a href="https://zenodo.org/records/10562677" target="_blank">on Zenodo</a>

  ```{note}
  The pipeline used to be a repository of its own (`msysbio/microbetag-nf`), the amplicon preparation a
  stand-alone image (`hariszaf/microbetag_prep`, up to v1.0.1), and the database and server had their own
  repositories too. They are all in the `microbetag` repository now.
  ```

  Get the repository with:

  ```bash
      git clone https://github.com/hariszaf/microbetag.git
  ```



## Install `MGG` Cytoscape 


  To start using *microbetag* and/or to visualize *microbetag*-annotated networks, you need first, to make sure you have **Cytoscape** on your system; if not, go ahead and [download Cytoscape](https://cytoscape.org/download.html). 

  Then, you need to install the *microbetag* app (`MGG`) from [Cytoscape App store](https://apps.cytoscape.org/apps/mgg).
  Make sure **you first launch Cytoscape** and then visit Cytoscape Appstore.
  If you have already visited the MGG page on Cytoscape Appstore, **launch Cytoscape and refresh the Cytoscape Appstore page**.
  You should now see the **Install** button.

  ![mgg install](_static/img/install_button_mgg.png)

  By clicking it, it will be automatically integrated on your Cytoscape. 
  If you visit Cytoscape Appstore and you have not lunched Cytoscape, you will see a *Download* button instead of the *Install*.
  As already mentioned, we suggest you lunch Cytoscape and refresh the page. 
  Otherwise, you can click the **Download** button and move manually the `.jar` file to the apps folder of your Cytoscape.

  You can also get `MGG` from within Cytoscape by clicking on the `Apps` tab of the main bar and then `App  Store > Show App Store` and typing `microbetag` on the box that pops up.

  Once the app is installed, you may click on the `Apps` tab, and you will find *MGG* there.

  ![mgg_overall](_static/img/app/mainMenu.png)



## Install the pipeline and the package

  To annotate a network using your own genomes, there are two parts:

  - the [**pipeline**](#the-pipeline) runs the external tools on your genomes
    (Prodigal, kofam, phenotrex, CarveMe/gapseq/ModelSEEDpy, FAPROTAX, FlashWeave, manta), and then the package;
  - the [**package**](#the-package) computes the pathway and seed complementarities
    from the outputs of those tools and builds the annotated network.

  To start from your genomes, install the pipeline: it runs the package in its own image, so you do not need to
  install the package too. If you already have the outputs of the tools
  (e.g. GEMs or KEGG annotations of your genomes), the package is enough.


### The pipeline

  You need [Nextflow](https://www.nextflow.io/); if not already available, follow the instructions
  <a href="https://www.nextflow.io/docs/latest/getstarted.html#installation" target="_blank">here</a>.

  You also need a containerization technology to run the tools:
  <a href="https://docs.docker.com/get-docker/" target="_blank">Docker</a>, or
  <a href="https://docs.sylabs.io/guides/3.0/user-guide/installation.html" target="_blank">Singularity</a> /
  <a href="https://apptainer.org/docs/admin/main/installation.html" target="_blank">Apptainer</a>.
  On an HPC system, Singularity or Apptainer is most likely already available, while Docker would not be an option.

  The pipeline is the `pipeline/` folder of the repository:

  ```bash
      git clone https://github.com/hariszaf/microbetag.git
      cd microbetag/pipeline
  ```

  Follow its instructions [here](./tutorials_local/nf.md).
  Depending on the steps you run, the pipeline also needs:

  - the [KOfam database](https://www.genome.jp/ftp/db/kofam/), for the KEGG ORTHOLOGY annotation;
  - a [Gurobi license](#install-gurobi-license), for GEM reconstruction with CarveMe or gap-filling with DNNGIOR
    (or [CPLEX](#install-cplex), for gapseq);
  - the GTDB 16S training set, for the taxonomy assignment of the [prep workflow](./tutorials_otf/prep.md).


### The package

  `microbetag` needs Python 3.10 or later. Install it from its source code:

  ```bash
      git clone https://github.com/hariszaf/microbetag.git
      cd microbetag/microbetag
      pip install .
  ```

  Add the `otf` extra (`pip install ".[otf]"`) for the on-the-fly version, which queries the database of a server,
  and the `test` extra to run the tests (`pytest`, from the `microbetag/` folder of the repository).

  Check that it is available:

  ```bash
      microbetag --version
  ```


### .. as a container

  The `microbetag` package is also available as a Docker image, which also works with Singularity/Apptainer:

  ```bash
      docker pull hariszaf/microbetag:<version>
  ```

  If `<version>` is left out, Docker pulls the latest version.
  See the tutorial on [running `microbetag` on a container](./tutorials_local/containers.md).

  ```{note}
  `microbetag` needs a **configuration YAML file**; get the template of the version you are using with
  `microbetag --write-template config.yml`.
  See [the configuration file](./tutorials_local/config.md).
  ```



## A server

  The server runs in containers, with Docker: the web app, with the package, the pipeline's tools and Nextflow,
  and the MySQL database. It has its own page: [your own server](./tutorials_local/server.md).


## Prepare amplicon data for the on-the-fly version

  For **amplicon** datasets that cannot be analyzed directly on-the-fly, you can run the computationally heavy tasks of

  - the network inference through FlashWeave, and
  - the taxonomy assignment against
    <a href="https://zenodo.org/records/6655692" target="_blank">a GTDB-specific 16S rRNA database</a>,
    which improves the matching of your taxa to the GTDB genomes of our server,

  locally, with the **prep workflow** of the [pipeline](#the-pipeline).
  See [the preparation steps](./tutorials_otf/prep.md).


## Further dependencies 

### Containerization technologies: Docker and Singularity/Apptainer


Most of `microbetag`'s modules are available as containers too. 

So far, we have tested them using:

* 🐳 [Docker](https://docs.docker.com/get-docker/) 

* ⚡[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html): specified for HPC systems


### Install `gurobi` license

  GEM reconstruction with CarveMe, and gap-filling of ModelSEEDpy reconstructions with DNNGIOR, need a solver;
  CarveMe requires one to run (Gurobi or CPLEX). These steps run in the [pipeline](./tutorials_local/nf.md);
  the package itself does not need a solver.

  The reconstruction approaches:

  - **using `carveme`:** works on both DNA and protein sequences and uses the
    <a href="http://bigg.ucsd.edu" target="_blank">BiGG identifiers</a>.

  - **using `modelseedpy`:** requires a RAST annotation of your genomes, which depends on a connection to the RAST server.
    It uses the <a href="https://modelseed.org" target="_blank">ModelSEED resource</a> and its identifiers,
    but it can be rather time-consuming and quite often unsuccessful, due to RAST-related issues.

  In the following sections we provide some links on how to get a Gurobi license (for academics):

#### .. on a local system

  We found 
  <a href="https://www.youtube.com/watch?v=oW6ma8rdZk8" target="_blank">this video</a> 
  (released on 2022) quite helpful on how to get a Gurobi license. 

  Your license is a `gurobi.lic` file. To check that Python can use the license, in an environment with `gurobipy`:

  ```python
      >>> import gurobipy
      >>> m = gurobipy.Model()
      Set parameter Username
      Academic license - for non-commercial use only - expires 2025-04-15
  ```

#### .. on a container

  To use Gurobi in a container, you need a different kind of license, called
  <a href="https://www.gurobi.com/downloads/" target="_blank">**Web License Service (WLS)**</a>.
  You may find the following
  <a href="https://support.gurobi.com/hc/en-us/community/posts/4406485885841-Installing-Gurobi-on-a-Docker-container-Ubuntu" target="_blank">link</a>
  useful on how to get one.

  The WLS license is again a `gurobi.lic` file; you need only the license, not a Gurobi installation.
  The [pipeline](./tutorials_local/nf.md) describes how to give it to the
  CarveMe and DNNGIOR containers.


### Install CPLEX

[IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio)

Jump to [IBM login](https://login.ibm.com/) page. 
If you do not have an account already, create one by clicking on the _Create an IBMid_


![alt text](image.png)



https://academic.ibm.com/a2mt/downloads/data_science#/


https://ronennir.medium.com/installing-cplex-optimization-studio-on-ubuntu-20-04-53e234ca4ec2


https://www.youtube.com/watch?v=o8plELhkazU


