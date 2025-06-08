---
title: Run from Python
layout: default
parent: Additional tutorials
nav_order: 4
description: "how to run microbetag from Python"
---


# 🎯 `microbetag` locally

  Running `microbetag` locally provides a great range of features and can be used as an interface for several tasks,
  network annotation related and not. 


  Due to its many dependencies and the variety of possible input and output configurations, 
  installing and running `microbetag` across different operating systems and platforms can be challenging. 
  Performance may vary depending on the specific setup too.

  If you encounter any issues, please don’t hesitate to reach out to us via our
  <a href="https://matrix.to/#/#microbetagcommunity:matrix.org" target="_blank">Matrix community</a>.

  Bug reports, feature suggestions, and general feedback—especially from the diverse ways in which `microbetag` may be used are the most welcome too.



## ✅ Check availability

  Once installed, `microbetag` should be available when running:

  ```{code}
    (microbetag) u0156635@gbw-l-l0074:~$ microbetag --help
    usage: microbetag [-h] [--config CONFIG] [-v]

    Microbetag CLI

    options:
    -h, --help            show this help message and exit
    --config CONFIG, -c CONFIG
                            Path to the configuration yaml file.
    -v, --version         Show Microbetag version
  ```

  If not, please follow instructions [here](../installation.md#install-microbetag-locally) on how to get it! 



## <img src="../_static/img/Python_Logo.png" style="width: 1em; height: 1em;"> Use `microbetag`'s features as a library

  You can now use `microbetag`'s features independently: either general tasks required for the actual annotation steps
  or annotation steps on their own. 
  For example, you may run individual genome annotation steps, build your own GENREs, 
  as well as the routines allowing you to export pathway and/or seed complementarities. 

  To go through them, you may have a look at the [API Reference](../autoapi/index) but also 
  on the <a href="https://github.com/msysbio/microbetag/tree/develop/tests" target="_blank">`unittests` of ours</a>.

  Input/output for these tests can be found in the 
  <a href="https://github.com/msysbio/microbetag/tree/develop/test_data" target="_blank">`test_data`</a>
  folder where each case is explained further. 
  
  In the [.. using custom genomes](./local.md) tutorial, we provide an example on how to run a complete 
  `microbetag`-pipeline, starting from an abundance table and a list of genome sequencing files, 
  and end up with an annotated network. 

  `microbetag` and its required environments are also available in containerized versions. 
  If you wish to install it as such, follow instructions [here](../installation#containerization-technologies-docker-and-singularityapptainer).
  If already installed, you may jump to the corresponding tutorial on [running `microbetag` on a container](./containers.md).




  ```{important}
  Remember that `microbetag` will use all files that have been produced from a previous run, 
  given the output directory stays the same! 
  That means that if you wish to re-calculate a part for any reason, 
  e.g. run FlashWeave again with a different set of parameters, you need to make sure
  that you first **delete** the previous output files of this step! 
  ```







<!-- 
  For example, one can easily run the whole _microbetag_ pipeline 
  using the corresponding [`config.yml`](./config.md) template to the `microbetag` version to be used, running:

  ```python
  from microbetag.config import Config
  from microbetag.microbetag import run_microbetag

  with open("config.yml", "r") as yaml_file:
    config = Config(yaml.safe_load(yaml_file), "config.yml")

  run_microbetag(config)
  ```

  or simply by

  ```bash
  microbetag --config config.yml
  ``` -->



## 🤖 Use `microbetag` API programmatically

  One may use the API routes described through Python. 
  For example, to get the _microbetagDB_ genomes related species with NCBI Taxonomy ID 853:

  ```python
  import requests
  url = "https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId/853"
  r = requests.get(url)
  r.status_code
  200
  r.json()
  {"537007": ["GCA_002222595.2"]}
  ```

  For more about how to use the _microbetagDB_ API you may check [here](../autoapi/api).
