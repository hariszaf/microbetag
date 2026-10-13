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


  If you encounter any issues, please don’t hesitate to reach out to us via our
  <a href="https://matrix.to/#/#microbetagcommunity:matrix.org" target="_blank">Matrix community</a>.

  Bug reports, feature suggestions, and general feedback—especially from the diverse ways in which `microbetag` may be used are the most welcome too.



## ✅ Check availability

  Once installed, `microbetag` should be available when running:

  ```{code}
    $ microbetag --help
    usage: microbetag [-h] [--config CONFIG] [-v]

    Microbetag CLI

    options:
    -h, --help            show this help message and exit
    --config CONFIG, -c CONFIG
                            Path to the configuration yaml file.
    -v, --version         Show Microbetag version
  ```

  If not, please follow instructions [here](../installation.md#the-package) on how to get it! 



## <img src="../_static/img/Python_Logo.png" style="width: 1em; height: 1em;"> Use `microbetag`'s features as a library

  You can now use `microbetag`'s features independently, for example the routines that export
  pathway and/or seed complementarities, or that build the annotated network.
  The same steps are available as commands (`microbetag-pcompl`, `microbetag-seed-sets`, `microbetag-annotate`, ...),
  which the [microbetag-nf](./nf.md) pipeline uses; each has `--help`.

  To go through them, you may have a look at the [API Reference](../autoapi/index) but also 
  on the <a href="https://github.com/msysbio/microbetag/tree/develop/microbetag/tests" target="_blank">`unittests` of ours</a>.

  Input/output for these tests can be found in the 
  <a href="https://github.com/msysbio/microbetag/tree/develop/microbetag/tests/data" target="_blank">`tests/data`</a>
  folder where each case is explained further. 
  
  In the [.. using custom genomes](./local.md) tutorial, we provide an example of a complete `microbetag` run,
  from a network and the outputs of the external tools to an annotated network.

  `microbetag` is also available as a container; see [running `microbetag` on a container](./containers.md).




  ```{important}
  Remember that `microbetag` will use all files that have been produced from a previous run, 
  given the output directory stays the same! 
  That means that if you wish to re-calculate a part for any reason, 
  e.g. with different GEMs, you need to make sure
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
