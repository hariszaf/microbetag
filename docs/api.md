---
layout: default
title: API documentation
nav_order: 4
---

# API documentation
{: .no_toc }

<!-- Just the Docs has some specific configuration parameters that can be defined in your Jekyll site's _config.yml file. -->
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



microbetagDB API provides programmatic access to the data. 

The base address to the API is https://msysbio.gbiomed.kuleuven.be/.

Below you may find the syntax to retrieve the various data and/or annotations included.


## Get genome ids for a NCBI Taxonoy Id

To check whether a species is present on microbetag, one may find its corresponding **NCBI Taxonomy Id** and search for 
related genomes present on the microbetag DB. 

For example, assuming we are interested in the *Blautia hansenii* DSM 20583 strain, we find from NCBI Taxonomy that its corresponding id is [537007](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=537007).

Using the `ncbiTaxId-to-genomeId` route we may get the related genomes on microbetag DB:

```bash
curl -X GET https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId/537007
```
that returns a list of genomes used in microbetag annotations:

```bash
{
  "853": [
    "GCA_002222595",
  ]
}
```

Similarly, in Python:

```python
>>> import requests
>>> url = "https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId/853"
>>> r = requests.get(url)
>>> r.status_code
200
>>> r.json()
{"537007": ["GCA_002222595.2"]}
```


## Get phenotypic traits of a species 

One may get the corresponding phenotypic traits of a genome present on the microbetag DB 
using the `/phen-traits/` route and the **GTDB genome id**.

```bash
curl -X GET https://msysbio.gbiomed.kuleuven.be/phen-traits/GCF_002222595.2
```

this would return (we show only a part of the outcome)

```bash
{
  "NOB": "NO",
  "NOBScore": "0.8078",
  "T3SS": "NO",
  "T3SSScore": "0.8391",
  "T6SS": "NO",
  "T6SSScore": "0.7836",
  "aSaccharolytic": "NO",
  "aSaccharolyticScore": "0.8672",
   ...
}
```

{: .note}
For a thorough description of the abbreviations used, have a look in the microbetag's [modules tab](modules.md#based-on-phendb).


{: .warning }
> 1. Currently, microbetag has annotations only for the GTDB representative genomes. Thus, genomes returned by the `ncbiTaxId-to-genomeId` route that come from other resources (e.g., MGnify, KEGG) do not have phenotypic tratis.
> 2. In case a GTDB genome returns an "Internal Server Error", please try again replacing the "GCA" with "GCF". 




## Get pathway complementarities for a pair of genomes 

In case you are interested in the complementarities of a specific pair of GTDB representative genomes, where `genome_A` is the beneficiary and `genome_B` the donor, 
one may use the `genome-complements` route, followed by `genome_A`, followed by `genome_B`. 

Here is an example: we are interested in the pathway complementarities of a Desulfurococcaceae archeon genome (GCA_011364525.1) with on of the Gram-negative bacterium *Malikia spinosa* (GCA_002980625.1) as its potential donor:


```bash
curl -X GET https://msysbio.gbiomed.kuleuven.be/genome-complements/GCA_011364525.1/GCA_002980625.1
```

Here is its partial output:

```bash
  [
    [
      "M00001",
      "K01803;K21071;K01624",
      "K25026;K15916;K21071;K01624;K01803;K00150;K00927;K15635;K01689;K00873"
    ]
  ],
  [
    [
      "M00002",
      "K01803",
      "K01803;K00150;K00927;K15635;K01689;K00873"
    ]
  ],
  [
    [
      "M00003",
      "K01803",
      "K01596;K01689;K15635;K00927;K00150;K01803;K01622"
    ]
  ],
  [
    [
      "M00016",
      "K00928;K00215;K01586;K00674;K00821;K01778",
      "K00928;K00133;K01714;K00215;K00674;K00821;K01439;K01778;K01586"
    ]
  ],
  [
    [
      "M00016",
      "K00928;K14267;K00215;K01586;K00674;K01778",
      "K00928;K00133;K01714;K00215;K00674;K14267;K01439;K01778;K01586"
    ]
  ]
```

Let us now describe its meaning. 
Each entry of the above output stands for a specific KEGG module






## Common errors

There are three possible types of client errors on API calls:

- 400 Bad requests.

- 404 Not found.



