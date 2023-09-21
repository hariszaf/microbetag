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

To check whether a species is present on microbetag, one may find its corresponding NCBI Taxonomy Id and search for 
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
using the `/phen-traits/` route and the GTDB genome id. 

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

{: .d-inline-block }

{: .warning }
> Currently, microbetag has annotations only for the GTDB representative genomes. Thus, genomes returned by the `ncbiTaxId-to-genomeId` route that come from other resources (e.g., MGnify, KEGG) do not have phenotypic tratis.


<!-- > 2. In case a GTDB genome returns an "Internal Server Error", please try again replacing the "GCA" with "GCF".  -->






{: .warning }
> [Google Analytics 4 will replace Universal Analytics](https://support.google.com/analytics/answer/11583528). On **July 1, 2023**, standard Universal Analytics properties will stop processing new hits. The earlier you migrate, the more historical data and insights you will have in Google Analytics 4.






## Common errors

There are three possible types of client errors on API calls:

- 400 Bad requests.

- 404 Not found.



