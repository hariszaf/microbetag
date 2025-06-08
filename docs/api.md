---
layout: default
title: Access *microbetagDB* through its API
---

# *microbetagDB* API

  *microbetagDB* API provides programmatic access to the data.
  Using the Application Programming Interface (API) you can access the *microbetagDB* directly to get information about 
  genome-based predicted traits of a specific taxon, potential pathway complementarities of a taxa pair etc. 

  The base address to the API is <a href="https://msysbio.gbiomed.kuleuven.be/" target="_blank">https://msysbio.gbiomed.kuleuven.be/</a>.

  Below you may find the syntax to retrieve the various data and/or annotations included.

  Remember that *microbetag* is a NCBI Taxonomy oriented resource. 
  That means that a "species" of interest is a NCBI Taxonomy ID.
  For example, if you are interested in *Bifidobacterium animalis*, you first need to go to the 
  <a href="https://www.ncbi.nlm.nih.gov/taxonomy/" target="_blank">NCBI Taxonomy portal</a> 
  and get its corresponding ID.
  However, once you do so you get a list of subspecies and strains available. 
  You can either use the species ID or one of a specific strain in your queries. 

  <!-- *microbetag* has a special feature called `get_children` as in some cases there is no genomic information for the species level, 
  but there is at lower levels. 
  For example, in case of *Bifidobacterium animalis*, *microbetagDB* has no genomes for its corresponding NCBI Taxonomy ID (28025) 
  but it does have for the *Bifidobacterium animalis* subsp. animalis ATCC 25527 (703613).  -->



  <!-- This example highlights the issues derived from the current taxonomy schemes as 
  Even more interestingly, according to GTDB 
  *Bifidobacterium canis* -->

## Get genome IDs for a NCBI Taxonomy ID

  To check whether a species is present on *microbetag*, one may find its corresponding **NCBI Taxonomy ID** 
  and search for related genomes present on the *microbetagDB*. 

  For example, assuming we are interested in the *Blautia hansenii* DSM 20583 strain, 
  we find from NCBI Taxonomy that its corresponding ID is 
  <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=537007" target="_blank">537007</a>.

  Using the `ncbiTaxId-to-genomeId` route we may get the related genomes on *microbetagDB*:

  ```bash
  curl https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId?ncbiTaxId=537007
  ```
  <!--   # curl http://localhost:1337/ncbiTaxId-to-genomeId?ncbiTaxId=537007 -->
  that returns a list of genomes used in *microbetag* annotations:

  ```bash
  {
    "537007": [
      "GCF_002222595.2"
    ]
  }
  ```

  In this case there is a genome available for this NCBI Taxonomy ID (GCF_002222595.2).

  If no genome is available, then you will get a `NotFoundError`.
  For example, for the *Bifidobacterium animalis* subsp. animalis IM386 (NCBI Tax ID: 1402194)
  there is no genome in *microbetagDB*:

  ```bash
  curl https://msysbio.gbiomed.kuleuven.be/ncbiTaxId-to-genomeId?ncbiTaxId=1402194
  {
    "error": "NCBI Taxonomy Id 1402194 was not found in microbetagDB."
  }

  ```


## 🧬 Genome-based predicted phenotypic traits

  Once you have identified the genomes related to your NCBI Taxonomy ID 
  under study using the `ncbiTaxId-to-genomeId` route, 
  you may get the corresponding phenotypic traits of that genome(s) 
  using the `/phen-traits/` route and the corresponding **genome ID**.

  For example, in case of *Blautia hansenii* DSM 20583 (NCBI Taxonomy ID: 537007) 
  we saw there is a genome on *microbetagDB*; to get its phenotypic traits we can simply run:

  ```bash
  curl https://msysbio.gbiomed.kuleuven.be/phen-traits?ncbiAssemblyId=GCF_002222595.2

  ```
  <!--   # curl http://localhost:1337/phen-traits?ncbiAssemblyId=GCF_002222595.2 -->
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

  ```{note}
  For a thorough description of the abbreviations used, have a look in the *microbetag*'s [modules tab](modules/modules.md#based-on-phendb).
  ```

  ```{warning}
  1. Currently, microbetag has annotations only for the GTDB representative genomes. 
  Thus, genomes returned by the `ncbiTaxId-to-genomeId` route that come from other resources 
  (e.g., MGnify, KEGG) do not have phenotypic tratis.

  1. All genomes have a 3-letter prefix that is either GCA or GCG. In case a GTDB genome 
  you are querying returns an "Internal Server Error", try again replacing that prefix; 
  e.g. if you initially had "GCA_002222595.2", try again with "GCF_002222595.2". 
  ```

  In case a genome ID is provided for which there are no phenotypic traits on *microbetagDB*, you will get a message explaining this:

  ```bash
  No Phen traits for the genome id asked.         
  Make sure you are asking for a GTDB v202 representative genome.
  ```


## 🧩 Get pathway complementarities

  In case you are interested in the complementarities of a specific pair of GTDB representative genomes, 
  where $species_A$ is the ***beneficiary***, and $species_B$ the ***donor***,
  one may use the `genome-complements` route, providing their corresponding IDs, and the type of the IDs provided.


  * 🔧 **Query Parameters**

  | Parameter       | Type     | Required | Description                                               |
  |:---------------:|:--------:|:--------:|:---------------------------------------------------------:|
  | `beneficiaryId` | `string` | ✅        | ID of the potential beneficiary (e.g. `GCA_011050685.1`) |
  | `donorId`       | `string` | ✅        | ID of the potential donor (e.g. `1402194`)               |
  | `typeCategory`  | `string` | ✅        | `ncbiGenomeIds \| ncbiTaxonomyIds \| patricGenomeIds`    |


  ```{warning}
  Both IDs need to be of the same type.
  ```


  * 🔍 **Example Request**

  Here is an example where _Desulfurococcaceae archeon_ genome (GCA_011364525.1) is the beneficiary taxon,
  and *Malikia spinosa* (GCA_002980625.1) its potential donor:

  ```bash
  curl "https://msysbio.gbiomed.kuleuven.be/pathway-complements?beneficiaryId=GCA_011050685.1&donorId=GCA_002289575.1&typeCategory=ncbiGenomeIds"
  ```
  <!--   # curl "http://localhost:1337/pathway-complements?beneficiaryId=GCA_011050685.1&donorId=GCA_002289575.1&typeCategory=ncbiGenomeIds" -->

  ```{warning}
  Please make sure you have the quotation marks in the beginning and the end of your URL.
  ```


  Here is its partial output:

  ```bash
  {
    "0": {
      "beneficiary-ncbi-tax-id": 182270,
      "complements": [
        {
          "beneficiary-genome": "GCA_011050685.1",
          "complements": [
            [
              [
                "M00001",
                "K21071",
                "K00845;K01810;K21071;K01624;K01803;K00134;K00927;K01834;K01689;K00873",
                "https://www.kegg.jp/kegg-bin/show_pathway?map00010/K00845%09%23EAD1DC/K01810%09%23EAD1DC/K01624%09%23EAD1DC/K01803%09%23EAD1DC/K00134%09%23EAD1DC/K00927%09%23EAD1DC/K01834%09%23EAD1DC/K01689%09%23EAD1DC/K00873%09%23EAD1DC/K21071%09%2300A898/"
              ]
            ],
          #....  MORE COMPLEMENT CASES 
          ],
          "donor-genome": "GCA_002289575.1"
        }
      ],
      "donor-ncbi-tax-id": 161920
    }
  }


  ```

  Let us now describe its content. 

  * The response returns an entry for each beneficiary entry ID

  * Each entry consists of the `beneficiary-ncbi-tax-id` and the `donor-ncbi-tax-id`, the NCBI Taxonomy IDs of the two taxa

  * For each such pair there is a`complements` list. Each entry of which stands for  
    a combination of `beneficiary-genome` - `donor-genome` IDs and has again a `complements` list with the actual pathway complementarities found 
    between the two genomes

  * Each entry in that list has 4 elements: 
    * the KEGG module the complement is part of
    * the KO terms(s) to be exchanged
    * he complete *alternative* (set of KOs that ensure the module)
    * a URL to a colored KEGG map that includes the module and highlights KOs of the benefciary and those to get from the donor

  * Under the `complements` key, you can find the different potential complements between those two genomes. 
  In the above chunk we only show the first two. 
  Each complement refers to a specific KEGG module, denoted in the `module` key.
  Each complement entry, has also a `complement` key with the exact KO terms that the donor would have to provide, so the donor would have a complete alternative of the module. 
  In the `complete-alternative` you can find all the KOs (both those the beneficiary carries on its own and those it would get from the donor) to have a complete alternative of the module under study. 
  Last, the `coloured-map` key provides you the link to the KEGG map showing the module under study and the KO terms to be used; beneficiary's KOs are coloured with pink while those to get from the donor with green.


  ```{warning}
  When using `curl` make sure you use quotes ("<URL>") around your URL.
  ```


## <img src="_static/img/app/seed.jpg" style="width: 1em; height: 1em; vertical-align: text-bottom;"> Get seed scores and complements

### Get competition and complementarity seed scores 
  
  When we calculate the seed scores, we consider both taxa as $species_A$ and $species_B$   
  (see on the 
  [Modules tab](./modules/modules.md#seed-scores-and-complements-based-on-genome-scale-draft-reconstructions-gems) 
  for more).


  Seed scores can be retrieved in a similar to the complements way.
  For example, let's check the scores between *Streptomyces* sp. AW19M42 (NCBI Taxonomy ID: 1379686) 
  and *Afipia clevelandensis* ATCC 49720 (NCBI Taxonomy ID: 883079). 
  By running:

  ```bash
  curl "https://msysbio.gbiomed.kuleuven.be/seed-scores?beneficiaryId=1379686&donorId=883079&typeCategory=ncbiTaxonomyIds"
  ```
  <!--     # curl "http://localhost:1337/seed-scores?beneficiaryId=1379686&donorId=883079&typeCategory=ncbiTaxonomyIds" -->

  we get

  ```bash
  {
    "beneficiary_maps": {
      "gc_and_patric_ids": {
        "GCF_000470535.1": "1379686.4"
      },
      "ncbi_tax_id": "1379686"
    },
    "donor_maps": {
      "gc_and_patric_ids": {
        "GCF_000336555.1": "883079.3"
      },
      "ncbi_tax_id": "883079"
    },
    "seed-scores": {
      "0": {
        "CompetitionScore": "0.51",
        "CooperationScore": "0.33",
        "PATRIC_A": "1379686.4",
        "PATRIC_B": "883079.3"
      },
      "1": {
        "CompetitionScore": "0.67",
        "CooperationScore": "0.2",
        "PATRIC_A": "883079.3",
        "PATRIC_B": "1379686.4"
      }
    }
  }
  ```

  where, in the first case `["seed-scores"][0]`, *Streptomyces* is considered as $speciesA$ 
  and *Afpia* as $speciesB$, and in case `[1]` the other way around. 


  If I run the same with the reverse order on the Tax IDs,
  then I get the same output only with a different order. 

  ```bash
  curl "https://msysbio.gbiomed.kuleuven.be/seed-scores?beneficiaryId=883079&donorId=1379686&typeCategory=ncbiTaxonomyIds"
  ```
  <!--     # curl "http://localhost:1337/seed-scores?beneficiaryId=883079&donorId=1379686&typeCategory=ncbiTaxonomyIds" -->

  ```{important}
    The function returns **two pairs of seed scores**, in which:
    * the first genome provided is considered as $speciesA$ for the seed score indices,
    * and the second one as $speciesB$.
  ```
  <!-- In its current version, our API is not clear enough, and you need to remember that 
  the **first entry** considers the first genome as $speciesA$, 
  and the second genome as $speciesB$, 
  while in the **second entry** it is the other way around. 
  This will be fixed in a future release.  -->

### Get seed complements 

  Seed complements can be retrieved for pairs of 3 different categories of IDs: 

  - NCBI Taxonomy IDs (`type_of_ids: ncbiTaxonomyIds`)
  - NCBI Genome accession IDs (`type_of_ids: ncbiGenomeIds`) and
  - PATRIC IDs (`type_of_ids: patricGenomeIds`)

  The main route of this feature is 
  ```bash
  https://msysbio.gbiomed.kuleuven.be/seed-complements/<beneficiary_id>/<donor_id>/<type_of_ids>
  ```

  For example, one might get seed complements for a pair of NCBI Taxonomy IDs like this:

  ```bash
  curl "https://msysbio.gbiomed.kuleuven.be/seed-complements?beneficiaryId=1379686&donorId=883079&typeCategory=ncbiTaxonomyIds"
  ```
  <!--     #curl "http://localhost:1337/seed-complements?beneficiaryId=1379686&donorId=883079&typeCategory=ncbiTaxonomyIds"
  -->

  ```{note}
  > This route has no default for your id type! If not provided by the user, the API will fail.
  ```


## Common errors

  There are two common types of client errors on API calls:

  - 400 Bad requests.

  In this case, you most probably are asking a malformed request syntax, invalid request message framing, or deceptive request routing.
  Check again your query and make sure you are u sing the right syntax.

  - 404 Not found.

  In this case, the server cannot find the requested resource. 

  You can have such errors also in cases you are asking for a genome/species/pair of such that is not part of the *microbetagDB*. 

  Keep in mind that you can always contact us through our [Matrix community](https://matrix.to/#/#microbetagcommunity:matrix.org) for more.
