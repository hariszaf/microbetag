---
layout: default
title: Modules
nav_order: 2
has_children: true
permalink: /docs/modules
---


# Modules 

<!-- {: .no_toc } 
Just the Docs has some specific configuration parameters that can be defined in your Jekyll site's _config.yml file.
{: .fs-6 .fw-300 } 
-->

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

microbetag gets as input either a co-occurrence network or an abundance table where either [Silva](https://www.arb-silva.de) or [GTDB](https://gtdb.ecogenomic.org) taxonomies have been used. 
When an abundance table is provided, microbetag firsts builds a co-occurrence network using [FlashWeave](https://github.com/meringlab/FlashWeave.jl) [1].

microbetag consists of 4 major modules: 

- literature oriented taxa funcitonal annotation using [FAPROTAX](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php) [2]
- functional annotation using an updated, local instance of [phenDB](https://phendb.org) using all representative genomes of GTDB and [`phenotrex`](https://phenotrex.readthedocs.io/en/latest/usage.html)
- pathway complementarity annotations between taxa have been found co-correlated; both taxa ara considered as potential donor and beneficiary (see [Pathway complementarity: an example](#pathway-complementarity-an-example) for more)
- Complementarity [3] and competition [4] scores between draft metabolic reconstructions of GTDB representative genomes, mapped to the input taxa using [PhyloMint](https://github.com/mgtools/PhyloMint) (see [Seed-based complementarity and competition scores]() for more) 




## Functional annotations


### Based on FAPROTAX


FAPROTAX maps taxa (e.g. genera or species) to metabolic or other ecologically relevant functions based on the literature on cultured representatives. 
For a thorough description of that, you may have a look to its [documentation site](http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Instructions). 

The taxonomy assigned to each OTU/ASV (amplicon data) or bin (shotgun data) on the abundance table provided by the user, is mapped to a list of functions one can check [here](faprotax-functions.md).


As an example, here is how the FAPROTAX outcome looks like for the case of denitrification function: 


![faprotax example denitrif](../../assets/images/faprotax_denitrification.png)






### Based on phenDB 



Under the [Traits predicted based on phenDB models](phen-traits.md) tab, we provide a description of each feature abbreviation, based on those from the [phenDB group](https://phendb.org/reports/modeloverview). 


The annotation is referring to the species under study. 
Each trait gets a "Yes" or "No" decision along with an accurracy score. 
For example `NOB` : *species under study is part of the clade of NOB*. 

Here is an example of how 2 GTDB genomes look like: 

![phen traits example](../../assets/images/phen_traits_fmt.png)






## Pathway complementarity: an example

pathways found in both taxa of an association are further explored to check whether the processes of each of the two taxa are complementary denoting a positive interaction. Likewise, if the same processes are found to occur in both taxa, a negative interaction will be derived.

![complementarity in kegg example](../assets/images/kegg_example.png){: width=80% }






## References

[1] Tackmann, J., Rodrigues, J.F.M. and von Mering, C., 2019. Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell systems, 9(3), pp.286-296, DOI: [10.1016/j.cels.2019.08.002](https://doi.org/10.1016/j.cels.2019.08.002).

[2] Louca, S., Parfrey, L.W., Doebeli, M. (2016) - Decoupling function and taxonomy in the global ocean microbiome. Science 353:1272-1277, DOI: [10.1126/science.aaf4507](https://doi.org/10.1126/science.aaf4507).

[3] Levy, R., Carr, R., Kreimer, A., Freilich, S., Borenstein, E. "NetCooperate: a network-based tool for inferring host-microbe and microbe-microbe cooperation." BMC Bioinformatics, 2015.

[4] Kreimer, A., Doron-Faigenboim, A., Borenstein, E., Freilich, S. "NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species." Bioinformatics, 2012.

