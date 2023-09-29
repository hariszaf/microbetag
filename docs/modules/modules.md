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

Once a network is availalbe, microbetag identifies the taxonomic level that has been assigned to each entry, for example 
`D_0__Bacteria; D_1__Firmicutes; D_2__Clostridia; D_3__Clostridiales; D_4__Ruminococcaceae; D_5__uncultured; D_6__uncultured rumen bacterium`
has reached the family level, while
`D_0__Bacteria; D_1__Actinobacteria; D_2__Coriobacteriia; D_3__Coriobacteriales; D_4__Coriobacteriaceae; D_5__Collinsella; D_6__uncultured bacterium`
is at the genus level.

The network annotation consists of 4 major modules: 

- **literature oriented** taxa funcitonal annotation using [**FAPROTAX**](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php) [2]
- **genomic oriented** functional annotation using an updated, local instance of [**phenDB**](https://phendb.org) using all representative genomes of GTDB and [`phenotrex`](https://phenotrex.readthedocs.io/en/latest/usage.html)
- **pathway complementarity** annotations between taxa have been found co-correlated; both taxa ara considered as potential donor and beneficiary (see [Pathway complementarity: an example](#pathway-complementarity-an-example) for more)
- **complementarity** [3] and **competition** [4] **seed scores** between draft metabolic reconstructions of GTDB representative genomes, mapped to the input taxa using [**PhyloMint**](https://github.com/mgtools/PhyloMint) (see [Seed-based complementarity and competition scores]() for more) 

FAPROTAX returns annotations for taxa (nodes) that they have been taxonomically annotated at the family or even at the order level in some cases, while the other 3 annotation types 
return annotations only at the strain and the species level.
**Nodes** that have species or strain taxonomic annotation are mapped to their closest representative GTDB genomes and based on those, they get phenDB-like functional annotations. 
**Edges** linking nodes that have been assigned at the species or strain level, i.e. both nodes of the association have a species/strain taxonomic annotation, are annotated using 
the pathway complementarity and the seed scores approaches. 

Below, you will find further background and examples of each annotation type. 


## Functional annotations


### Based on FAPROTAX


[**FAPROTAX**](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php) [2] maps taxa (e.g. genera or species) 
to metabolic or other ecologically relevant functions based on the literature on cultured representatives. 
It currently comprises more than 7600 annotation rules, covering ~4700 prokaryotic clades. 
Each annotation rule comes with literature citations and can thus be independently verified.
16S rRNA oriented approaches (e.g., PICRUSt, Tax4Fun etc) estimate community gene content based on available sequenced genomes.
Contrary, FAPROTAX estimates metabolic phenotypes based on experimental evidence.

The taxonomy assigned to each OTU/ASV (amplicon data) or bin (shotgun data) on the abundance table provided by the user, is mapped to a list of functions one can check [here](faprotax-functions.md).

As an example, here is how the FAPROTAX outcome looks like for the case of **denitrification** function: 


![faprotax example denitrif](../../assets/images/faprotax_denitrification.png)

All ASVs present in this file are related to the **denetrification** function.
Numbers represent the ASV abundance in each sample. 

microbetag runs FAPROTAX agains the abundance table and parses the subtables ([`seqId_faprotax_functions_assignment`](https://github.com/msysbio/microbetagApp/blob/main/services/web/microbetag/scripts/utils.py#L181)) 
to annotate each node with the corresponding function. 



### Based on phenDB 



[`phenotrex`](https://phenotrex.readthedocs.io/en/latest/usage.html) enables phenotypic trait prediction on user's metagenomic genomes/bins etc.

Phenotrex classifiers were re-trained using the genomes provided by phenDB for each model. 
For example, for the acetic acid production case, the [corresponding webpage of phenDB](https://phendb.org/reports/modeldetails?model_id=16) pointed to the set of genomes that had been originally used. 
These genomes were recovered and the classifiers were re-trained using 



Under the [Traits predicted based on phenDB models](phen-traits.md) tab, we provide a description of each feature abbreviation, based on those from the [phenDB group](https://phendb.org/reports/modeloverview). 


The annotation is referring to the species under study. 
Each trait gets a "Yes" or "No" decision along with an accurracy score. 
For example `NOB` : *species under study is part of the clade of NOB*. 

Here is an example of how 2 GTDB genomes look like: 

![phen traits example](../../assets/images/phen_traits_fmt.png)










## Pathway complementarity: an example

pathways found in both taxa of an association are further explored to check whether the processes of each of the two taxa are complementary denoting a positive interaction. Likewise, if the same processes are found to occur in both taxa, a negative interaction will be derived.

![complementarity in kegg example](../../assets/images/kegg_example.png){: width=80% }








## Seed scores based on genome-scale draft reconstructions 

Based on Borenstein et al () a metabolic network's “seed set”—the set of compounds that, based on the network topology, are exogenously acquired"


The fraction of A’s seed set that is also in B’s seed set, normalized by the weighted sum of the confidence score


As described in the PhyloMint papaer: the Complementarity Index is calculated as the fraction of the seed set of the genome-scale reconstruction of species A, that is found within B’s metabolic network but not part of B’s seed set, normalized by the number of A’s seed set in B’s entire metabolic network [34, 45]. MIComplementarity represents the potential for A’s to utilize the potential metabolic output of B



![seed concept](../../assets/images/seed_concept_example.png)

Node A is a seed, as it cannot be activated by any other node in the network.

In contrast, nodes F, G, and H are all marked as seeds, but can be seen to be interdependent: Activating one of these nodes would activate the rest, but at least one must be active to activate the rest. These nodes form a "seed group".

To quantify the relevance of each identified seed, we assign each seed a "confidence level", ranging from 0 to 1. A confidence level of 0 would correspond to a non-seed node, while a 1 would correspond to a seed that cannot be activated by another node. Seeds which belong to a seed group with more than 1 seed are given a fractional confidence level, the inverse of the number of seeds in the group. Nodes F, G, and H would then each have a confidence level of 1/3.



![netcooperate seed example](../../assets/images/seed_network_example.png)

In this example, we have a pair of simple networks. The top network has two seed groups, one of size one and the other size three (circled with a double line). The bottom network has three seed groups, two singletons, and one two-node seed group (circled with a double line). Seeds of one network present in the other are circled (solid lines: present, but not a seed of, the second network; dashed lines: present as a seed of the second network). The Biosynthetic Support Score of the top network on the bottom (i.e. treating the top as a parasite of the bottom) is 1.0. Note that all seed groups (but not all seeds) are present in the bottom network, and that node F is a seed of both networks. The Metabolic Complementarity Index of the top network on the bottom (i.e. treating both networks as co-occurring microbes) is 0.5. Because node F is a seed of the bottom network, it is not complementary to the top network’s seed set.




## References

[1] Tackmann, J., Rodrigues, J.F.M. and von Mering, C., 2019. Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell systems, 9(3), pp.286-296, DOI: [10.1016/j.cels.2019.08.002](https://doi.org/10.1016/j.cels.2019.08.002).

[2] Louca, S., Parfrey, L.W., Doebeli, M. (2016) - Decoupling function and taxonomy in the global ocean microbiome. Science 353:1272-1277, DOI: [10.1126/science.aaf4507](https://doi.org/10.1126/science.aaf4507).

[3] Levy, R., Carr, R., Kreimer, A., Freilich, S., Borenstein, E. "NetCooperate: a network-based tool for inferring host-microbe and microbe-microbe cooperation." BMC Bioinformatics, 2015.

[4] Kreimer, A., Doron-Faigenboim, A., Borenstein, E., Freilich, S. "NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species." Bioinformatics, 2012.


