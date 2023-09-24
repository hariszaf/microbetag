---
layout: default
title: Modules
nav_order: 2
---


# Modules 
{: .no_toc }


Just the Docs has some specific configuration parameters that can be defined in your Jekyll site's _config.yml file.
{: .fs-6 .fw-300 }

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
- pathway complementarity annotations between taxa have been found co-correlated; both taxa ara considered as potential donor and beneficiary (see [Pathway complementarity: an example]() for more)
- Complementarity [3] and competition [4] scores between draft metabolic reconstructions of GTDB representative genomes, mapped to the input taxa using [PhyloMint](https://github.com/mgtools/PhyloMint) (see [Seed-based complementarity and competition scores]() for more) 




## Functional annotations


### Based on phenDB 



We provide a description of each feature abbreviation, based on those from the [phenDB group](https://phendb.org/reports/modeloverview). The annotation is referring to the species under study. For example `NOB` $\rightarrow$ *species under study is part of the clade of NOB*. 

{: .important-title }
> Abbreviations' descriptions
>
> | **Abbreviation** | Description |
> |`NOB`| part of the clade of nitrite-oxidizing bacteria |
> |`T3SS`|  expresses a Type III secretion system |
> |`T6SS`| species under study expresses a Type VI secretion system |
> |`aSaccharolytic` | species under study has an asaccharolytic lifestyle|
> |`aceticAcid` | species under study produces acetic acid |
> |`aob` | species under study is part of the clade of ammonia-oxidizing bacteria |
> |`autoCo2` | species under study is capable of growth with CO2 as sole carbon source |
> |`butanol` | species under study produces butanol |
> |`butyricAcid` | species under study produces butyric acid |
> |`dGlucose` | species under study utilizes D-glucose |
> |`dLacticAcid` | species under study produces D lactic acid |
> |`ethanol` | species  producing ethanol |
> |`fermentative` |  has a Fermentative lifestyle |
> |`fixingN2` |  capable of fixing N2 |
> |`formicAcid` |  producing formic_acid |
> |`gtdbId` | GTDB id of the genome used for this analysis |
> |`halophilic` | halophilic lifestyle |
> |`hydrogen` | producing hydrogen |
> |`indole` | producing indole |
> |`isobutyricAcid` | producing isobutyric_acid |
> |`isovalericAcid` | producing isovaleric_acid |
> |`lLacticAcid` | producing L_lactic_acid |
> |`methanotroph` |  capable of growth with methane as the sole carbon source |
> |`nonFermentative` | has a nonfermentative lifestyle |
> |`phototrophy` | is able of phototrophy |
> |`psychrophilic` | species  has a psychrophilic lifestyle |
> |`rAcetoin` |  is producing R_acetoin |
> |`saccharolytic` |  has a Saccharolytic lifestyle |
> `succinicAcid` |  is producing succinic_acid |
> |`sulfateReducer` | is capable of using sulfate as terminal electron acceptor |
> |`symbiont` |  has an obligate intracellular lifestyle |
> |`thermophilic` | species has a thermophilic lifestyle|
>







![complementarity in kegg example](../assets/images/kegg_example.png){: width=80% }






## References

[1] Tackmann, J., Rodrigues, J.F.M. and von Mering, C., 2019. Rapid inference of direct interactions in large-scale ecological networks from heterogeneous microbial sequencing data. Cell systems, 9(3), pp.286-296, DOI: [10.1016/j.cels.2019.08.002](https://doi.org/10.1016/j.cels.2019.08.002).
[2] Louca, S., Parfrey, L.W., Doebeli, M. (2016) - Decoupling function and taxonomy in the global ocean microbiome. Science 353:1272-1277, DOI: [10.1126/science.aaf4507](https://doi.org/10.1126/science.aaf4507).
[3] Levy, R., Carr, R., Kreimer, A., Freilich, S., Borenstein, E. "NetCooperate: a network-based tool for inferring host-microbe and microbe-microbe cooperation." BMC Bioinformatics, 2015.
[4] Kreimer, A., Doron-Faigenboim, A., Borenstein, E., Freilich, S. "NetCmpt: a network-based tool for calculating the metabolic competition between bacterial species." Bioinformatics, 2012.


