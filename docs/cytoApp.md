---
layout: default
title: Cytoscape App
nav_order: 3
---

# Cytoscape App
{: .no_toc }

---




![microbetag CyApp](../assets/images/cyApp.png){: width=25% }



## How to run

### starting from an abundance table 



[Download ][1] dataset used for this example.

[1]:{{ site.url }}/microbetag/download/dada2_use_case.tsv



microbetag requires for a 7-level taxonomy scheme; 
for example 
```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;Caldanaerobius;Caldanaerobius polysaccharolyticus
```

in case an entry reaches only to a higher taxonomic level, microbetag fills the entry with NA values

for example

```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae
```

would become

```bash
Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;NA;NA;NA
```


{: .important-title}
> Curate your taxonomies! 
> 
> If you have a taxonomy that "skips" a level,  




if you already have a network, then you need to provide both the network and the abundance table ! 


fuzzywuzzy uses a threshold - that is set relatively high (90) so there are no false positives. 
in order not to loose species level annotations, please have a look so you do not any unessary characters on your taxonomies, e.g. `[Salmonella] infantis` 
would get a lower score that `Salmonella infantis`, so removing `[` and `]` characters would benefit. 





To check 


Sensitive vs fast mode
• Implementation of
conditional independence:
– Sensitive mode: partial
correlations on abundances,
assumes multivariate normal
distribution (weak assumption)
– Fast mode: mutual information
on presence/absences


HE mode
• FlashWeave can optionally ignore
zeros (‘structural zeros’) to deal
with heterogeneous samples 


multi-habitat or -protocol data sets with ideally at least thousands of samples;

sensitive=false for faster, but more coarse-grained associations

## example
Here we need a step-by-step with screenshots and/or videos of the features of 

