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

microbetag requires for a 7-level taxonomy scheme; 
for example 

Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;Caldanaerobius;Caldanaerobius polysaccharolyticus

in case an entry reaches only to a higher taxonomic level, microbetag fills the entry with NA values

for example

Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;

would become

Bacteria;Firmicutes;Thermoanaerobacteria;Thermoanaerobacterales;Thermoanaerobacteraceae;NA;NA;NA


{: .important-title}
> Curate your taxonomies! 
> 
> If you have a taxonomy that "skips" a level,  







