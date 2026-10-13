


The `microbetag` folder is a `git submodule` of the `microbetag` stand-alone tool:

```{bash}
git submodule add --branch fix-phylomint https://github.com/hariszaf/microbetag.git microbetag
```


In the `mappings` folder, we keep taxonomy mappings (between NCBI, GTDB, Silva etc.), as well as, 
seed and non-seed sets of the GTDB GEMs reconstructed with ModelSEED and mappings 
between ModelSEED and KEGG ontologies. 



In the `.dockerignore` we ignore the `ext_data/` of the sub-module repo, where locally we have the `kofam` database (~6.8GB).
Make sure you mount the `kofam` database.



