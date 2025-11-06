
Two main workflows: 

1. catalogue precomputation ([`precalc.nf`](./precalc.nf))
2. Network annotation using a selected precomputation set ([`mtg.nf`](./mtg.nf))

You may run modules under `modules/` individually, for example:

```
nextflow run modules/phenotrex/main.nf -params-file params/phenotrex.yaml -entry PHENOTREX

# or

nextflow run modules/phenotrex/main.nf -config modules/phenotrex/phenotrex.config -entry PHENOTREX

```

> In case of the [`mtg_annotate`](./modules/mtg_annotate/mtg_annotate.nf), use only the `-config` case. 

You can edit the parameters file you will find under [`params`](./params) accordingly, to better tune modules based on the needs of your data.




## Notes

Inside a workflow/process scope, `def x = ...` is interpreted as:

> Declare a local variable named `ko_list_ch` AND also shadow/override `Channel` symbol resolution.

Then `Channel.fromPath` is misinterpreted as:

> There is a variable named `Channel` defined here



# Since we can have up to 5 parallel sessions with the WLS license of Gurobi,
# which is the only way to go on with a Docker container, we set by defalt max_forks to 5.
# Then, our input files will be split in chunks of max_forks size for parallel processing.


```
nextflow run workflows/precalc.nf -params-file workflows/precalc.yaml  -entry MICROBETAG_PRECALC
```


### Modules 

```
nextflow run modules/faprotax/main.nf -params-file modules/faprotax/faprotax.yaml 

nextflow run modules/flashweave/main.nf -params-file modules/flashweave/flashweave.yaml 

nextflow run modules/gem_recon/main.nf -params-file modules/gem_recon/gem_recon.yaml 

# Before the kofam 
nextflow run modules/prodigal/main.nf -params-file modules/prodigal/prodigal.yaml 

nextflow run modules/kofam/main.nf -params-file modules/kofam/kofam.yaml 

nextflow run modules/manta/main.nf -params-file modules/manta/manta.yaml 

nextflow run modules/pathway_compl/main.nf -params-file modules/pathway_compl/pathway_compl.yaml 

nextflow run modules/phenotrex/main.nf -params-file modules/phenotrex/phenotrex.yaml 

nextflow run modules/seed_compl/main.nf -params-file modules/seed_compl/seed_compl.yaml
```

Last, but not least, the most tricky of the modules, the network annotation ! 



### Subworkflows

```
nextflow run subworkflows/pathway_complementarity/main.nf \
    -params-file subworkflows/pathway_complementarity/pc.yaml \
    -entry PATHWAY_COMPLEMENTARITY
```

```
nextflow run subworkflows/seed_complementarity/main.nf \
    -params-file subworkflows/seed_complementarity/sc.yaml \
    -entry SEED_COMPLEMENTARITY
```


### Workflows

```
nextflow run workflows/precalc.nf \
    -params-file workflows/precalc.yaml \
    -entry MICROBETAG_PRECALC
```


```
nextflow run workflows/net_annotate.nf \
    -params-file workflows/net_annotate.yaml \
    -entry MICROBETAG_ANNOTATE 
```




