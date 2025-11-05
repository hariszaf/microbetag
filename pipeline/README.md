
Two main workflows: 

1. catalogue precomputation ([`precalc.nf`](./precalc.nf))
2. Network annotation using a selected precomputation set ([`mtg.nf`](./mtg.nf))

You may run modules under `modules/` individually, for example:

``
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
