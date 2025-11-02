
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



