# Test data

| Folder | Used by | Contents |
|---|---|---|
| `microbetag/` | `integration/test_microbetag.py` | Hessler data set: 7 bins, their co-occurrence network, abundance table and the outputs of the external tools (KO annotation, GEMs, phenotrex predictions, FAPROTAX sub-tables, manta clustered network), with `config.yml` |
| `path_compl/` | `unit/test_path_compl.py`, `unit/test_cli.py` | KO annotation of 7 bins (`ko_merged_7bins.txt.gz`) |
| `seed_compl/` | `unit/test_seed_compl.py`, `unit/test_cli.py` | two GEMs and their reference seed sets, confidences and scores |
| `seed_compl/iML1515.xml.gz` | `unit/test_seed_inference.py` | *E. coli* K-12 MG1655, the curated iML1515 model of [BiGG](http://bigg.ucsd.edu/models/iML1515) (Monk et al., 2017): the positive control of the seed sets |
| `manta/` | `unit/test_manta.py` | networks with an abundance table or a sequence-to-taxonomy map, to build manta's input |

The inputs used to test the external tools (Prodigal, kofam, phenotrex, FlashWeave, CarveMe,
ModelSEEDpy, FAPROTAX) moved with those steps to the pipeline (`pipeline/tests`).
