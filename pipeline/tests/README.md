# Test data

```
tests/
  data/        inputs
  expected/    reference outputs; runs never write here
```

| Folder | Contents | Used by |
|---|---|---|
| `data/genomes` | the example genome catalogue (3 genomes) | `tests/params/precalc.yaml`, prodigal, phenotrex |
| `data/abd_data` | example abundance tables and metadata | `tests/params/annotate.yaml`, FAPROTAX, FlashWeave, manta |
| `data/prep` | an amplicon table with sequences, its metadata, and the output of the former `microbetag_prep` tool | `tests/params/prep.yaml`, taxonomy |
| `data/media` | gapseq's glucose minimal medium (`MM_glu.csv`, ModelSEED ids, from gapseq's `dat/media`), and the same in BiGG ids (`MM_glu_bigg.csv`, through MetaNetX, and chromate as cro4) | gap-filling: ModelSEEDpy, gapseq, DNNGIOR |
| `data/kofam_database` | `ko_list`, and a subset of the KOfam profiles for tests (`profiles_test.tar.gz`; unpack it first) | kofam |
| `data/<tool>` | inputs of the single-tool tests, moved from the microbetag repository (see `data/README.md`) | |
| `expected/precalc` | the pre-calculations of `precalc.nf` on `data/genomes` | `tests/params/annotate.yaml`, single modules |
| `expected/annotate` | the annotated network of `annotate.nf` with `expected/precalc` and `data/abd_data/abund.tsv` | |

The params files of the tests (`params/*.yaml`: the workflows and sub-workflows; `../modules/*/*.yaml`) read from here
and write to `tests/results/`, which is not tracked. To check a change, run a workflow and compare its `tests/results/` with
`tests/expected/`, e.g.:

```bash
nextflow run main.nf -params-file tests/params/annotate.yaml
# same edges as the reference network; the weights may differ in their last decimals
diff <(cut -f1,2 tests/results/annotate/flashweave/fw_net.txt | sort) <(cut -f1,2 tests/expected/annotate/flashweave/fw_net.txt | sort)
```

Databases you download (e.g. the GTDB 16S reference for `prep`) go in `databases/`, which is not tracked either.
