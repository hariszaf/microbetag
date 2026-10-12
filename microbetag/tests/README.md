# microbetag tests

```
tests/
  conftest.py     shared fixtures
  unit/           single microbetag steps
  integration/    a complete microbetag run
  data/           inputs and reference outputs (see data/README.md)
```

The tests need only the package; the external tools are tested in
[microbetag-nf](https://github.com/msysbio/microbetag-nf).

## Running the tests

Install the package with the `test` extra, then run pytest from the **`microbetag/` folder** of the repository,
where `pyproject.toml` is:

```bash
pip install -e ".[test]"
pytest                                    # everything
pytest tests/unit                         # single steps
pytest tests/integration                  # the complete run
pytest -k seed                            # tests whose name contains "seed"
pytest tests/unit/test_cli.py -v -s       # verbose, showing microbetag's log output
```

## Trying single steps by hand

The commands that the Nextflow pipeline uses can be run on the test data directly.
Run them in a scratch folder, so that their outputs stay out of the repository;
`MTG` points to the package folder (`microbetag/` of the repository):

```bash
MTG=$(pwd)                     # run this from microbetag/, where pyproject.toml is
mkdir -p ~/mtg-try && cd ~/mtg-try

# Pathway complementarities from a 3-column KO table
microbetag-pcompl alts.json pc.json $MTG/tests/data/path_compl/input_files/ko_merged_7bins.txt.gz false 2

# Seed set, non-seed set and confidence of a GEM
microbetag-seed-sets $MTG/tests/data/seed_compl/input_files/bin_101.xml BiGG bin_101_sets.json

# Seed scores and complements of one species against the others
microbetag-seed-compls bin_101 \
    $MTG/tests/data/seed_compl/output_files/nonSeedSetDic.json \
    $MTG/tests/data/seed_compl/output_files/confidenceDic.json \
    false scores.tsv compls.json
```

Every command has `--help`, e.g. `microbetag-annotate --help`.

To run the whole of microbetag, copy the example, since microbetag uses the folder of the
configuration file as its working directory:

```bash
cp -r $MTG/tests/data/microbetag ~/mtg-try/full && cd ~/mtg-try/full
microbetag --config config.yml
```

Results are written to `output/`: the pathway and seed complementarities and the annotated network
(`mtag_net_<date>.cx2`), which you can open in Cytoscape.

## Writing tests

- Use plain pytest functions and fixtures; do not depend on the order tests run in.
- Never write into `tests/data`. Write outputs under `tmp_path`, or use the `data_copy("<name>")` fixture
  when microbetag writes next to its config file.
- Assert on the results, not only that nothing crashed.
