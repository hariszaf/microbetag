# microbetag tests

```
tests/
  conftest.py     shared fixtures and skip markers
  unit/           Python only; run in CI on every push and pull request
  integration/    need external tools, conda environments or large data
  data/           inputs and reference outputs (see data/README.md)
```

## Setup

Install microbetag and its tools once, from the root of the repository
(add `--phenotrex`, `--modelseed`, `--dnngior` or `--all` for the optional environments):

```bash
bash install.sh
conda activate microbetag
```

`install.sh` installs the package in editable mode with the `test` extra, so pytest is included
and changes to the code are picked up without reinstalling.
If you only need the unit tests, `pip install -e ".[test]"` is enough.

## Running the tests

Run pytest from the **root of the repository**, where `pyproject.toml` is:

```bash
pytest                                    # everything; tests whose requirements are missing are skipped
pytest tests/unit                         # Python-only tests, what CI runs
pytest tests/integration                  # tests that need external tools or conda environments
pytest tests/integration/test_prodigal.py # a single file
pytest -k seed                            # tests whose name contains "seed"
pytest tests/unit/test_cli.py -v -s       # verbose, showing microbetag's log output
```

A test whose requirements are missing is **skipped** with a hint on how to install them, e.g.

```
SKIPPED tests/integration/test_phenotrex.py: conda env mtg-phenotrex not found; run `bash install.sh --phenotrex`
```

## Trying single steps by hand

The commands that the Nextflow pipeline uses can be run on the test data directly.
Run them in a scratch folder, so that their outputs stay out of the repository;
`MTG` points to the root of the repository:

```bash
MTG=$(pwd)                     # run this from the root of the repository
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

To run the whole pipeline, copy the example, since microbetag uses the folder of the
configuration file as its working directory:

```bash
cp -r $MTG/tests/data/microbetag ~/mtg-try/full && cd ~/mtg-try/full
tar xzf input_files.tar.gz && gunzip input_files/ko_merged.txt.gz
microbetag --config config_from_net.yml
```

Results are written to `mtg_complete_output/`. This configuration builds GEMs with CarveMe,
which needs a Gurobi license; set `seed_complementarity: value: false` in the YAML to run without it.

| Test | Needs |
|---|---|
| `unit/*` | `pip install ".[test]"`; `test_manta.py::test_run_manta` also needs manta (included in `bash install.sh`) |
| `integration/test_prodigal.py` | `bash install.sh` |
| `integration/test_kegg_annotation.py` | `bash install.sh` (uses the kofam subset in `ext_data/kofam_database`) |
| `integration/test_flashweave.py` | `bash install.sh` |
| `integration/test_carve.py` | `bash install.sh` and a Gurobi license |
| `integration/test_phenotrex.py` | `bash install.sh --phenotrex` |
| `integration/test_modelseed.py` | `bash install.sh --modelseed --dnngior`, a Gurobi license, and RASTtk for genomes |
| `integration/test_build_cx2.py` | expected to fail until `data/build_cx2/output_files` is regenerated |
| `integration/test_microbetag.py` | `bash install.sh --phenotrex`; complete run, takes a while |

## Writing tests

- Use plain pytest functions and fixtures; do not depend on the order tests run in.
- Never write into `tests/data`. Write outputs under `tmp_path`, or use the `data_copy("<name>")` fixture
  when microbetag writes next to its config file.
- Declare requirements with markers instead of failing:
  `@pytest.mark.requires_tools("prodigal")`, `@pytest.mark.requires_conda_env("mtg-phenotrex")`,
  `@pytest.mark.requires_data("microbetag/input_files")`, or the `gurobi_license` fixture.
- Assert on the results, not only that nothing crashed.
- `load_config(path_to_yaml)` from `conftest.py` builds microbetag's own `Config` from a test YAML.
  microbetag uses the folder of the YAML as its working directory, so keep a test's inputs next to it.
