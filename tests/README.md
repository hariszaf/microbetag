# microbetag tests

```
tests/
  conftest.py     shared fixtures and skip markers
  unit/           Python only; run in CI on every push and pull request
  integration/    need external tools, conda environments or large data
  data/           inputs and reference outputs (see data/README.md)
```

## Running

```bash
conda activate microbetag
pytest                      # everything; tests whose requirements are missing are skipped
pytest tests/unit           # what CI runs
pytest tests/integration -k prodigal
```

A test whose requirements are missing is **skipped** with a hint on how to install them, e.g.

```
SKIPPED tests/integration/test_phenotrex.py: conda env mtg-phenotrex not found; run `bash install.sh --phenotrex`
```

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
