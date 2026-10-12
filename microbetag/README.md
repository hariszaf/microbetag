# microbetag: the Python package

Annotates microbial co-occurrence networks with phenotypic traits, and pathway and seed complementarities,
computed from the outputs of external tools (KEGG annotations, genome-scale models, phenotrex predictions).
It does not run the tools themselves: the [pipeline](../pipeline) does.

```
microbetag/
├── pyproject.toml   package metadata and commands
├── Dockerfile       the hariszaf/microbetag image
├── src/microbetag/  the package
└── tests/           unit and integration tests
```

Install, from this folder:

```bash
pip install .            # or ".[otf]" for the on-the-fly version (microbetagDB), ".[test]" for the tests
microbetag --template config.yml
microbetag --config config.yml
```

Run the tests:

```bash
pip install -e ".[test]"
pytest
```

Build the image:

```bash
docker build -t hariszaf/microbetag:<version> .
```

See the [repository's README](../README.md) and the [documentation](https://microbetag.readthedocs.io/).
