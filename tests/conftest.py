"""
Shared fixtures and markers for the microbetag tests.

- tests/unit:        Python only, run in CI on every push/PR.
- tests/integration: need external tools (install.sh), conda environments or large data;
                     they are skipped with a hint when what they need is missing.

Tests never write into tests/data: copy what you need with the `data_copy` fixture and
write outputs under pytest's `tmp_path`.
"""

import shutil
from pathlib import Path

import pytest

DATA_DIR  = Path(__file__).parent / "data"
KOFAM_DIR = Path(__file__).parent.parent / "ext_data" / "kofam_database"

# Tools that install.sh can't provide
TOOL_HINTS = {
    "rast-create-genome": "install the BV-BRC command-line tools (RASTtk): "
                          "https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html",
}


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_tools(*names): skip unless the external tools are on PATH (bash install.sh)",
    )
    config.addinivalue_line(
        "markers",
        "requires_conda_env(name): skip unless the conda environment exists (bash install.sh --<option>)",
    )
    config.addinivalue_line(
        "markers",
        "requires_data(*paths): skip unless the paths exist under tests/data",
    )


def pytest_runtest_setup(item):
    for marker in item.iter_markers("requires_tools"):
        missing = [t for t in marker.args if shutil.which(t) is None]
        if missing:
            hint = TOOL_HINTS.get(missing[0], "run `bash install.sh` and `conda activate microbetag`")
            pytest.skip(f"{', '.join(missing)} not found; {hint}")

    for marker in item.iter_markers("requires_conda_env"):
        from microbetag import deps
        name = marker.args[0]
        if name not in deps._conda_envs():
            pytest.skip(f"conda env {name} not found; run `bash install.sh {deps.CONDA_ENVS[name][1]}`")

    for marker in item.iter_markers("requires_data"):
        missing = [p for p in marker.args if not (DATA_DIR / p).exists()]
        if missing:
            pytest.skip(f"test data missing: {', '.join(missing)}")


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Read-only test data (tests/data)."""
    return DATA_DIR


@pytest.fixture
def data_copy(tmp_path):
    """
    Copy tests/data/<name> into a temporary directory and return its path.

    Use it for tests whose outputs are written next to their inputs or config file,
    so that a run never changes tracked files.
    """
    def _copy(name: str) -> Path:
        dest = tmp_path / name
        shutil.copytree(DATA_DIR / name, dest)
        return dest

    return _copy


@pytest.fixture(scope="session")
def gurobi_license():
    """Skip unless gurobipy can create a model, i.e. a valid Gurobi license is installed."""
    gp = pytest.importorskip("gurobipy", reason="gurobipy not installed; pip install 'microbetag[genres]'")
    try:
        gp.Model().dispose()
    except gp.GurobiError as e:
        pytest.skip(f"no valid Gurobi license ({e}); see https://www.gurobi.com/academia/")


def load_config(config_file: Path, **overrides):
    """
    microbetag's own Config from a test YAML; paths in it are relative to the YAML's folder.
    `overrides` replace top-level YAML entries, e.g. kofam_database={"dir_path": "/abs/path"}.
    """
    import yaml
    from microbetag.config import Config

    with open(config_file) as f:
        conf = yaml.safe_load(f)
    conf.update(overrides)
    return Config(conf, str(config_file))
