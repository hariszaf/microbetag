"""
Shared fixtures for the microbetag tests.

- tests/unit:        single microbetag steps
- tests/integration: a complete microbetag run

Tests never write into tests/data: copy what you need with the `data_copy` fixture and
write outputs under pytest's `tmp_path`.

Options:
  --draw          draw the graphs of the seed-set tests (tests/unit/test_get_seeds.py)
  --graph-dir     where to save them (default: test_graphs)
"""

import shutil
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


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


def pytest_addoption(parser):
    parser.addoption(
        "--draw",
        action="store_true",
        default=False,
        help="Draw graph for the test.",
    )

    parser.addoption(
        "--graph-dir",
        action="store",
        default="test_graphs",
        help="Directory where graphs will be saved",
    )


@pytest.fixture
def draw_enabled(request):
    return request.config.getoption("--draw")


@pytest.fixture
def graph_dir(request):
    return Path(request.config.getoption("--graph-dir"))
