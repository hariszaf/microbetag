from pathlib import Path

import pytest


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
