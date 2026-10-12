"""Network inference with FlashWeave, with and without a metadata file."""

import os

import pytest

from conftest import load_config
from microbetag.tools import run_flashweave

pytestmark = pytest.mark.requires_tools("julia")


@pytest.fixture
def config(data_copy):
    pytest.importorskip("julia", reason="pyjulia not installed; run `bash install.sh`")
    return load_config(data_copy("flashweave") / "config_test_fw.yml")


def test_run_flashweave(config):
    config.metadata      = "false"
    config.metadata_file = None

    run_flashweave(config=config)

    assert os.path.getsize(config.network) > 0


def test_run_flashweave_with_metadata(config):
    config.network = os.path.join(config.outdir, "network_metadata.edgelist")

    run_flashweave(config=config)

    assert os.path.getsize(config.network) > 0
