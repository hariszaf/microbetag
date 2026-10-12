"""
KEGG ORTHOLOGY annotation with hmmsearch against a small subset of the kofam profiles
(ext_data/kofam_database/profiles_test.tar.gz and ko_list_test), on proteins predicted by Prodigal.
"""

import os
import shutil
import tarfile
from pathlib import Path

import pytest

from microbetag.tools import kegg_annotation, run_prodigal
from microbetag.utils import bin_kos_to_file, ko_list_parser, merge_ko

pytestmark = pytest.mark.requires_tools("prodigal", "hmmsearch")

KOFAM_TEST = Path(__file__).resolve().parents[2] / "ext_data" / "kofam_database"


@pytest.fixture
def kofam_db(tmp_path):
    """A kofam database folder (ko_list + profiles/) with the test subset."""
    db = tmp_path / "kofam"
    db.mkdir()
    shutil.copy(KOFAM_TEST / "ko_list_test", db / "ko_list")
    with tarfile.open(KOFAM_TEST / "profiles_test.tar.gz") as tar:
        tar.extractall(db)
    (db / "profiles_test").rename(db / "profiles")
    return db


def test_kegg_annotation(tmp_path, data_dir, kofam_db):
    bin_id = "bin_101"
    run_prodigal(str(data_dir / "prodigal" / "input_files" / f"{bin_id}.fa"), bin_id, str(tmp_path))

    hmmout_dir = tmp_path / "hmmout"
    bin_kos_dir = hmmout_dir / bin_id
    bin_kos_dir.mkdir(parents=True)

    kegg_annotation(
        faa=str(tmp_path / f"{bin_id}.faa"),
        basename=bin_id,
        out_dir=str(hmmout_dir),
        db_dir=str(kofam_db),
        ko_dic=ko_list_parser(str(kofam_db / "ko_list")),
        threads=2,
    )
    bin_kos_to_file(hmmout_dir=str(bin_kos_dir), bin_id=bin_id)

    ko_merged = tmp_path / "ko_merged.txt"
    merge_ko(str(hmmout_dir), str(ko_merged))

    assert os.path.getsize(ko_merged) > 0
    lines = ko_merged.read_text().splitlines()
    assert any(line.startswith(bin_id) and "\tK" in line for line in lines[1:])
