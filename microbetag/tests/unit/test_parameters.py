"""
The parameters of the configuration (templates/parameters.yaml): the template, --describe, and that the code reads
only parameters described there.
"""

import re
from pathlib import Path

import pytest
import yaml

from microbetag import parameters
from microbetag.microbetag import main

SRC = Path(parameters.__file__).parent


def test_every_parameter_is_described():
    for name, p in parameters.all_parameters().items():
        assert p.get("description"), name
        assert sum(f in p for f in parameters.FIELDS) == 1, name


def test_template_has_every_parameter():
    spec     = parameters.load_parameters()
    template = yaml.safe_load(parameters.render_template(spec))
    params   = {n: p for s in spec["sections"] for n, p in s["parameters"].items()}
    assert list(template) == list(params)
    for name, p in params.items():
        assert template[name] == {parameters.field(p): parameters.template_value(p)}
    assert template["seed_complementarity"] == {"value": True}
    assert template["metadata_file"] == {"file_path": None}


def test_defaults():
    # What Config uses for a parameter left out, or set to null
    assert parameters.get_value({}, "threads") == 2
    assert parameters.get_value({"threads": None}, "threads") == 2
    assert parameters.get_value({"threads": {"value": None}}, "threads") == 2
    assert parameters.get_value({"threads": {"value": 8}}, "threads") == 8
    assert parameters.get_value({}, "seed_complementarity") is False
    assert parameters.get_value({}, "pc_percent") == 1


def test_docs_tables():
    md = parameters.render_markdown()
    for name in parameters.all_parameters():
        assert f"| `{name}` |" in md
    assert "| `seed_complementarity` | `value` | `false` | `true` |" in md


def test_code_reads_only_described_parameters():
    read = set()
    for py in SRC.rglob("*.py"):
        text  = py.read_text()
        read |= set(re.findall(r'(?:conf|yaml)\.get\("(\w+)"', text))
        read |= set(re.findall(r'get_value\([\w.]+, "(\w+)"', text))
    read -= {"__config_dir__"}
    assert read and read <= set(parameters.all_parameters())


def test_unknown_parameters():
    assert parameters.unknown_parameters({"edge_list": {}, "edgelist": {}, "__config_dir__": "."}) == ["edgelist"]


def test_cli_template_and_describe(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["--write-template", str(tmp_path / "config.yml")])
    assert "metadata_file" in yaml.safe_load((tmp_path / "config.yml").read_text())
    capsys.readouterr()

    with pytest.raises(SystemExit):
        main(["--describe", "metadata_file"])
    out = capsys.readouterr().out
    assert out.startswith("metadata_file  (file_path, default: not set)") and "metavar" in out

    with pytest.raises(SystemExit, match="Unknown parameter"):
        main(["--describe", "metadata"])
