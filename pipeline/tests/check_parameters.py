#!/usr/bin/env python3
"""
Check parameters.yaml against the pipeline:

- every parameter has a description, its workflows, and a value or parameters (a group);
- every params.<name> of the workflows (main.nf, workflows/, and the modules they include) is in parameters.yaml,
  but for those only the standalone workflow blocks of the modules read (MODULE_ONLY);
- each test params file (tests/params/<workflow>.yaml) sets only known parameters.

Run from pipeline/:  python3 tests/check_parameters.py
"""

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent

# Read only by the workflow blocks of modules run on their own (nextflow run modules/<module>/main.nf), with the
# params files of the modules
MODULE_ONLY = {"genome_taxonomy", "precalc", "abundance_table", "cyjs_net", "faprotax_sub_tables", "manta_net"}

# Not parameters: methods of the params map
NOT_PARAMS = {"containsKey"}


def check_entry(name, p, errors, top=True):
    if not p.get("description"):
        errors.append(f"{name}: no description")
    if top and "workflows" not in p and name not in SPEC.get("internal", {}):
        errors.append(f"{name}: no workflows")
    if "parameters" in p:
        for k, v in p["parameters"].items():
            check_entry(f"{name}.{k}", v, errors, top=False)
    elif "value" not in p:
        errors.append(f"{name}: no value (default)")


def unknown_in(group, value, prefix):
    found = []
    for k, v in (value or {}).items():
        p = group.get("parameters", {}).get(k)
        if p is None:
            found.append(f"{prefix}.{k}")
        elif "parameters" in p and isinstance(v, dict):
            found += unknown_in(p, v, f"{prefix}.{k}")
    return found


SPEC = yaml.safe_load((ROOT / "parameters.yaml").read_text())


def main():
    errors = []
    params = {n: p for s in SPEC["sections"] for n, p in s["parameters"].items()}
    params.update(SPEC.get("internal", {}))

    for name, p in params.items():
        check_entry(name, p, errors)
        for wf in p.get("workflows", []):
            if wf not in SPEC["workflows"]:
                errors.append(f"{name}: unknown workflow {wf}")

    sources = [ROOT / "main.nf", *ROOT.glob("workflows/**/*.nf"), *ROOT.glob("modules/**/*.nf")]
    read = set()
    for nf in sources:
        read |= set(re.findall(r"(?<![\w/])params\.([A-Za-z_]\w*)", nf.read_text()))
    for name in sorted(read - set(params) - MODULE_ONLY - NOT_PARAMS):
        errors.append(f"params.{name} is read by the pipeline but not in parameters.yaml")

    for wf in SPEC["workflows"]:
        f = ROOT / "tests" / "params" / f"{wf}.yaml"
        given = yaml.safe_load(f.read_text()) or {}
        for name, value in given.items():
            if name not in params:
                errors.append(f"{f.relative_to(ROOT)}: unknown parameter {name}")
            elif "parameters" in params[name] and isinstance(value, dict):
                errors += [f"{f.relative_to(ROOT)}: unknown parameter {u}" for u in unknown_in(params[name], value, name)]

    for e in errors:
        print(e)
    print(f"{len(params)} parameters, {len(errors)} problems")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
