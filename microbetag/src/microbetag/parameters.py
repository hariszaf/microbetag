"""
The parameters of the microbetag configuration, from templates/parameters.yaml: the one place where they are
described. The configuration template, `microbetag --describe` and the check for unknown parameters come from it.
"""

import textwrap
from importlib import resources

import yaml

FIELDS = ("value", "file_path", "dir_path")
WIDTH  = 110


def load_parameters() -> dict:
    return yaml.safe_load(resources.files("microbetag").joinpath("templates", "parameters.yaml").read_text())


def all_parameters(spec: dict = None) -> dict:
    """{name: its entry}, those of the template and the internal ones."""
    spec = spec or load_parameters()
    params = {name: p for section in spec["sections"] for name, p in section["parameters"].items()}
    params.update(spec.get("internal", {}))
    return params


def field(param: dict) -> str:
    """The key under which the configuration has the parameter's value (value, file_path or dir_path)."""
    return next(f for f in FIELDS if f in param)


def _comment(text: str, width: int = WIDTH) -> str:
    lines = []
    for paragraph in text.strip().split("\n\n"):
        lines += textwrap.wrap(" ".join(paragraph.split()), width - 2) or [""]
        lines.append("")
    return "\n".join(f"# {line}".rstrip() for line in lines[:-1])


def render_template(spec: dict = None) -> str:
    """The configuration template: every parameter, with its default, and its description as a comment."""
    spec = spec or load_parameters()
    rule = "# " + "-" * (WIDTH - 2)
    out  = [f"# {spec['title']} (microbetag --template config.yml)", "#", _comment(spec["intro"]), ""]
    for section in spec["sections"]:
        out += [rule, f"# {section['title']}", rule, ""]
        for name, p in section["parameters"].items():
            desc = p["description"] + (" Required." if p.get("required") else "")
            if p.get("example") is not None:
                desc += f" E.g. {p['example']}"
            key, default = field(p), p[field(p)]
            out += [_comment(desc), yaml.safe_dump({name: {key: default}}, sort_keys=False).rstrip(), ""]
    return "\n".join(out)


def describe(name: str = None) -> str:
    """The description of a parameter, or of all of them."""
    params = all_parameters()
    if name is not None and name not in params:
        raise KeyError(f"Unknown parameter: {name}. Known: {', '.join(params)}")
    lines = []
    for n, p in params.items():
        if name is None or n == name:
            default = p[field(p)]
            lines.append(f"{n}  ({field(p)}, default: {default if default is not None else 'not set'})")
            lines.append(textwrap.indent(textwrap.fill(" ".join(p["description"].split()), WIDTH - 4), "    "))
    return "\n".join(lines)


def unknown_parameters(conf: dict) -> list:
    """Parameters of a configuration that microbetag does not know (e.g. a typo)."""
    known = set(all_parameters())
    return sorted(k for k in conf if not k.startswith("__") and k not in known)
