"""
The parameters of the microbetag configuration, from templates/parameters.yaml: the one place where they are
described. The configuration template, `microbetag --describe`, the defaults of Config, its check for unknown
parameters and the docs page of the parameters come from it.

This module imports nothing of microbetag, so that docs/conf.py can load it on its own.
"""

import textwrap
from importlib import resources

import yaml

FIELDS = ("value", "file_path", "dir_path")
WIDTH  = 110


def load_parameters(path: str = None) -> dict:
    if path is not None:
        with open(path) as f:
            return yaml.safe_load(f)
    return yaml.safe_load(resources.files("microbetag").joinpath("templates", "parameters.yaml").read_text())


_PARAMETERS = None


def default(name: str):
    """What microbetag uses for a parameter that the configuration leaves out."""
    global _PARAMETERS
    if _PARAMETERS is None:
        _PARAMETERS = all_parameters()
    p = _PARAMETERS[name]
    return p[field(p)]


def get_value(conf: dict, key: str):
    """
    The 'value' of a parameter of the configuration (conf[key]['value']), or its default if the configuration leaves
    it out or sets it to null.
    """
    value = (conf.get(key) or {}).get("value")
    return default(key) if value is None else value


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


def template_value(param: dict):
    """What the template writes for a parameter: its `template` value if any, else its default."""
    return param["template"] if "template" in param else param[field(param)]


def _show(value) -> str:
    return "not set" if value is None else str(value).lower() if isinstance(value, bool) else str(value)


def render_template(spec: dict = None) -> str:
    """The configuration template: every parameter, with its template value, and its description as a comment."""
    spec = spec or load_parameters()
    rule = "# " + "-" * (WIDTH - 2)
    out  = [f"# {spec['title']} (microbetag --write-template config.yml)", "#", _comment(spec["intro"]), ""]
    for section in spec["sections"]:
        out += [rule, f"# {section['title']}", rule, ""]
        for name, p in section["parameters"].items():
            desc = p["description"] + (" Required." if p.get("required") else "")
            if p.get("example") is not None:
                desc += f" E.g. {p['example']}"
            out += [_comment(desc), yaml.safe_dump({name: {field(p): template_value(p)}}, sort_keys=False).rstrip(), ""]
    return "\n".join(out)


def describe(name: str = None) -> str:
    """The description of a parameter, or of all of them."""
    params = all_parameters()
    if name is not None and name not in params:
        raise KeyError(f"Unknown parameter: {name}. Known: {', '.join(params)}")
    lines = []
    for n, p in params.items():
        if name is None or n == name:
            info = f"{field(p)}, default: {_show(p[field(p)])}" + (", required" if p.get("required") else "")
            lines.append(f"{n}  ({info})")
            lines.append(textwrap.indent(textwrap.fill(" ".join(p["description"].split()), WIDTH - 4), "    "))
    return "\n".join(lines)


def unknown_parameters(conf: dict) -> list:
    """Parameters of a configuration that microbetag does not know (e.g. a typo)."""
    known = set(all_parameters())
    return sorted(k for k in conf if not k.startswith("__") and k not in known)


def render_markdown(spec: dict = None) -> str:
    """The parameters as Markdown tables, one per section of the template, for the docs."""
    spec = spec or load_parameters()

    def cell(text):
        return " ".join(str(text).split()).replace("|", "\\|")

    out = []
    for section in spec["sections"]:
        out += [f"### {section['title']}", "", "| Parameter | Key | Default | In the template | Description |",
                "|---|---|---|---|---|"]
        for name, p in section["parameters"].items():
            desc = cell(p["description"]) + (" **Required.**" if p.get("required") else "")
            if p.get("example") is not None:
                desc += f" E.g. `{p['example']}`."
            shown = template_value(p)
            out.append(f"| `{name}` | `{field(p)}` | `{_show(p[field(p)])}` | "
                       f"{'`' + _show(shown) + '`' if shown is not None else ''} | {desc} |")
        out.append("")
    out += ["### Set by microbetag", "",
            "Not in the template: the on-the-fly server and the pre-calculations set them.", "",
            "| Parameter | Key | Default | Description |", "|---|---|---|---|"]
    for name, p in spec.get("internal", {}).items():
        out.append(f"| `{name}` | `{field(p)}` | `{_show(p[field(p)])}` | {cell(p['description'])} |")
    return "\n".join(out) + "\n"
