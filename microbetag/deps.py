"""
Checks for optional dependencies, with a how-to-install hint in the error message.

The external tools and conda environments are installed by install.sh (see envs/*.yml),
so every hint points to the install.sh option that provides the missing piece.
"""
import shutil
import subprocess
from functools import lru_cache

INSTALL_DOCS = "https://microbetag.readthedocs.io/en/latest/installation.html"

# conda environment -> (feature that needs it, install.sh option)
CONDA_ENVS = {
    "mtg-phenotrex": ("phenotypic trait prediction (phenotrex)", "--phenotrex"),
    "mtg-modelseed": ("GEM reconstruction with ModelSEEDpy", "--modelseed"),
    "mtg-dnngior"  : ("gap-filling with DNNGIOR", "--dnngior"),
}

# Tools that ship with the `microbetag` conda environment (envs/microbetag.yml)
TOOLS = {
    "prodigal"    : "ORF prediction",
    "hmmsearch"   : "KEGG ORTHOLOGY annotation",
    "diamond"     : "GEM reconstruction with CarveMe",
    "carve"       : "GEM reconstruction with CarveMe",
    "FragGeneScan": "ORF prediction for ModelSEEDpy reconstructions",
    "julia"       : "network inference with FlashWeave",
    "manta"       : "network clustering with manta",
}


class MissingDependencyError(RuntimeError):
    """Raised when an optional dependency needed for the requested step is not installed."""


def _message(what: str, feature: str, how: str) -> str:
    return (
        f"{what} is not available, but it is needed for {feature}.\n"
        f"To install it, run from the microbetag repository:\n\n"
        f"    {how}\n\n"
        f"See {INSTALL_DOCS} for details."
    )


@lru_cache(maxsize=None)
def _conda_envs() -> frozenset:
    """Names of the existing conda environments (empty if conda is not available)."""
    if shutil.which("conda") is None:
        return frozenset()
    res = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    return frozenset(
        line.split()[0] for line in res.stdout.splitlines()
        if line.strip() and not line.startswith("#")
    )


def require_conda_env(name: str) -> None:
    """Raise MissingDependencyError if the conda environment `name` does not exist."""
    if name in _conda_envs():
        return
    feature, option = CONDA_ENVS[name]
    raise MissingDependencyError(
        _message(f"The '{name}' conda environment", feature, f"bash install.sh {option}")
    )


def missing_tool_error(tool: str) -> MissingDependencyError:
    """Error for a missing tool of the `microbetag` conda environment."""
    feature = TOOLS.get(tool, "this step")
    return MissingDependencyError(
        _message(
            f"'{tool}'",
            feature,
            "bash install.sh\n\n"
            "and run microbetag from that environment:\n\n"
            "    conda activate microbetag",
        )
    )


def init_flashweave():
    """Start Julia and load FlashWeave, with an install hint if either is missing."""
    how = "bash install.sh            # FlashWeave is installed unless --no-flashweave is given"
    feature = "network inference with FlashWeave"
    try:
        from julia.api import Julia
    except ImportError:
        raise MissingDependencyError(
            _message("The 'julia' Python package (pyjulia)", feature, how)
        ) from None

    if shutil.which("julia") is None:
        raise missing_tool_error("julia")

    try:
        jl = Julia(compiled_modules=False)
        jl.using("FlashWeave")
    except Exception as e:
        raise MissingDependencyError(
            _message("The FlashWeave Julia package", feature, how) + f"\n\nOriginal error: {e}"
        ) from e
    return jl
