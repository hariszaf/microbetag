"""Tests for the install hints of optional dependencies (microbetag.deps)."""

import pytest

from microbetag import deps


def test_missing_conda_env_points_to_install_option(monkeypatch):
    monkeypatch.setattr(deps, "_conda_envs", lambda: frozenset({"base", "microbetag"}))
    with pytest.raises(deps.MissingDependencyError, match=r"bash install.sh --phenotrex"):
        deps.require_conda_env("mtg-phenotrex")


def test_existing_conda_env_passes(monkeypatch):
    monkeypatch.setattr(deps, "_conda_envs", lambda: frozenset({"mtg-dnngior"}))
    deps.require_conda_env("mtg-dnngior")


def test_missing_tool_names_feature():
    err = deps.missing_tool_error("prodigal")
    assert "ORF prediction" in str(err)
    assert "conda activate microbetag" in str(err)


def test_flashweave_without_julia(monkeypatch):
    monkeypatch.setattr(deps.shutil, "which", lambda _: None)
    with pytest.raises(deps.MissingDependencyError, match="FlashWeave"):
        deps.init_flashweave()
