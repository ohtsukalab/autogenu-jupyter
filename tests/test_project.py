from pathlib import Path

import nbformat

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_notebooks_are_valid():
    for notebook_path in REPOSITORY_ROOT.glob("*.ipynb"):
        notebook = nbformat.read(notebook_path, as_version=4)
        nbformat.validate(notebook)


def test_collection_distribution_is_not_a_dependency():
    requirements = (REPOSITORY_ROOT / "requirements.txt").read_text(encoding="utf-8")
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "\ncollection\n" not in f"\n{requirements}\n"
    assert '"collection"' not in pyproject
