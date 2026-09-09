from pathlib import Path

import nbformat

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_notebooks_are_valid():
    for notebook_path in REPOSITORY_ROOT.glob("*.ipynb"):
        notebook = nbformat.read(notebook_path, as_version=4)
        nbformat.validate(notebook)


def test_notebook_display_math_uses_standalone_delimiters():
    for notebook_path in REPOSITORY_ROOT.glob("*.ipynb"):
        notebook = nbformat.read(notebook_path, as_version=4)
        for cell_index, cell in enumerate(notebook.cells):
            if cell.cell_type != "markdown":
                continue
            for line_number, line in enumerate(cell.source.splitlines(), start=1):
                assert "$$" not in line or line.strip() == "$$", (
                    f"{notebook_path.name}: markdown cell {cell_index}, "
                    f"line {line_number} must place each $$ delimiter on its own line"
                )


def test_hexacopter_installs_generated_bindings_from_a_code_cell():
    notebook = nbformat.read(REPOSITORY_ROOT / "hexacopter.ipynb", as_version=4)
    install_cells = [
        cell
        for cell in notebook.cells
        if "ag.install_python_interface(" in cell.source
    ]
    assert len(install_cells) == 1
    assert install_cells[0].cell_type == "code"


def test_collection_distribution_is_not_a_dependency():
    requirements = (REPOSITORY_ROOT / "requirements.txt").read_text(encoding="utf-8")
    pyproject = (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "\ncollection\n" not in f"\n{requirements}\n"
    assert '"collection"' not in pyproject
