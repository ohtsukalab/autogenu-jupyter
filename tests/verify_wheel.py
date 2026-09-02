"""Validate the built wheel without importing from the source checkout."""

import os
import subprocess
import sys
import tempfile
import venv
import zipfile
from pathlib import Path

EXPECTED_MODULES = {
    "autogenu/__init__.py",
    "autogenu/animator.py",
    "autogenu/autogenu.py",
    "autogenu/install_python_interface.py",
    "autogenu/integrator.py",
    "autogenu/logger.py",
    "autogenu/plotter.py",
    "autogenu/symutils.py",
}


def python_in(environment):
    if os.name == "nt":
        return environment / "Scripts" / "python.exe"
    return environment / "bin" / "python"


def main():
    distribution_dir = Path(sys.argv[1]).resolve()
    wheels = sorted(distribution_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"Expected exactly one wheel in {distribution_dir}, found {wheels}")
    wheel = wheels[0]

    with zipfile.ZipFile(wheel) as archive:
        members = set(archive.namelist())
    missing = EXPECTED_MODULES - members
    if missing:
        raise RuntimeError(f"Wheel is missing package files: {sorted(missing)}")
    if any("__pycache__" in member or member.endswith(".pyc") for member in members):
        raise RuntimeError("Wheel contains Python cache files")

    metadata_command = (
        "import importlib.metadata as metadata; "
        "requirements = metadata.requires('autogenu-jupyter'); "
        "assert any(r.startswith('ipykernel') and 'vscode' in r for r in requirements); "
        "assert any(r.startswith('jupyterlab') and 'jupyter' in r for r in requirements); "
        "assert any(r.startswith('matplotlib') and 'plot' in r for r in requirements)"
    )

    with tempfile.TemporaryDirectory(prefix="autogenu-wheel-") as temporary:
        temporary = Path(temporary)
        environment = temporary / "venv"
        venv.EnvBuilder(with_pip=True).create(environment)
        python = python_in(environment)
        subprocess.run(
            [str(python), "-m", "pip", "install", str(wheel)],
            check=True,
        )
        subprocess.run([str(python), "-c", metadata_command], check=True)
        clean_environment = os.environ.copy()
        clean_environment.pop("PYTHONPATH", None)

        vscode_environment = temporary / "vscode-venv"
        venv.EnvBuilder(with_pip=True).create(vscode_environment)
        vscode_python = python_in(vscode_environment)
        subprocess.run(
            [str(vscode_python), "-m", "pip", "install", f"{wheel}[vscode]"],
            check=True,
        )
        subprocess.run(
            [
                str(vscode_python),
                "-c",
                (
                    "import ipykernel, matplotlib, seaborn; "
                    "from autogenu import Plotter; "
                    "assert Plotter"
                ),
            ],
            cwd=temporary,
            env=clean_environment,
            check=True,
        )
        subprocess.run(
            [
                str(python),
                "-c",
                (
                    "import importlib.metadata, pathlib, sys, autogenu; "
                    "path = pathlib.Path(autogenu.__file__).resolve(); "
                    "assert 'site-packages' in str(path), path; "
                    "assert 'matplotlib' not in sys.modules; "
                    "requirements = importlib.metadata.requires('autogenu-jupyter'); "
                    "core = [r for r in requirements if 'extra ==' not in r]; "
                    "assert all(not r.startswith(('jupyter', 'matplotlib', "
                    "'notebook', 'seaborn')) for r in core), core; "
                    "print(path)"
                ),
            ],
            cwd=temporary,
            env=clean_environment,
            check=True,
        )


if __name__ == "__main__":
    main()
