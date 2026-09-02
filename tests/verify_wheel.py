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

    with tempfile.TemporaryDirectory(prefix="autogenu-wheel-") as temporary:
        temporary = Path(temporary)
        environment = temporary / "venv"
        venv.EnvBuilder(with_pip=True).create(environment)
        python = python_in(environment)
        subprocess.run(
            [str(python), "-m", "pip", "install", str(wheel)],
            check=True,
        )
        clean_environment = os.environ.copy()
        clean_environment.pop("PYTHONPATH", None)
        subprocess.run(
            [
                str(python),
                "-c",
                (
                    "import pathlib, autogenu; "
                    "path = pathlib.Path(autogenu.__file__).resolve(); "
                    "assert 'site-packages' in str(path), path; "
                    "print(path)"
                ),
            ],
            cwd=temporary,
            env=clean_environment,
            check=True,
        )


if __name__ == "__main__":
    main()
