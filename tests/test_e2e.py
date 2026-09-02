import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

import autogenu


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
OCP_NAME = "e2e_minimal"


def _python_in(venv_dir):
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


@pytest.mark.e2e
def test_generate_build_run_install_and_import(tmp_path, monkeypatch):
    """Exercise the complete generated-code workflow with a minimal OCP."""
    generated_dir = REPOSITORY_ROOT / "generated" / OCP_NAME
    shutil.rmtree(generated_dir, ignore_errors=True)
    monkeypatch.chdir(REPOSITORY_ROOT)

    try:
        generator = autogenu.AutoGenU(OCP_NAME, nx=1, nu=1)
        x = generator.define_x()
        u = generator.define_u()
        generator.set_functions(
            f=[u[0]],
            C=[],
            h=[],
            L=(x[0] ** 2 + u[0] ** 2) / 2,
            phi=x[0] ** 2 / 2,
        )
        generator.generate_ocp_definition()
        generator.set_nlp_type(autogenu.NLPType.SingleShooting)
        generator.set_horizon_params(Tf=0.1)
        generator.set_solver_params(0.01, 5, 1.0e-6, 100.0, 3)
        generator.set_initialization_params([0.0], 1.0e-8, 20)
        generator.set_simulation_params(0.0, [0.1], 0.02)
        generator.generate_main()
        generator.generate_python_bindings()
        generator.generate_cmake()

        generator.build_main(vectorize=False, remove_build_dir=True)
        executable = generator.get_executable_path()
        result = subprocess.run(
            [str(executable)],
            cwd=generated_dir / "build",
            check=True,
            capture_output=True,
            text=True,
        )
        assert "error" not in result.stderr.lower()

        generator.build_python_interface(vectorize=False, remove_build_dir=True)

        venv_dir = tmp_path / "consumer-venv"
        subprocess.run(
            [sys.executable, "-m", "venv", "--system-site-packages", str(venv_dir)],
            check=True,
        )
        consumer_python = _python_in(venv_dir)
        installer = REPOSITORY_ROOT / "autogenu" / "install_python_interface.py"
        subprocess.run(
            [str(consumer_python), str(installer), str(generated_dir), OCP_NAME],
            check=True,
        )
        subprocess.run(
            [
                str(consumer_python),
                "-c",
                (
                    "import cgmres.common; "
                    f"import cgmres.{OCP_NAME}; "
                    "print('generated bindings: OK')"
                ),
            ],
            cwd=tmp_path,
            check=True,
        )
    finally:
        shutil.rmtree(generated_dir, ignore_errors=True)
