import os
import shutil
import subprocess
from pathlib import Path

import pytest

import autogenu

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
OCP_NAME = "sanitizer_minimal"

pytestmark = pytest.mark.skipif(
    os.environ.get("CGMRES_RUN_SANITIZERS") != "1",
    reason="set CGMRES_RUN_SANITIZERS=1 to run sanitizer tests",
)


@pytest.mark.e2e
def test_generated_simulation_with_asan_and_ubsan(monkeypatch):
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
        generator.generate_cmake()
        generator.build_main(
            vectorize=False,
            remove_build_dir=True,
            warnings_as_errors=True,
            sanitizers=True,
        )

        result = subprocess.run(
            [str(generator.get_executable_path())],
            cwd=generated_dir / "build",
            check=True,
            capture_output=True,
            text=True,
        )
        assert "error" not in result.stderr.lower()
    finally:
        shutil.rmtree(generated_dir, ignore_errors=True)
