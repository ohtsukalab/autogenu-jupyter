import hashlib
import json
import os
from pathlib import Path

import autogenu

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = REPOSITORY_ROOT / "tests" / "snapshots" / "minimal_generation.json"


def _generate_minimal_problem(work_directory):
    previous_directory = Path.cwd()
    os.chdir(work_directory)
    try:
        generator = autogenu.AutoGenU("snapshot_minimal", nx=1, nu=1)
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
        return work_directory / "generated" / "snapshot_minimal"
    finally:
        os.chdir(previous_directory)


def _snapshot_manifest(generated_directory):
    manifest = {}
    for generated_file in sorted(generated_directory.rglob("*")):
        if generated_file.is_file():
            relative_path = generated_file.relative_to(generated_directory).as_posix()
            manifest[relative_path] = hashlib.sha256(generated_file.read_bytes()).hexdigest()
    return manifest


def test_minimal_generation_matches_snapshot(tmp_path):
    generated_directory = _generate_minimal_problem(tmp_path)
    actual = _snapshot_manifest(generated_directory)

    if os.environ.get("UPDATE_SNAPSHOTS") == "1":
        SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
        SNAPSHOT.write_text(json.dumps(actual, indent=2) + "\n", encoding="utf-8")

    expected = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    assert actual == expected
