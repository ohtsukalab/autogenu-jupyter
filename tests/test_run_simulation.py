from pathlib import Path
from unittest.mock import call

from autogenu import AutoGenU
from autogenu import autogenu as autogenu_module


def test_run_simulation_uses_build_root_for_multi_config_executable(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    generator = AutoGenU("sample", nx=1, nu=1)
    build_dir = Path(generator.get_ocp_build_dir()).resolve()
    executable = build_dir / "Release" / "sample.exe"
    executable.parent.mkdir(parents=True)
    executable.touch()
    monkeypatch.setattr(generator, "get_executable_path", lambda: executable)

    calls = []
    monkeypatch.setattr(
        autogenu_module.subprocess,
        "run",
        lambda command, cwd, check: calls.append(
            call(command, cwd=cwd, check=check)
        ),
    )

    generator.run_simulation()

    assert calls == [call([str(executable)], cwd=build_dir, check=True)]
    assert Path(generator.get_ocp_log_dir()).is_dir()
