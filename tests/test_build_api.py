from pathlib import Path
from unittest.mock import call

import autogenu.autogenu as build_api


def test_generator_aliases_and_native_auto():
    assert build_api._cmake_generator_args("Auto") == []
    assert build_api._cmake_generator_args("MSYS") == ["-G", "MSYS Makefiles"]
    assert build_api._cmake_generator_args("MinGW") == ["-G", "MinGW Makefiles"]
    assert build_api._cmake_generator_args("Ninja") == ["-G", "Ninja"]


def test_build_cpp_uses_source_build_and_config(monkeypatch, tmp_path):
    calls = []

    def fake_run(command, check):
        calls.append(call(command, check=check))

    monkeypatch.setattr(build_api.subprocess, "run", fake_run)
    build_dir = tmp_path / "project" / "build"
    build_api.build_cpp(
        "Ninja", build_dir, ["-DFEATURE=ON"], config="Debug", parallel=2
    )
    assert calls == [
        call(
            [
                "cmake", "-S", str(build_dir.parent.resolve()),
                "-B", str(build_dir.resolve()), "-G", "Ninja", "-DFEATURE=ON",
            ],
            check=True,
        ),
        call(
            [
                "cmake", "--build", str(build_dir.resolve()),
                "--config", "Debug", "--parallel", "2",
            ],
            check=True,
        ),
    ]


def test_find_executable_supports_multi_config(monkeypatch, tmp_path):
    monkeypatch.setattr(build_api.platform, "system", lambda: "Windows")
    executable = tmp_path / "Release" / "sample.exe"
    executable.parent.mkdir()
    executable.touch()
    assert build_api.find_executable(tmp_path, "sample") == executable
