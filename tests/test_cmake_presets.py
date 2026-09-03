import json
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_cmake_presets_cover_supported_quality_modes():
    presets = json.loads((REPOSITORY_ROOT / "CMakePresets.json").read_text())
    expected = {"dev", "strict", "clang-tidy", "sanitizers"}

    for section in ("configurePresets", "buildPresets", "testPresets"):
        names = {preset["name"] for preset in presets[section]}
        assert expected <= names


def test_cmake_smoke_test_is_registered():
    cmake = (REPOSITORY_ROOT / "CMakeLists.txt").read_text()
    assert "include(CTest)" in cmake
    assert "add_executable(cgmres_smoke" in cmake
    assert "add_test(NAME cgmres.smoke" in cmake
