from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_eigen_stack_limit_is_configured_before_eigen_include():
    types_header = (REPOSITORY_ROOT / "include" / "cgmres" / "types.hpp").read_text()
    definition = types_header.index("#define EIGEN_STACK_ALLOCATION_LIMIT 0")
    eigen_include = types_header.index('#include "cgmres/thirdparty/eigen/Eigen/Core"')

    assert "#ifndef EIGEN_STACK_ALLOCATION_LIMIT" in types_header
    assert definition < eigen_include


def test_generated_cmake_has_cross_compiler_warning_policy():
    cmake_template = (
        REPOSITORY_ROOT / "autogenu" / "templates" / "CMakeLists.txt.in"
    ).read_text()

    for expected in (
        "CGMRES_WARNINGS_AS_ERRORS",
        "/W4",
        "/WX",
        "/wd4100",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Werror",
        "-Wno-unused-parameter",
    ):
        assert expected in cmake_template


def test_generated_cmake_has_gcc_and_clang_sanitizer_policy():
    cmake_template = (
        REPOSITORY_ROOT / "autogenu" / "templates" / "CMakeLists.txt.in"
    ).read_text()

    for expected in (
        "CGMRES_ENABLE_SANITIZERS",
        "-fsanitize=address,undefined",
        "-fno-omit-frame-pointer",
    ):
        assert expected in cmake_template
