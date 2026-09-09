from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_HEADERS = REPOSITORY_ROOT / "include" / "cgmres" / "python"


def test_bindings_use_functions_instead_of_multiline_macros():
    headers = {
        header.name: header.read_text()
        for header in sorted(PYTHON_HEADERS.glob("*.hpp"))
    }
    all_bindings = "\n".join(headers.values())

    assert "#define DEFINE_PYBIND11" not in all_bindings
    assert "\\\n" not in all_bindings

    expected_functions = {
        "horizon.hpp": "bind_horizon",
        "multiple_shooting_cgmres_solver.hpp": "bind_multiple_shooting_cgmres_solver",
        "single_shooting_cgmres_solver.hpp": "bind_single_shooting_cgmres_solver",
        "solver_settings.hpp": "bind_solver_settings",
        "timer.hpp": "bind_timer",
        "zero_horizon_ocp_solver.hpp": "bind_zero_horizon_ocp_solver",
    }
    for header, function in expected_functions.items():
        assert function in headers[header]
        assert headers[header].startswith("#pragma once\n")
        assert "namespace cgmres" in headers[header]
