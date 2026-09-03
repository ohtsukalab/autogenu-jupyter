import math

import pytest

import autogenu


def _generator_with_functions(*, inequality=False):
    generator = autogenu.AutoGenU("validation_problem", nx=2, nu=1)
    x = generator.define_x()
    u = generator.define_u()
    h = [x[0] - 1] if inequality else []
    generator.set_functions(
        f=[x[1], u[0]],
        C=[],
        h=h,
        L=(x[0] ** 2 + x[1] ** 2 + u[0] ** 2) / 2,
        phi=(x[0] ** 2 + x[1] ** 2) / 2,
    )
    return generator


@pytest.mark.parametrize("name", ["two-link", "2link", "class", "for", ""])
def test_ocp_name_must_be_a_cross_language_identifier(name):
    with pytest.raises(ValueError, match="valid non-keyword C\\+\\+ and Python"):
        autogenu.AutoGenU(name, nx=1, nu=1)


@pytest.mark.parametrize("argument,value", [("nx", 0), ("nu", -1)])
def test_dimensions_must_be_positive_integers(argument, value):
    kwargs = {"ocp_name": "valid", "nx": 1, "nu": 1}
    kwargs[argument] = value
    with pytest.raises(ValueError, match=rf"{argument} must be at least 1"):
        autogenu.AutoGenU(**kwargs)


def test_boolean_is_not_accepted_as_an_integer_dimension():
    with pytest.raises(TypeError, match="nx must be an integer; got bool"):
        autogenu.AutoGenU("valid", nx=True, nu=1)


def test_function_dimension_error_reports_expected_and_actual_sizes():
    generator = autogenu.AutoGenU("valid", nx=2, nu=1)
    with pytest.raises(ValueError, match="f must contain 2 values; got 1"):
        generator.set_functions(f=[0], C=[], h=[], L=0, phi=0)


def test_unknown_and_duplicate_user_variables_are_rejected():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    generator.define_scalar_var("mass")

    with pytest.raises(ValueError, match="mass.*already defined"):
        generator.define_array_var("mass", 2)
    with pytest.raises(ValueError, match="unknown scalar variable 'gravity'"):
        generator.set_scalar_var("gravity", 9.81)
    with pytest.raises(ValueError, match="unknown array variable 'weights'"):
        generator.set_array_var("weights", [1.0])


def test_array_value_dimension_is_validated():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    generator.define_array_var("weights", 2)
    with pytest.raises(ValueError, match="values for 'weights' must contain 2 values; got 1"):
        generator.set_array_var("weights", [1.0])


def test_model_parameters_allow_symbolic_cpp_constants():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    generator.define_scalar_var("angle")
    generator.define_array_var("reference", 2)
    generator.set_scalar_var("angle", "M_PI")
    generator.set_array_var("reference", ["M_PI", 0.0])


@pytest.mark.parametrize(
    "setter,args,message",
    [
        ("set_horizon_params", (0.0,), "Tf must be greater than 0.0"),
        ("set_horizon_params", (math.inf,), "Tf must be finite"),
        ("set_solver_params", (0.0, 10, 1e-8, 100.0, 3), "sampling_time"),
        ("set_solver_params", (0.01, 0, 1e-8, 100.0, 3), "N must be at least 1"),
        ("set_simulation_params", (0.0, [0.0], 1.0), "initial_state must contain 2"),
    ],
)
def test_numeric_configuration_errors_name_the_invalid_parameter(
    setter, args, message
):
    generator = _generator_with_functions()
    with pytest.raises(ValueError, match=message):
        getattr(generator, setter)(*args)


def test_control_bound_index_and_interval_are_validated():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    with pytest.raises(ValueError, match="uindex must be between 0 and 0; got 1"):
        generator.add_control_input_bounds(1, -1.0, 1.0, 0.1)
    with pytest.raises(ValueError, match="umin must be less than umax"):
        generator.add_control_input_bounds(0, 1.0, 1.0, 0.1)


def test_fb_epsilon_dimension_and_values_are_validated():
    generator = _generator_with_functions(inequality=True)
    with pytest.raises(ValueError, match="FB_epsilon must contain 1 values; got 0"):
        generator.set_FB_epsilon([])
    with pytest.raises(ValueError, match=r"FB_epsilon\[0\] must be at least 0.0"):
        generator.set_FB_epsilon([-1.0])


def test_generate_main_reports_all_missing_configuration():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    with pytest.raises(RuntimeError) as error:
        generator.generate_main()

    message = str(error.value)
    for setter in (
        "set_functions()",
        "set_nlp_type()",
        "set_horizon_params()",
        "set_solver_params()",
        "set_initialization_params()",
        "set_simulation_params()",
    ):
        assert setter in message


def test_generate_ocp_requires_functions_and_fb_epsilon(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    with pytest.raises(RuntimeError, match=r"call set_functions\(\) first"):
        generator.generate_ocp_definition()

    x = generator.define_x()
    u = generator.define_u()
    generator.set_functions(f=[u[0]], C=[], h=[x[0]], L=0, phi=0)
    with pytest.raises(ValueError, match="FB_epsilon must contain 1 values; got 0"):
        generator.generate_ocp_definition()


def test_generate_ocp_reports_unset_array_parameters(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    generator.define_array_var("weights", 2)
    generator.set_functions(f=[0], C=[], h=[], L=0, phi=0)

    with pytest.raises(RuntimeError, match=r"'weights'.*set_array_var\(\)"):
        generator.generate_ocp_definition()


def test_nlp_type_requires_enum_value():
    generator = autogenu.AutoGenU("valid", nx=1, nu=1)
    with pytest.raises(TypeError, match="nlp_type must be an NLPType value"):
        generator.set_nlp_type("SingleShooting")
