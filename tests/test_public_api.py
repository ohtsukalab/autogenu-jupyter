import autogenu

EXPECTED_PUBLIC_API = {
    "AutoGenU",
    "NLPType",
    "Logger",
    "forward_euler",
    "RK4",
    "install_python_interface",
    "generate_docs",
    "open_docs",
    "Plotter",
}


def test_public_api_is_explicit_and_stable():
    assert set(autogenu.__all__) == EXPECTED_PUBLIC_API
    assert len(autogenu.__all__) == len(EXPECTED_PUBLIC_API)


def test_internal_build_helpers_are_not_top_level_exports():
    assert "build_cpp" not in autogenu.__all__
    assert not hasattr(autogenu, "build_cpp")
    assert not hasattr(autogenu, "subprocess")
    assert not hasattr(autogenu, "sympy")


def test_example_animators_are_not_top_level_exports():
    for name in ("TwoLinkArm", "CartPole", "Hexacopter", "MobileRobot"):
        assert name not in autogenu.__all__
        assert not hasattr(autogenu, name)
