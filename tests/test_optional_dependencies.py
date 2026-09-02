import sys

import autogenu


def test_core_import_does_not_eagerly_import_plotting_dependencies():
    assert "matplotlib" not in sys.modules
    assert "seaborn" not in sys.modules
    assert autogenu.AutoGenU


def test_plotting_api_is_available_lazily():
    assert autogenu.Plotter
    assert autogenu.CartPole
