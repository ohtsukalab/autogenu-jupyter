"""Public AutoGenU API with plotting components loaded on demand."""

from importlib import import_module

from .autogenu import *
from .install_python_interface import *
from .integrator import *
from .logger import *

_OPTIONAL_EXPORTS = {
    "Plotter": (".plotter", "plot"),
    "TwoLinkArm": (".animator", "plot"),
    "CartPole": (".animator", "plot"),
    "Hexacopter": (".animator", "plot"),
    "MobileRobot": (".animator", "plot"),
}


def __getattr__(name):
    """Load plotting helpers only when they are requested."""
    optional_export = _OPTIONAL_EXPORTS.get(name)
    if optional_export is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, extra = optional_export
    try:
        value = getattr(import_module(module_name, __name__), name)
    except ImportError as error:
        raise ImportError(
            f"{name} requires optional plotting dependencies. "
            f'Install them with: python -m pip install ".[{extra}]"'
        ) from error
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_OPTIONAL_EXPORTS))
