"""Public AutoGenU API with plotting components loaded on demand."""

from importlib import import_module

from .autogenu import AutoGenU, NLPType, generate_docs, open_docs
from .install_python_interface import install_python_interface
from .integrator import RK4, forward_euler
from .logger import Logger

__all__ = [
    "AutoGenU",
    "NLPType",
    "Logger",
    "forward_euler",
    "RK4",
    "install_python_interface",
    "generate_docs",
    "open_docs",
    "Plotter",
]

_OPTIONAL_EXPORTS = {
    "Plotter": (".plotter", "plot"),
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
