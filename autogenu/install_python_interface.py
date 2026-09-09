import argparse
import importlib
import shutil
import sys
import sysconfig
from os import PathLike
from pathlib import Path
from typing import Optional, Union

Pathish = Union[str, PathLike]


def _active_site_packages() -> Path:
    """Return the platform-specific site-packages of the running Python."""
    return Path(sysconfig.get_path("platlib")).resolve()


def install_python_interface(
    project_root_dir: Pathish,
    ocp_name: str,
    install_prefix: Optional[Pathish] = None,
) -> Path:
    """Install generated bindings into the active Python environment.

    If ``install_prefix`` is omitted, the bindings are installed into the
    platform-specific site-packages directory of the running interpreter. In
    an activated virtual environment, this is the virtual environment's own
    site-packages directory.
    """
    if install_prefix is None:
        install_prefix = _active_site_packages()
    else:
        install_prefix = Path(install_prefix).expanduser().resolve()

    project_root_dir = Path(project_root_dir).resolve()
    install_destination = install_prefix / "cgmres"
    build_dir = project_root_dir / "build" / "python"

    def collect_files(directory, patterns):
        return [
            path
            for pattern in patterns
            for path in directory.rglob(pattern)
            if path.is_file()
        ]

    binding_patterns = ("*.so", "*.dylib", "*.pyd")
    ocp_bindings = collect_files(build_dir / ocp_name, binding_patterns)
    common_bindings = collect_files(build_dir / "common", binding_patterns)
    if not ocp_bindings or not common_bindings:
        raise FileNotFoundError(
            "Generated Python bindings were not found. "
            "Run build_python_interface() before installing them."
        )

    ocp_python_files = collect_files(project_root_dir / "python" / ocp_name, ("*.py",))
    common_python_files = collect_files(project_root_dir / "python" / "common", ("*.py",))

    install_destination.mkdir(parents=True, exist_ok=True)
    (install_destination / "__init__.py").touch(exist_ok=True)
    for module_name, bindings, python_files in (
        (ocp_name, ocp_bindings, ocp_python_files),
        ("common", common_bindings, common_python_files),
    ):
        module_destination = install_destination / module_name
        module_destination.mkdir(parents=True, exist_ok=True)
        for source in (*bindings, *python_files):
            shutil.copy2(source, module_destination / source.name)

    importlib.invalidate_caches()
    print(f"Python interfaces have been installed at {install_destination}")
    print(f"Interpreter: {sys.executable}")
    return install_destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Install generated cgmres Python bindings")
    parser.add_argument("project_root_dir")
    parser.add_argument("ocp_name")
    parser.add_argument("install_prefix", nargs="?", default=None)
    args = parser.parse_args()
    install_python_interface(args.project_root_dir, args.ocp_name, args.install_prefix)
