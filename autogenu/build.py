"""Cross-platform CMake build primitives used by :class:`AutoGenU`."""

import platform
import shutil
import subprocess
from os import PathLike
from pathlib import Path
from typing import List, Optional, Sequence, Union

Pathish = Union[str, PathLike[str]]


def cmake_generator_args(generator: str) -> List[str]:
    """Translate legacy aliases while allowing every CMake generator."""
    aliases = {"MSYS": "MSYS Makefiles", "MinGW": "MinGW Makefiles"}
    if not generator or generator == "Auto":
        return []
    return ["-G", aliases.get(generator, generator)]


def build_cpp(
    generator: str,
    build_dir: Pathish,
    build_options: Sequence[str],
    config: str = "Release",
    parallel: Optional[int] = None,
) -> Path:
    """Configure and build a CMake project, raising immediately on failure."""
    build_dir = Path(build_dir).resolve()
    source_dir = build_dir.parent
    build_dir.mkdir(parents=True, exist_ok=True)
    configure_command = [
        "cmake",
        "-S",
        str(source_dir),
        "-B",
        str(build_dir),
        *cmake_generator_args(generator),
        *build_options,
    ]
    build_command = ["cmake", "--build", str(build_dir), "--config", config]
    if parallel is not None:
        build_command.extend(["--parallel", str(parallel)])
    print("Configure command:", *configure_command)
    subprocess.run(configure_command, check=True)
    print("Build command:", *build_command)
    subprocess.run(build_command, check=True)
    return build_dir


def find_executable(
    build_dir: Pathish, target_name: str, config: str = "Release"
) -> Path:
    """Locate an executable from single- or multi-configuration generators."""
    build_dir = Path(build_dir).resolve()
    executable_name = target_name + (".exe" if platform.system() == "Windows" else "")
    candidates = [build_dir / executable_name, build_dir / config / executable_name]
    candidates.extend(build_dir.glob("*/" + executable_name))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Executable '{executable_name}' was not found below '{build_dir}'.")


def remove_build_directory(
    project_dir: Pathish, directory_name: str = "build"
) -> None:
    """Remove a generated build directory without invoking a platform shell."""
    shutil.rmtree(Path(project_dir) / directory_name, ignore_errors=True)
