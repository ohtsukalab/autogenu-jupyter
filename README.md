# AutoGenU for Jupyter

[![build](https://github.com/ohtsukalab/autogenu-jupyter/actions/workflows/build.yaml/badge.svg?branch=master)](https://github.com/ohtsukalab/autogenu-jupyter/actions/workflows/build.yaml)
[![doxygen](https://github.com/ohtsukalab/autogenu-jupyter/actions/workflows/doxygen.yaml/badge.svg)](https://github.com/ohtsukalab/autogenu-jupyter/actions/workflows/doxygen.yaml)


## Introduction
This project provides the continuation/GMRES method (C/GMRES method) based solvers for nonlinear model predictive control (NMPC) and an automatic code generator for NMPC, called `AutoGenU`.

The following C/GMRES based solvers are provided: 
- `MultipleShootingCGMRESSolver` : The multiple shooting based C/GMRES method with condensing of the state and costate directions.
- `SingleShootingCGMRESSolver` : The original C/GMRES method (single shooting).

## Requirement
- C++17 compiler (GCC, Clang, or MSVC)
- CMake 4, git
- Python 3.8 or later, SymPy, and NumPy for the core code-generation API
- Jupyter, VS Code kernel, and plotting packages are available as optional extras
- ffmpeg (to generate animations in the example notebooks)
- Doxygen (optional, to generate C++ docs)


## Usage
### 1. Setup requirements
Please confirm that you clone this repository as 
```
git clone https://github.com/ohtsukalab/autogenu-jupyter --recursive
```
Otherwise, please do the following command:
```
git submodule update --init --recursive
```
Create and activate a virtual environment, then install the Python package via
```
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
```
The default installation is intentionally minimal and installs only NumPy and
SymPy. Choose an extra for the environment you use:
```
# VS Code notebooks: kernel support plus plotting
python -m pip install ".[vscode]"

# JupyterLab or Jupyter Notebook plus plotting
python -m pip install ".[jupyter]"

# Plotting helpers without a notebook frontend
python -m pip install ".[plot]"

# Contributor environment (tests, packaging tools, and notebooks)
python -m pip install ".[dev]"
```

In VS Code connected to WSL, select
`.venv/bin/python` with **Notebook: Select Notebook Kernel**. Confirm the
selected kernel from a notebook cell with:
```python
import sys
print(sys.executable)
```

### 2. Code generation
`AutoGenU.ipynb` generates the following source files under your setting state equation, constraints, cost function, and parameters: 
- `ocp.hpp` : A definition of the optimal control problem (OCP).
- `main.cpp` : An executablb of the closed-loop simulation.
- `CMakeLists.txt` : Scripts to build C++ projects. 
- Files in `python` directory : Source files of Python interface via pybind11.

You can generate these files, run simulations, plot results, and install the Python interfaces through `AutoGenU.ipynb`.

The build API uses CMake consistently on Linux, macOS, and Windows:
```python
# Let CMake select the native generator. On Windows this normally uses MSVC.
ag.build_main(generator="Auto", config="Release", parallel=2)

# Explicit generators such as Ninja are also supported.
ag.build_python_interface(generator="Ninja", config="Release")
```
The legacy `MSYS` and `MinGW` generator names remain available. Build failures
raise `subprocess.CalledProcessError`, and `ag.get_executable_path()` locates
executables produced by both single- and multi-configuration generators.

### Public Python API

The supported top-level API is explicitly defined by `autogenu.__all__` and
contains only problem-independent functionality: `AutoGenU`, `NLPType`,
integration and logging helpers, documentation helpers, and the generic
`Plotter`. Internal CMake helpers and example-specific animators are not
exported at the package top level.

Example-specific animation helpers remain available from their module when
needed by the bundled examples:
```python
from autogenu.animator import CartPole, Hexacopter, MobileRobot, TwoLinkArm
```

Advanced users can access the low-level, cross-platform build primitives from
the dedicated module:
```python
from autogenu.build import build_cpp, cmake_generator_args, find_executable
```
Application code should normally use `AutoGenU.build_main()` and
`AutoGenU.build_python_interface()` instead.


### 3. Python bindings
Python bindings are built and installed via `.ipynb` files. Activate the
virtual environment before starting Jupyter; the bindings are installed into
that environment's `site-packages` directory by default:
```
source .venv/bin/activate
python -m pip install ".[jupyter]"
jupyter lab
```
No manual `PYTHONPATH` setting is required. The interfaces can be imported as
```
import cgmres.common # this includes horizon, solver settings, etc.
import cgmres.your_ocp_name # this includes OCP definition and NMPC solvers 
```


### 4. Install header-only `cgmres` C++ library
Aside from the notebook for the code-generation, the C++ `cgmres` library, which is a header-only library, can be installed by running
```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=YOUR_INSTALL_DESTINATION
make install 
```
at the project root directory of `autogenu-jupyter`.  
Then you can build the NMPC code with the generated `ocp.hpp` file and without `.ipynb` notebook files.   
The examples are found in `examples/cpp` directory.  


### 5. Install `autogenu` Python module
The Python module `autogenu` can be installed by running
```
python3 -m pip install .
```
at the project root directory of `autogenu-jupyter`.
Further, if you install have installed header-only `cgmres` C++ library as step 4, then you can run `.ipynb` files for the code generation in everywhere.


### Documentation
C++ API documentation of `cgmres` library is found at https://ohtsukalab.github.io/autogenu-jupyter/.   
Python interfaces are almost the same as the C++ API, so please refere to https://ohtsukalab.github.io/autogenu-jupyter/ even for Python interfaces as well as the [tips for conversions between C++ and Python](https://ohtsukalab.github.io/autogenu-jupyter/md__github_workspace_examples_python__r_e_a_d_m_e.html).


## Demos
Demos are presented in `cartpole.ipynb`, `pendubot.ipynb`, `hexacopter.ipynb`, and `mobilerobot.ipynb`. You can obtain the following simulation results jusy by runnig these `.ipynb` files. The details of the each OCP formulations are described in each `.ipynb` files.

<img src="https://raw.githubusercontent.com/wiki/mayataka/CGMRES/images/cartpole.gif" width="300"> &nbsp;
<img src="https://raw.githubusercontent.com/wiki/mayataka/CGMRES/images/pendubot.gif" width="300"> 

<img src="https://raw.githubusercontent.com/wiki/mayataka/CGMRES/images/hexacopter.gif" width="450">  

<img src="https://raw.githubusercontent.com/wiki/mayataka/CGMRES/images/mobilerobot.gif" width="450"> 


## License
MIT

## Citing autogenu-jupyter

We'd appriciate if you use cite the following conference paper:

```
@inproceedings{katayama2020autogenu,
  title={Automatic code generation tool for nonlinear model predictive control with {J}upyter},
  author={Sotaro Katayama and Toshiyuki Ohtsuka},
  booktitle={{The 21st IFAC World Congress 2020}},
  pages={7033-7040},
  year={2020}}
```

## References
1. [T. Ohtsuka A continuation/GMRES method for fast computation of nonlinear receding horizon control, Automatica, Vol. 40, No. 4, pp. 563-574 (2004)](https://doi.org/10.1016/j.automatica.2003.11.005)
2. [C. T. Kelly, Iterative methods for linear and nonlinear equations, Frontiers in Apllied Mathematics, SIAM (1995)](https://doi.org/10.1137/1.9781611970944)
3. [Y. Shimizu, T. Ohtsuka, M. Diehl, A real‐time algorithm for nonlinear receding horizon control using multiple shooting and continuation/Krylov method, International Journal of Robust and Nonlinear Control, Vol. 19, No. 8, pp. 919-936 (2008)](https://doi.org/10.1002/rnc.1363)
