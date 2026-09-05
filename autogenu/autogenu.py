import keyword
import math
import numbers
import os
import platform
import re
import shutil
import subprocess
import sys
from collections import namedtuple
from enum import Enum, auto
from os import PathLike
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple, Union

from sympy.core.expr import Expr
from sympy.core.symbol import Symbol, symbols
from sympy.functions.elementary.miscellaneous import sqrt

from . import build as build_api
from . import symutils
from .install_python_interface import install_python_interface
from .template_renderer import write_generated_file

SymbolicExpression = Union[Expr, int, float]
ModelParameter = Union[SymbolicExpression, str]
Pathish = Union[str, PathLike]

_CPP_KEYWORDS = {
    "alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand",
    "bitor", "bool", "break", "case", "catch", "char", "char8_t",
    "char16_t", "char32_t", "class", "compl", "concept", "const",
    "consteval", "constexpr", "constinit", "const_cast", "continue",
    "co_await", "co_return", "co_yield", "decltype", "default", "delete",
    "do", "double", "dynamic_cast", "else", "enum", "explicit", "export",
    "extern", "false", "float", "for", "friend", "goto", "if", "inline",
    "int", "long", "mutable", "namespace", "new", "noexcept", "not",
    "not_eq", "nullptr", "operator", "or", "or_eq", "private", "protected",
    "public", "register", "reinterpret_cast", "requires", "return", "short",
    "signed", "sizeof", "static", "static_assert", "static_cast", "struct",
    "switch", "template", "this", "thread_local", "throw", "true", "try",
    "typedef", "typeid", "typename", "union", "unsigned", "using",
    "virtual", "void", "volatile", "wchar_t", "while", "xor", "xor_eq",
}
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _require_identifier(name, value):
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string; got {type(value).__name__}")
    if (
        not _IDENTIFIER_PATTERN.fullmatch(value)
        or keyword.iskeyword(value)
        or value in _CPP_KEYWORDS
    ):
        raise ValueError(
            f"{name} must be a valid non-keyword C++ and Python identifier; "
            f"got {value!r}"
        )


def _require_integer(name, value, *, minimum=1):
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{name} must be an integer; got {type(value).__name__}")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}; got {value!r}")


def _require_finite_number(name, value, *, minimum=None, strict=False):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real number; got {type(value).__name__}")
    if not math.isfinite(float(value)):
        raise ValueError(f"{name} must be finite; got {value!r}")
    if minimum is not None:
        valid = value > minimum if strict else value >= minimum
        if not valid:
            relation = "greater than" if strict else "at least"
            raise ValueError(f"{name} must be {relation} {minimum}; got {value!r}")


def _require_model_parameter(name, value):
    if isinstance(value, str):
        if not value.strip():
            raise ValueError(f"{name} must not be an empty C++ expression")
        return
    if isinstance(value, Expr):
        if value.is_finite is False:
            raise ValueError(f"{name} must be finite; got {value!r}")
        return
    _require_finite_number(name, value)


def _require_sequence(name, value):
    if isinstance(value, (str, bytes)) or not hasattr(value, "__len__"):
        raise TypeError(f"{name} must be a sized sequence; got {type(value).__name__}")


def _require_boolean(name, value):
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool; got {type(value).__name__}")


def _require_length(name, value, expected):
    _require_sequence(name, value)
    actual = len(value)
    if actual != expected:
        raise ValueError(f"{name} must contain {expected} values; got {actual}")


class ScalarVariable:
    def __init__(self, symbol: Symbol, name: str, value=0.0):
        self.symbol = symbol
        self.name = name 
        self.value = value

class ArrayVariable:
    def __init__(self, symbol, name: str, size: int, values=None):
        _require_integer("size", size)
        self.symbol = symbol
        self.name = name 
        self.size = size
        self.values = [] if values is None else values

class ControlInputBound:
    def __init__(self, uindex: int, umin, umax, dummy_weight):
        _require_integer("uindex", uindex, minimum=0)
        _require_finite_number("umin", umin)
        _require_finite_number("umax", umax)
        if umin >= umax:
            raise ValueError(f"umin must be less than umax; got {umin!r} >= {umax!r}")
        _require_finite_number("dummy_weight", dummy_weight, minimum=0.0)
        self.uindex = uindex 
        self.umin = umin
        self.umax = umax
        self.dummy_weight = dummy_weight

SymbolicFunctions = namedtuple('SymbolicFunctions', ['f', 'phix', 'hx', 'hu'])

class NLPType(Enum):
    SingleShooting = auto()
    MultipleShooting = auto()

HorizonParams = namedtuple('HorizonParams', ['Tf', 'alpha'])

SolverParams = namedtuple('SolverParams', ['sampling_time', 'N', 'finite_difference_epsilon', 'zeta', 'kmax'])

InitializationParams = namedtuple('InitializationParams', ['solution_initial_guess', 'tolerance', 'max_iteraions'])

SimulationParams = namedtuple('SimulationParams', ['initial_time', 'initial_state', 'simulation_length'])


class AutoGenU(object):
    """ Automatic C++ code generator for the C/GMRES methods. 

        Args: 
            ocp_name: The name of the optimal control problem (OCP). 
                The directory 'generated/ocp_name' is made and C++ source files 
                are generated in the directory.
            nx: The dimension of the state of the system. 
            nu: The dimension of the control input of the system. 
    """
    def __init__(self, ocp_name: str, nx: int, nu: int) -> None:
        _require_identifier("ocp_name", ocp_name)
        _require_integer("nx", nx)
        _require_integer("nu", nu)
        self.__ocp_name = ocp_name
        self.__nx = nx
        self.__nu = nu
        self.__nc = 0
        self.__nh = 0
        self.__scalar_vars = []
        self.__array_vars = []
        self.__ubounds = []
        self.__FB_epsilon = []
        self.__symbolic_functions = None
        self.__nlp_type = None
        self.__horizon_params = None
        self.__solver_params = None
        self.__initialization_params = None
        self.__simulation_params = None

    def get_ocp_name(self) -> str:
        return self.__ocp_name

    def get_ocp_dir(self) -> str:
        return os.path.join(os.getcwd(), os.path.abspath('generated'), self.__ocp_name)

    def get_ocp_pybind_dir(self) -> str:
        return os.path.join(os.getcwd(), os.path.abspath('generated'), self.__ocp_name, 'python')

    def get_ocp_build_dir(self) -> str:
        return os.path.join(os.getcwd(), self.get_ocp_dir(), 'build')

    def get_ocp_log_dir(self) -> str:
        return os.path.join(os.getcwd(), self.get_ocp_dir(), 'log')

    def _require_configured(
        self, action: str, requirements: Sequence[Tuple[Any, str]]
    ) -> None:
        missing = [setter for value, setter in requirements if value is None]
        if missing:
            calls = ", ".join(f"{setter}()" for setter in missing)
            raise RuntimeError(f"{action} requires prior configuration: call {calls} first")

    def define_t(self) -> Symbol:
        """ Returns symbolic scalar variable 't'.
        """
        return Symbol('t')

    def define_x(self) -> Tuple[Symbol, ...]:
        """ Returns symbolic vector variable 'x' whose size is nx.
        """
        return symbols('x[0:%d]' %(self.__nx))

    def define_u(self) -> Tuple[Symbol, ...]:
        """ Returns symbolic vector variable 'u' whose size is nu.
        """
        return symbols('u[0:%d]' %(self.__nu))

    def define_scalar_var(self, var_name: str) -> Symbol:
        """ Returns symbolic variable whose name is var_name. 

            Args:
                var_name: Name of the scalar variable.
        """
        _require_identifier("var_name", var_name)
        defined_names = {var.name for var in self.__scalar_vars + self.__array_vars}
        if var_name in defined_names:
            raise ValueError(f"var_name must be unique; {var_name!r} is already defined")
        var_symbol = Symbol(var_name)
        self.__scalar_vars.append(ScalarVariable(var_symbol, var_name))
        return var_symbol

    def define_scalar_vars(self, *var_name_list: str) -> List[Symbol]:
        """ Returns symbolic variables whose names are given by 
            var_name_list. 

            Args:
                var_name_list: Names of the scalar variables.
        """
        var_symbols = []
        for var_name in var_name_list:
            var_symbol = self.define_scalar_var(var_name)
            var_symbols.append(var_symbol)
        return var_symbols

    def define_array_var(
        self, var_name: str, size: int
    ) -> Tuple[Symbol, ...]:
        """ Returns symbolic vector variable whose names is var_name. 

            Args:
                var_name: Name of the array variable.
                size: Size of the array variable.
        """
        _require_identifier("var_name", var_name)
        _require_integer("size", size)
        defined_names = {var.name for var in self.__scalar_vars + self.__array_vars}
        if var_name in defined_names:
            raise ValueError(f"var_name must be unique; {var_name!r} is already defined")
        array_var = symbols(var_name+'[0:%d]' %(size))
        self.__array_vars.append(ArrayVariable(array_var, var_name, size, []))
        return array_var

    def set_FB_epsilon(self, FB_epsilon: Sequence[float]) -> None:
        """ Set reguralization term of the semi-smooth Fischer-Burumeister (FB) 
            method. Set the array whose size is dimension of the inequality 
            constraints considered by semi-smooth FB method.

            Args:
                FB epsilon: Array of the reguralization term. 
        """
        _require_sequence("FB_epsilon", FB_epsilon)
        if self.__symbolic_functions is not None:
            _require_length("FB_epsilon", FB_epsilon, self.__nh)
        for index, eps in enumerate(FB_epsilon):
            _require_finite_number(f"FB_epsilon[{index}]", eps, minimum=0.0)
        self.__FB_epsilon = list(FB_epsilon)

    def set_scalar_var(self, name: str, value: ModelParameter) -> None:
        """ Set the value of the scalar variable you defied. 

            Args:
                name: Name of the scalar variable.
                value: Value of the scalar variable.
        """
        _require_identifier("name", name)
        _require_model_parameter("value", value)
        for scalar_var in self.__scalar_vars:
            if name == scalar_var.name:
                scalar_var.value = value
                return
        raise ValueError(f"unknown scalar variable {name!r}; define it before setting it")

    def set_scalar_vars(
        self, *name_and_value_list: Sequence[ModelParameter]
    ) -> None:
        """ Set the values of the scalar variables you defied. 

            Args:
                name_and_value_lis: A list composed of the name of 
                the scalar variable and value of the scalar variable.
        """
        for name_and_value in name_and_value_list:
            _require_length("name_and_value", name_and_value, 2)
            name = name_and_value[0]
            if not isinstance(name, str):
                raise TypeError(
                    "name_and_value[0] must be a string; "
                    f"got {type(name).__name__}"
                )
            value = name_and_value[1]
            self.set_scalar_var(name, value)

    def set_array_var(self, name: str, values: Sequence[ModelParameter]) -> None:
        """ Set the value of the array variable you defied. 

            Args:
                name: Name of the arrray variable.
                values: Values of the arry variable. This size is used as  
                        the size of the array variable.
        """
        _require_identifier("name", name)
        _require_sequence("values", values)
        for array_var in self.__array_vars:
            if name == array_var.name:
                _require_length(f"values for {name!r}", values, array_var.size)
                for index, value in enumerate(values):
                    _require_model_parameter(f"values[{index}]", value)
                array_var.values = list(values)
                return
        raise ValueError(f"unknown array variable {name!r}; define it before setting it")

    def set_functions(
        self,
        f: Sequence[SymbolicExpression],
        C: Sequence[SymbolicExpression],
        h: Sequence[SymbolicExpression],
        L: SymbolicExpression,
        phi: SymbolicExpression,
    ) -> None:
        """ Sets functions that defines the optimal control problem.

            Args: 
                f: The state equation. The dimension must be nx.
                C: The equality consrtaints. If there are no equality 
                    constraints, set the empty list.
                h: The inequality consrtaints considered by semi-smooth 
                    Fischer-Burumeister method. If there are no such inequality 
                    constraints, set the empty list.
                L: The stage cost.
                phi: The terminal cost.
        """
        _require_length("f", f, self.__nx)
        _require_sequence("C", C)
        _require_sequence("h", h)
        if L is None:
            raise ValueError("L must be a scalar symbolic expression; got None")
        if phi is None:
            raise ValueError("phi must be a scalar symbolic expression; got None")
        self.__nc = len(C)
        self.__nh = len(h)
        x = symbols('x[0:%d]' %(self.__nx))
        u = symbols('u[0:%d]' %(self.__nu+self.__nc+self.__nh))
        lmd = symbols('lmd[0:%d]' %(self.__nx))
        hamiltonian = L + sum(lmd[i] * f[i] for i in range(self.__nx))
        hamiltonian += sum(u[self.__nu+i] * C[i] for i in range(self.__nc))
        nuc = self.__nu + self.__nc
        hamiltonian += sum(u[nuc+i] * h[i] for i in range(self.__nh))
        hx = symutils.diff_scalar_func(hamiltonian, x)
        hu = symutils.diff_scalar_func(hamiltonian, u)
        fb_eps = symbols('fb_eps[0:%d]' %(self.__nh))
        for i in range(self.__nh):
            hu[nuc+i] = sqrt(u[nuc+i]**2 + h[i]**2 + fb_eps[i]) - (u[nuc+i] - h[i])
        phix = symutils.diff_scalar_func(phi, x)
        self.__symbolic_functions = SymbolicFunctions(f, phix, hx, hu)

    def add_control_input_bounds(
        self, uindex: int, umin: float, umax: float, dummy_weight: float
        ) -> None:
        """ Adds the bax constraints on the control input that is condensed in 
            linear problem. 

            Args: 
                uindex: The index of the constrianed control input element. 
                umin: The minimum value of the constrianed control input. 
                umax: The minimum value of the constrianed control input. 
                dummy_weight: An weight to stabilize the numerical computation.
        """
        _require_integer("uindex", uindex, minimum=0)
        if uindex >= self.__nu:
            raise ValueError(
                f"uindex must be between 0 and {self.__nu - 1}; got {uindex}"
            )
        _require_finite_number("umin", umin)
        _require_finite_number("umax", umax)
        if umin >= umax:
            raise ValueError(f"umin must be less than umax; got {umin!r} >= {umax!r}")
        _require_finite_number("dummy_weight", dummy_weight, minimum=0.0)
        find_same_index = False
        for ub in self.__ubounds:
            if ub.uindex == uindex:
                find_same_index = True
                ub.umin = umin
                ub.umax = umax
                ub.dummy_weight = dummy_weight
        if not find_same_index:
            self.__ubounds.append(ControlInputBound(uindex, umin, umax, dummy_weight))

    def set_nlp_type(self, nlp_type: NLPType) -> None:
        """ Sets solver types of the C/GMRES methods. 

            Args: 
                nlp_type: The solver type. Choose from 
                NLPType.SingleShooting and NLPType.MultipleShooting, 
        """
        if not isinstance(nlp_type, NLPType):
            raise TypeError(
                "nlp_type must be an NLPType value "
                f"(SingleShooting or MultipleShooting); got {nlp_type!r}"
            )
        self.__nlp_type = nlp_type

    def set_horizon_params(self, Tf: float, alpha: float = 0.0) -> None:
        """ Sets parameters of the horizon of NMPC. If alpha > 0, then the 
            length of the horzion at time t is given by Tf * (1-exp(-alpha*t)). 
            If alpha is not positive, the it is given by Tf.

            Args: 
                Tf, alpha: Parameter about the length of the horizon of NMPC.
        """
        _require_finite_number("Tf", Tf, minimum=0.0, strict=True)
        _require_finite_number("alpha", alpha)
        self.__horizon_params = HorizonParams(Tf, alpha)

    def set_solver_params(
            self, sampling_time: float, N: int, finite_difference_epsilon: float,
            zeta: float, kmax: int
        ) -> None:
        """ Sets parameters of the NMPC solvers based on the C/GMRES method. 

            Args: 
                sampling_time: The sampling period of NMPC
                N: The number of the grid for the discretization
                    of the horizon of NMPC.
                finite_difference_epsilon: The small positive value for 
                    finitei difference approximation used in the FD-GMRES. 
                zeta: A stabilization parameter of the C/GMRES method. It may 
                    work well if you set as zeta=1/sampling_period.
                kmax: Maximam number of the iteration of the Krylov 
                    subspace method for the linear problem. 
        """
        _require_finite_number(
            "sampling_time", sampling_time, minimum=0.0, strict=True
        )
        _require_integer("N", N)
        _require_finite_number(
            "finite_difference_epsilon",
            finite_difference_epsilon,
            minimum=0.0,
            strict=True,
        )
        _require_finite_number("zeta", zeta, minimum=0.0, strict=True)
        _require_integer("kmax", kmax)
        self.__solver_params = SolverParams(sampling_time, N, finite_difference_epsilon, zeta, kmax)

    def set_initialization_params(
            self, solution_initial_guess: Sequence[float], tolerance: float=1.0e-04,
            max_iterations: int=100
        ) -> None:
        """ Set parameters for the initialization of the C/GMRES solvers. 

            Args: 
                solution_initial_guess: The initial guess of the solution of the 
                    initialization. Size must be the nu + dimensions of C and 
                    h.
                torelance: The residual torelance of the 
                    initialization solved by Newton's method. The Newton 
                    iteration terminates if the optimality error is smaller than 
                    this value.
                max_iteration: The maximum number of the Newton iteration. 
        """
        expected_size = self.__nu + self.__nc + self.__nh
        _require_length("solution_initial_guess", solution_initial_guess, expected_size)
        for index, value in enumerate(solution_initial_guess):
            _require_finite_number(f"solution_initial_guess[{index}]", value)
        _require_finite_number("tolerance", tolerance, minimum=0.0)
        _require_integer("max_iterations", max_iterations, minimum=0)
        self.__initialization_params = InitializationParams(
            list(solution_initial_guess), tolerance, max_iterations
        )

    def set_simulation_params(
            self, initial_time: float, initial_state: Sequence[float],
            simulation_length: float
        ) -> None:
        """ Set parameters for numerical simulation. 

            Args: 
                initial_time: The time parameter at the beginning of the 
                    simulation. 
                initial_state: The state of the system at the beginning of the 
                    simulation. 
                simulation_length: The length of the numerical simulation. 
        """
        _require_finite_number("initial_time", initial_time)
        _require_length("initial_state", initial_state, self.__nx)
        for index, value in enumerate(initial_state):
            _require_finite_number(f"initial_state[{index}]", value)
        _require_finite_number(
            "simulation_length", simulation_length, minimum=0.0, strict=True
        )
        self.__simulation_params = SimulationParams(
            initial_time, list(initial_state), simulation_length
        )

    def generate_ocp_definition(
        self,
        simplification: bool = False,
        common_subexpression_elimination: bool = False,
    ) -> None:
        """ Generates the C++ source file in which the equations to solve the 
            optimal control problem are described. Before call this method, 
            set_functions() must be called.

            Args: 
                simplification: The flag for simplification. If True, the 
                    Symbolic functions are simplified. Default is False.
                common_subexpression_elimination: The flag for common subexpression elimination. If True, 
                    common subexpressions are eliminated. Default is False.
        """
        _require_boolean("simplification", simplification)
        _require_boolean(
            "common_subexpression_elimination", common_subexpression_elimination
        )
        self._require_configured(
            "generate_ocp_definition()",
            [(self.__symbolic_functions, "set_functions")],
        )
        assert self.__symbolic_functions is not None
        unset_arrays = [
            array_var.name
            for array_var in self.__array_vars
            if len(array_var.values) != array_var.size
        ]
        if unset_arrays:
            names = ", ".join(repr(name) for name in unset_arrays)
            raise RuntimeError(
                "generate_ocp_definition() requires values for array variables "
                f"{names}; call set_array_var() first"
            )
        if self.__nh > 0:
            _require_length("FB_epsilon", self.__FB_epsilon, self.__nh)
        os.makedirs(self.get_ocp_pybind_dir(), exist_ok=True)
        os.makedirs(os.path.join(self.get_ocp_pybind_dir(), self.__ocp_name), exist_ok=True)
        os.makedirs(os.path.join(self.get_ocp_pybind_dir(), 'common'), exist_ok=True)
        if simplification:
            symutils.simplify(self.__symbolic_functions.f)
            symutils.simplify(self.__symbolic_functions.hx)
            symutils.simplify(self.__symbolic_functions.hu)
            symutils.simplify(self.__symbolic_functions.phix)
        f_model_h = open(
            os.path.join(self.get_ocp_dir(), "ocp.hpp"),
            "w",
            encoding="utf-8",
            newline="\n",
        )
        f_model_h.write('// This file was automatically generated by autogenu-jupyter (https://github.com/ohtsukalab/autogenu-jupyter). \n')
        f_model_h.write('// The autogenu-jupyter copyright holders make no ownership claim of its contents. \n\n')
        f_model_h.write(
            '#ifndef CGMRES__OCP_'+str(self.__ocp_name).upper()+'_HPP_ \n'
        )
        f_model_h.write(
            '#define CGMRES__OCP_'+str(self.__ocp_name).upper()+'_HPP_ \n'
        )
        f_model_h.writelines([
""" 
#define _USE_MATH_DEFINES

#include <cmath>
#include <array>
#include <iostream>

#include "cgmres/types.hpp"
#include "cgmres/detail/macros.hpp"

namespace cgmres {

/// 
""" 
        ])
        f_model_h.write('/// @class OCP_'+self.__ocp_name+'\n')
        f_model_h.write('/// @brief Definition of the optimal control problem (OCP) of '+self.__ocp_name+'.\n')
        f_model_h.write('/// \n')
        f_model_h.write('class OCP_'+self.__ocp_name+' {')
        f_model_h.writelines([
""" 
public:
  ///
  /// @brief Dimension of the state. 
  ///
"""
        ])
        f_model_h.write(
            '  static constexpr int nx = '+str(self.__nx)+';\n'
        )
        f_model_h.writelines([
""" 
  ///
  /// @brief Dimension of the control input. 
  ///
"""
        ])
        f_model_h.write(
            '  static constexpr int nu = '
            +str(self.__nu)+';\n'
        )
        f_model_h.writelines([
""" 
  ///
  /// @brief Dimension of the equality constraints. 
  ///
"""
        ])
        f_model_h.write(
            '  static constexpr int nc = '
            +str(self.__nc+self.__nh)+';\n'
        )
        f_model_h.writelines([
""" 
  ///
  /// @brief Dimension of the Fischer-Burmeister function (already counded in nc). 
  ///
"""
        ])
        f_model_h.write(
            '  static constexpr int nh = '
            +str(self.__nh)+';\n'
        )
        f_model_h.writelines([
""" 
  ///
  /// @brief Dimension of the concatenation of the control input and equality constraints. 
  ///
  static constexpr int nuc = nu + nc;

  ///
  /// @brief Dimension of the bound constraints on the control input. 
  ///
"""
        ])
        f_model_h.write(
            '  static constexpr int nub = '
            +str(len(self.__ubounds))+';\n\n'
        )
        f_model_h.writelines([
            '  double '+scalar_var.name+' = '
            +str(scalar_var.value)+';\n' for scalar_var in self.__scalar_vars
        ])
        f_model_h.write('\n')
        for array_var in self.__array_vars:
            f_model_h.write(
                '  std::array<double, '+str(array_var.size)+'> '+array_var.name+' = {'
            )
            for i in range(array_var.size-1):
                f_model_h.write(str(array_var.values[i])+', ')
            f_model_h.write(str(array_var.values[array_var.size-1])+'};\n')
        if len(self.__ubounds) > 0:
            nub = len(self.__ubounds)
            f_model_h.write('\n  static constexpr std::array<int, nub> ubound_indices = {')
            for i in range(nub-1):
                f_model_h.write(str(self.__ubounds[i].uindex)+', ')
            f_model_h.write(str(self.__ubounds[nub-1].uindex)+'};\n')
            f_model_h.write('  std::array<double, nub> umin = {')
            for i in range(nub-1):
                f_model_h.write(str(self.__ubounds[i].umin)+', ')
            f_model_h.write(str(self.__ubounds[nub-1].umin)+'};\n')
            f_model_h.write('  std::array<double, nub> umax = {')
            for i in range(nub-1):
                f_model_h.write(str(self.__ubounds[i].umax)+', ')
            f_model_h.write(str(self.__ubounds[nub-1].umax)+'};\n')
            f_model_h.write('  std::array<double, nub> dummy_weight = {')
            for i in range(nub-1):
                f_model_h.write(str(self.__ubounds[i].dummy_weight)+', ')
            f_model_h.write(str(self.__ubounds[nub-1].dummy_weight)+'};\n')
        if self.__nh > 0:
            f_model_h.write('\n  std::array<double, nh> fb_eps = {')
            for i in range(self.__nh-1):
                f_model_h.write(str(self.__FB_epsilon[i])+', ')
            f_model_h.write(str(self.__FB_epsilon[self.__nh-1])+'};\n')
        f_model_h.write('\n  void disp(std::ostream& os) const {\n')
        f_model_h.write('    os << "OCP_'+self.__ocp_name+':" << std::endl;\n')
        f_model_h.write('    os << "  nx:  " << nx << std::endl;\n')
        f_model_h.write('    os << "  nu:  " << nu << std::endl;\n')
        f_model_h.write('    os << "  nc:  " << nc << std::endl;\n')
        f_model_h.write('    os << "  nh:  " << nh << std::endl;\n')
        f_model_h.write('    os << "  nuc: " << nuc << std::endl;\n')
        f_model_h.write('    os << "  nub: " << nub << std::endl;\n')
        f_model_h.write('    os << std::endl;\n')
        f_model_h.writelines([
            '    os << "  '+scalar_var.name+': " << '+scalar_var.name+' << std::endl;\n' for scalar_var in self.__scalar_vars
        ])
        f_model_h.write('    os << std::endl;\n')
        f_model_h.write('    Eigen::IOFormat fmt(4, 0, ", ", "", "[", "]");\n')
        f_model_h.write('    Eigen::IOFormat intfmt(1, 0, ", ", "", "[", "]");\n')
        f_model_h.writelines([
            '    os << "  '+array_var.name+': " << Map<const VectorX>('+array_var.name+'.data(), '+array_var.name+'.size()).transpose().format(fmt) << std::endl;\n' for array_var in self.__array_vars
        ])
        if len(self.__ubounds) > 0:
            nub = len(self.__ubounds)
            f_model_h.write('    os << std::endl;\n')
            f_model_h.write('    os << "  ubound_indices: " << Map<const VectorXi>(ubound_indices.data(), ubound_indices.size()).transpose().format(intfmt) << std::endl;\n')
            f_model_h.write('    os << "  umin: " << Map<const VectorX>(umin.data(), umin.size()).transpose().format(fmt) << std::endl;\n')
            f_model_h.write('    os << "  umax: " << Map<const VectorX>(umax.data(), umax.size()).transpose().format(fmt) << std::endl;\n')
            f_model_h.write('    os << "  dummy_weight: " << Map<const VectorX>(dummy_weight.data(), dummy_weight.size()).transpose().format(fmt) << std::endl;\n')

        if self.__nh > 0:
            f_model_h.write('    os << std::endl;\n')
            f_model_h.write('    os << "  fb_eps: " << Map<const VectorX>(fb_eps.data(), fb_eps.size()).transpose().format(fmt) << std::endl;\n')
        f_model_h.write('  }\n\n')
        f_model_h.write('  friend std::ostream& operator<<(std::ostream& os, const OCP_'+self.__ocp_name+'& ocp) { \n')
        f_model_h.write('    ocp.disp(os);\n')
        f_model_h.write('    return os;\n')
        f_model_h.write('  }\n\n')
        f_model_h.writelines([
"""
  ///
  /// @brief Synchrozies the internal parameters of this OCP with the external references.
  /// This method is called at the beginning of each MPC update.
  ///
  void synchronize() {
  }

  ///
  /// @brief Computes the state equation dx = f(t, x, u).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Control input.
  /// @param[out] dx Evaluated value of the state equation.
  /// @remark This method is intended to be used inside of the cgmres solvers and does not check size of each argument. 
  /// Use the overloaded method if you call this outside of the cgmres solvers. 
  ///
  void eval_f(const double t, const double* x, const double* u, 
              double* dx) const {
""" 
        ])
        symutils.write_symfunc(f_model_h, self.__symbolic_functions.f, 'dx', common_subexpression_elimination)
        f_model_h.writelines([
""" 
  }

  ///
  /// @brief Computes the partial derivative of terminal cost with respect to state, 
  /// i.e., phix = dphi/dx(t, x).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[out] phix Evaluated value of the partial derivative of terminal cost.
  /// @remark This method is intended to be used inside of the cgmres solvers and does not check size of each argument. 
  /// Use the overloaded method if you call this outside of the cgmres solvers. 
  ///
  void eval_phix(const double t, const double* x, double* phix) const {
""" 
        ])
        symutils.write_symfunc(f_model_h, self.__symbolic_functions.phix, 'phix', common_subexpression_elimination)
        f_model_h.writelines([
""" 
  }

  ///
  /// @brief Computes the partial derivative of the Hamiltonian with respect to state, 
  /// i.e., hx = dH/dx(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. 
  /// @param[in] lmd Costate. 
  /// @param[out] hx Evaluated value of the partial derivative of the Hamiltonian.
  /// @remark This method is intended to be used inside of the cgmres solvers and does not check size of each argument. 
  /// Use the overloaded method if you call this outside of the cgmres solvers. 
  ///
  void eval_hx(const double t, const double* x, const double* u, 
               const double* lmd, double* hx) const {
""" 
        ])
        symutils.write_symfunc(f_model_h, self.__symbolic_functions.hx, 'hx', common_subexpression_elimination)
        f_model_h.writelines([
""" 
  }

  ///
  /// @brief Computes the partial derivative of the Hamiltonian with respect to control input and the equality constraints, 
  /// i.e., hu = dH/du(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State.
  /// @param[in] u Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. 
  /// @param[in] lmd Costate. 
  /// @param[out] hu Evaluated value of the partial derivative of the Hamiltonian.
  /// @remark This method is intended to be used inside of the cgmres solvers and does not check size of each argument. 
  /// Use the overloaded method if you call this outside of the cgmres solvers. 
  ///
  void eval_hu(const double t, const double* x, const double* u, 
               const double* lmd, double* hu) const {
""" 
        ])
        symutils.write_symfunc(f_model_h, self.__symbolic_functions.hu, 'hu', common_subexpression_elimination)
        f_model_h.writelines([
""" 
  }

  ///
  /// @brief Computes the state equation dx = f(t, x, u).
  /// @param[in] t Time.
  /// @param[in] x State. Size must be nx.
  /// @param[in] u Control input. Size must be nu.
  /// @param[out] dx Evaluated value of the state equation. Size must be nx.
  ///
  template <typename VectorType1, typename VectorType2, typename VectorType3>
  void eval_f(const double t, const MatrixBase<VectorType1>& x, 
              const MatrixBase<VectorType2>& u, 
              const MatrixBase<VectorType3>& dx) const {
    if (x.size() != nx) {
      throw std::invalid_argument("[OCP]: x.size() must be " + std::to_string(nx));
    }
    if (u.size() != nu) {
      throw std::invalid_argument("[OCP]: u.size() must be " + std::to_string(nu));
    }
    if (dx.size() != nx) {
      throw std::invalid_argument("[OCP]: dx.size() must be " + std::to_string(nx));
    }
    eval_f(t, x.derived().data(), u.derived().data(), CGMRES_EIGEN_CONST_CAST(VectorType3, dx).data());
  }

  ///
  /// @brief Computes the partial derivative of terminal cost with respect to state, 
  /// i.e., phix = dphi/dx(t, x).
  /// @param[in] t Time.
  /// @param[in] x State. Size must be nx.
  /// @param[out] phix Evaluated value of the partial derivative of terminal cost. Size must be nx.
  ///
  template <typename VectorType1, typename VectorType2>
  void eval_phix(const double t, const MatrixBase<VectorType1>& x, 
                 const MatrixBase<VectorType2>& phix) const {
    if (x.size() != nx) {
      throw std::invalid_argument("[OCP]: x.size() must be " + std::to_string(nx));
    }
    if (phix.size() != nx) {
      throw std::invalid_argument("[OCP]: phix.size() must be " + std::to_string(nx));
    }
    eval_phix(t, x.derived().data(), CGMRES_EIGEN_CONST_CAST(VectorType2, phix).data());
  }

  ///
  /// @brief Computes the partial derivative of the Hamiltonian with respect to the state, 
  /// i.e., hx = dH/dx(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State. Size must be nx.
  /// @param[in] uc Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. Size must be nuc. 
  /// @param[in] lmd Costate.  Size must be nx.
  /// @param[out] hx Evaluated value of the partial derivative of the Hamiltonian. Size must be nx.
  ///
  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_hx(const double t, const MatrixBase<VectorType1>& x, 
               const MatrixBase<VectorType2>& uc, 
               const MatrixBase<VectorType3>& lmd, 
               const MatrixBase<VectorType4>& hx) const {
    if (x.size() != nx) {
      throw std::invalid_argument("[OCP]: x.size() must be " + std::to_string(nx));
    }
    if (uc.size() != nuc) {
      throw std::invalid_argument("[OCP]: uc.size() must be " + std::to_string(nuc));
    }
    if (lmd.size() != nx) {
      throw std::invalid_argument("[OCP]: lmd.size() must be " + std::to_string(nx));
    }
    if (hx.size() != nuc) {
      throw std::invalid_argument("[OCP]: hx.size() must be " + std::to_string(nx));
    }
    eval_hx(t, x.derived().data(), uc.derived().data(), lmd.derived().data(), CGMRES_EIGEN_CONST_CAST(VectorType4, hx).data());
  }

  ///
  /// @brief Computes the partial derivative of the Hamiltonian with respect to control input and the equality constraints, 
  /// i.e., hu = dH/du(t, x, u, lmd).
  /// @param[in] t Time.
  /// @param[in] x State. Size must be nx.
  /// @param[in] uc Concatenatin of the control input and Lagrange multiplier with respect to the equality constraints. Size must be nuc. 
  /// @param[in] lmd Costate. Size must be nx. 
  /// @param[out] hu Evaluated value of the partial derivative of the Hamiltonian. Size must be nuc.
  ///
  template <typename VectorType1, typename VectorType2, typename VectorType3, typename VectorType4>
  void eval_hu(const double t, const MatrixBase<VectorType1>& x, 
               const MatrixBase<VectorType2>& uc, 
               const MatrixBase<VectorType3>& lmd, 
               const MatrixBase<VectorType4>& hu) const {
    if (x.size() != nx) {
      throw std::invalid_argument("[OCP]: x.size() must be " + std::to_string(nx));
    }
    if (uc.size() != nuc) {
      throw std::invalid_argument("[OCP]: uc.size() must be " + std::to_string(nuc));
    }
    if (lmd.size() != nx) {
      throw std::invalid_argument("[OCP]: lmd.size() must be " + std::to_string(nx));
    }
    if (hu.size() != nuc) {
      throw std::invalid_argument("[OCP]: hu.size() must be " + std::to_string(nuc));
    }
    eval_hu(t, x.derived().data(), uc.derived().data(), lmd.derived().data(), CGMRES_EIGEN_CONST_CAST(VectorType4, hu).data());
  }

};

} // namespace cgmres

#endif // CGMRES_OCP_HPP_
""" 
        ])
        f_model_h.close()
        print('\'ocp.hpp\', the definition of the OCP, is generated at', self.get_ocp_dir())

    def generate_main(self) -> None:
        """Generate the closed-loop simulation source from a packaged template."""
        self._require_configured(
            "generate_main()",
            [
                (self.__symbolic_functions, "set_functions"),
                (self.__nlp_type, "set_nlp_type"),
                (self.__horizon_params, "set_horizon_params"),
                (self.__solver_params, "set_solver_params"),
                (self.__initialization_params, "set_initialization_params"),
                (self.__simulation_params, "set_simulation_params"),
            ],
        )
        assert self.__symbolic_functions is not None
        assert self.__nlp_type is not None
        assert self.__horizon_params is not None
        assert self.__solver_params is not None
        assert self.__initialization_params is not None
        assert self.__simulation_params is not None
        solution_size = self.__nu + self.__nc + self.__nh
        _require_length(
            "solution_initial_guess",
            self.__initialization_params.solution_initial_guess,
            solution_size,
        )
        _require_length(
            "initial_state", self.__simulation_params.initial_state, self.__nx
        )
        if self.__nh > 0:
            _require_length("FB_epsilon", self.__FB_epsilon, self.__nh)

        if self.__nlp_type == NLPType.SingleShooting:
            solver_header = "single_shooting_cgmres_solver.hpp"
            solver_initialization = (
                "  cgmres::SingleShootingCGMRESSolver<cgmres::OCP_"
                + self.__ocp_name
                + ", N, kmax> mpc(ocp, horizon, settings);\n"
                "  mpc.set_uc(initializer.ucopt());\n"
                "  mpc.init_dummy_mu();"
            )
        elif self.__nlp_type == NLPType.MultipleShooting:
            solver_header = "multiple_shooting_cgmres_solver.hpp"
            solver_initialization = (
                "  cgmres::MultipleShootingCGMRESSolver<cgmres::OCP_"
                + self.__ocp_name
                + ", N, kmax> mpc(ocp, horizon, settings);\n"
                "  mpc.set_uc(initializer.ucopt());\n"
                "  mpc.init_x_lmd(t0, x0);\n"
                "  mpc.init_dummy_mu();"
            )
        else:
            raise NotImplementedError("Unsupported NLP type")

        write_generated_file(
            os.path.join(self.get_ocp_dir(), "main.cpp"),
            "main.cpp.in",
            solver_header=solver_header,
            ocp_name=self.__ocp_name,
            Tf=self.__horizon_params.Tf,
            alpha=self.__horizon_params.alpha,
            sampling_time=self.__solver_params.sampling_time,
            zeta=self.__solver_params.zeta,
            finite_difference_epsilon=self.__solver_params.finite_difference_epsilon,
            max_iter=self.__initialization_params.max_iteraions,
            opterr_tol=self.__initialization_params.tolerance,
            initial_time=self.__simulation_params.initial_time,
            state_size=len(self.__simulation_params.initial_state),
            initial_state=", ".join(map(str, self.__simulation_params.initial_state)),
            kmax_init=min(self.__solver_params.kmax, solution_size),
            solution_size=solution_size,
            solution_initial_guess=", ".join(
                map(str, self.__initialization_params.solution_initial_guess)
            ),
            N=self.__solver_params.N,
            kmax=min(self.__solver_params.kmax, self.__solver_params.N * solution_size),
            solver_initialization=solver_initialization,
            simulation_length=self.__simulation_params.simulation_length,
        )
        print("'main.cpp', the closed-loop simulation code, is generated at", self.get_ocp_dir())

    def generate_python_bindings(self) -> None:
        self._require_configured(
            "generate_python_bindings()",
            [
                (self.__symbolic_functions, "set_functions"),
                (self.__solver_params, "set_solver_params"),
            ],
        )
        assert self.__symbolic_functions is not None
        assert self.__solver_params is not None
        f_pybind11 = open(
            os.path.join(self.get_ocp_pybind_dir(), self.__ocp_name, "ocp.cpp"),
            "w",
            encoding="utf-8",
            newline="\n",
        )
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/ohtsukalab/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/types.hpp"
#include "ocp.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

""" 
        ])
        f_pybind11.write('using OCP = OCP_'+str(self.__ocp_name)+';\n')
        f_pybind11.writelines([
"""
PYBIND11_MODULE(ocp, m) { 
  py::class_<OCP>(m, "OCP")
    .def(py::init<>())  
    .def("clone", [](const OCP& self) { 
       auto copy = self; 
       return copy; 
     }) 
    .def("synchronize", &OCP::synchronize)
    .def("eval_f", [](const OCP& self, const Scalar t,  
                      const VectorX& x, const VectorX& u) { 
        Vector<OCP::nx> dx(Vector<OCP::nx>::Zero());
        self.eval_f(t, x, u, dx); 
        return dx;
     }, py::arg("t"), py::arg("x"), py::arg("u"))
    .def("eval_phix", [](const OCP& self, const Scalar t, const VectorX& x) {
        Vector<OCP::nx> phix(Vector<OCP::nx>::Zero());
        self.eval_phix(t, x, phix);
        return phix;
     }, py::arg("t"), py::arg("x"))
    .def("eval_hx", [](const OCP& self, const Scalar t, 
                       const VectorX& x, const VectorX& u, const VectorX& lmd) {
        Vector<OCP::nx> hx(Vector<OCP::nx>::Zero());
        self.eval_hx(t, x, u, lmd, hx);
        return hx;
     }, py::arg("t"), py::arg("x"), py::arg("u"), py::arg("lmd"))
    .def("eval_hu", [](const OCP& self, const Scalar t,
                       const VectorX& x, const VectorX& u, const VectorX& lmd) {
        Vector<OCP::nuc> hu(Vector<OCP::nuc>::Zero());
        self.eval_hu(t, x, u, lmd, hu);
        return hu;
     }, py::arg("t"), py::arg("x"), py::arg("u"), py::arg("lmd"))
""" 
        ])
        for scalar_var in self.__scalar_vars:
            f_pybind11.write('    .def_readwrite("'+scalar_var.name+'", &OCP::'+scalar_var.name+')\n')
        for array_var in self.__array_vars:
            name = array_var.name
            size = array_var.size
            f_pybind11.write('    .def_property("'+name+'", \n')
            f_pybind11.write('      [](const OCP& self) { return Map<const VectorX>(self.'+name+'.data(), self.'+name+'.size()); },\n')
            f_pybind11.write('      [](OCP& self, const VectorX& v) { \n')
            f_pybind11.write('        if (v.size() != '+str(size)+') { \n')
            f_pybind11.write('          throw std::invalid_argument("[OCP]: \''+name+'.size()\' must be "+std::to_string('+str(size)+')); \n')
            f_pybind11.write('        } Map<VectorX>(self.'+name+'.data(), self.'+name+'.size()) = v; })\n')
        if len(self.__ubounds) > 0:
            f_pybind11.writelines([
"""
    .def_readonly_static("ubound_indices", &OCP::ubound_indices)
    .def_property("umin", 
      [](const OCP& self) { return Map<const VectorX>(self.umin.data(), self.umin.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != static_cast<Eigen::Index>(self.umin.size())) {
          throw std::invalid_argument("[OCP]: 'umin.size()' must be "+std::to_string(self.umin.size()));
        } Map<VectorX>(self.umin.data(), self.umin.size()) = v; })
    .def_property("umax", 
      [](const OCP& self) { return Map<const VectorX>(self.umax.data(), self.umax.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != static_cast<Eigen::Index>(self.umax.size())) {
          throw std::invalid_argument("[OCP]: 'umax.size()' must be "+std::to_string(self.umax.size()));
        } Map<VectorX>(self.umax.data(), self.umax.size()) = v; })
    .def_property("dummy_weight", 
      [](const OCP& self) { return Map<const VectorX>(self.dummy_weight.data(), self.dummy_weight.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != static_cast<Eigen::Index>(self.dummy_weight.size())) {
          throw std::invalid_argument("[OCP]: 'dummy_weight.size()' must be "+std::to_string(self.dummy_weight.size()));
        } Map<VectorX>(self.dummy_weight.data(), self.dummy_weight.size()) = v; })
""" 
            ])
        if self.__nh > 0:
            f_pybind11.writelines([
"""
    .def_property("fb_eps", 
      [](const OCP& self) { return Map<const VectorX>(self.fb_eps.data(), self.fb_eps.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != self.fb_eps.size()) { 
          throw std::invalid_argument("[OCP]: 'fb_eps.size()' must be "+std::to_string(self.fb_eps.size()));
        } Map<VectorX>(self.fb_eps.data(), self.fb_eps.size()) = v; })
""" 
            ])
        f_pybind11.writelines([
"""
    .def_readonly_static("nx", &OCP::nx)
    .def_readonly_static("nu", &OCP::nu)
    .def_readonly_static("nc", &OCP::nc)
    .def_readonly_static("nh", &OCP::nh)
    .def_readonly_static("nuc", &OCP::nuc)
    .def_readonly_static("nub", &OCP::nub)
    .def("__str__", [](const OCP& self) { 
        std::stringstream ss; 
        ss << self; 
        return ss.str(); 
      }); 
}

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()

        solution_size = self.__nc + self.__nu + self.__nh
        solver_modules = [
            (
                "zero_horizon_ocp_solver",
                "cgmres/zero_horizon_ocp_solver.hpp",
                "cgmres/python/zero_horizon_ocp_solver.hpp",
                "constexpr int kmax_init = "
                + str(min(self.__solver_params.kmax, solution_size))
                + ";\n",
                "bind_zero_horizon_ocp_solver<OCP_"
                + self.__ocp_name
                + ", kmax_init>(m);",
            ),
            (
                "single_shooting_cgmres_solver",
                "cgmres/single_shooting_cgmres_solver.hpp",
                "cgmres/python/single_shooting_cgmres_solver.hpp",
                "constexpr int N = "
                + str(self.__solver_params.N)
                + ";\nconstexpr int kmax = "
                + str(min(self.__solver_params.kmax, self.__solver_params.N * solution_size))
                + ";\n",
                "bind_single_shooting_cgmres_solver<OCP_"
                + self.__ocp_name
                + ", N, kmax>(m);",
            ),
            (
                "multiple_shooting_cgmres_solver",
                "cgmres/multiple_shooting_cgmres_solver.hpp",
                "cgmres/python/multiple_shooting_cgmres_solver.hpp",
                "constexpr int N = "
                + str(self.__solver_params.N)
                + ";\nconstexpr int kmax = "
                + str(min(self.__solver_params.kmax, self.__solver_params.N * solution_size))
                + ";\n",
                "bind_multiple_shooting_cgmres_solver<OCP_"
                + self.__ocp_name
                + ", N, kmax>(m);",
            ),
        ]
        for module_name, core_header, binding_header, constants, bind_call in solver_modules:
            write_generated_file(
                os.path.join(self.get_ocp_pybind_dir(), self.__ocp_name, module_name + ".cpp"),
                "binding_module.cpp.in",
                core_header=core_header,
                binding_header=binding_header,
                ocp_include='#include "ocp.hpp"\n',
                constants=constants,
                module_name=module_name,
                bind_call=bind_call,
            )

        common_modules = [
            (
                "horizon",
                "cgmres/horizon.hpp",
                "cgmres/python/horizon.hpp",
                "bind_horizon(m);",
            ),
            (
                "solver_settings",
                "cgmres/solver_settings.hpp",
                "cgmres/python/solver_settings.hpp",
                "bind_solver_settings(m);",
            ),
            ("timer", "cgmres/timer.hpp", "cgmres/python/timer.hpp", "bind_timer(m);"),
        ]
        for module_name, core_header, binding_header, bind_call in common_modules:
            write_generated_file(
                os.path.join(self.get_ocp_pybind_dir(), "common", module_name + ".cpp"),
                "binding_module.cpp.in",
                core_header=core_header,
                binding_header=binding_header,
                ocp_include="",
                constants="",
                module_name=module_name,
                bind_call=bind_call,
            )

        write_generated_file(
            os.path.join(self.get_ocp_pybind_dir(), self.__ocp_name, "__init__.py"),
            "package_init.py.in",
            imports=(
                "from .ocp import *\n"
                "from .zero_horizon_ocp_solver import *\n"
                "from .single_shooting_cgmres_solver import *\n"
                "from .multiple_shooting_cgmres_solver import *"
            ),
        )
        write_generated_file(
            os.path.join(self.get_ocp_pybind_dir(), "common", "__init__.py"),
            "package_init.py.in",
            imports=(
                "from .horizon import *\n"
                "from .solver_settings import *\n"
                "from .timer import *"
            ),
        )
        print("pybind11 source codes are generated at", self.get_ocp_pybind_dir())



    def generate_cmake(self) -> None:
        """Generate CMake project files from packaged templates."""
        write_generated_file(
            os.path.join(self.get_ocp_dir(), "CMakeLists.txt"),
            "CMakeLists.txt.in",
            ocp_name=self.__ocp_name,
        )
        write_generated_file(
            os.path.join(self.get_ocp_pybind_dir(), self.__ocp_name, "CMakeLists.txt"),
            "bindings.CMakeLists.txt.in",
            module_targets=(
                "pybind11_add_cgmres_module(ocp)\n"
                "pybind11_add_cgmres_module(zero_horizon_ocp_solver)\n"
                "pybind11_add_cgmres_module(single_shooting_cgmres_solver)\n"
                "pybind11_add_cgmres_module(multiple_shooting_cgmres_solver)"
            ),
            install_destination="cgmres/" + self.__ocp_name,
        )
        write_generated_file(
            os.path.join(self.get_ocp_pybind_dir(), "common", "CMakeLists.txt"),
            "bindings.CMakeLists.txt.in",
            module_targets=(
                "pybind11_add_cgmres_module(solver_settings)\n"
                "pybind11_add_cgmres_module(horizon)\n"
                "pybind11_add_cgmres_module(timer)"
            ),
            install_destination="cgmres/common",
        )
        print("CMakeLists.txt are generated at", self.get_ocp_pybind_dir())

    def git_submodule_update(self) -> None:
        """ Updates git submodules
        """
        print('Update git submodules...')
        subprocess.run(
            ['git', 'submodule', 'update', '--init', '--recursive'], 
            cwd=os.getcwd(), 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=(platform.system()=='Windows')
        )
        print('Successfully updated git submodules\n')

    def build_main(self, generator: str='Auto', vectorize: bool=True,
                   remove_build_dir: bool=False, config: str='Release',
                   parallel: Optional[int]=None, warnings_as_errors: bool=False,
                   sanitizers: bool=False) -> None:
        """ Builds execute file to run numerical simulation. 

            Args: 
                generator: An optional variable for Windows user to choose the
                    generator. If 'MSYS', then 'MSYS Makefiles' is used. If 
                    'MinGW', then 'MinGW Makefiles' is used. The default value 
                    is 'Auto' and the generator is selected automatically. If 
                    sh.exe exists in your PATH, MSYS is choosed, and otherwise 
                    MinGW is used. If different value from 'MSYS' and 'MinGW', 
                    generator is selected automatically.
                vectorize: If True, vectorization ('-march=native' compile option) is enabled.
                    Default is True.
                remove_build_dir: If true, the existing build directory is 
                    removed and if False, the build directory is not removed.
                    Need to be set True is you change CMake configuration, e.g., 
                    if you change the generator. The default value is False.
                warnings_as_errors: Treat compiler warnings as build errors.
                    The default value is False.
                sanitizers: Enable AddressSanitizer and UndefinedBehaviorSanitizer.
                    Requires GCC or Clang. The default value is False.
        """
        if remove_build_dir:
            build_api.remove_build_directory(self.get_ocp_dir())
        build_dir = self.get_ocp_build_dir()
        os.makedirs(build_dir, exist_ok=True)
        if vectorize:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=ON', '-DBUILD_MAIN=ON', '-DBUILD_PYTHON_INTERFACE=OFF']
        else:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=OFF', '-DBUILD_MAIN=ON', '-DBUILD_PYTHON_INTERFACE=OFF']
        if warnings_as_errors:
            build_options.append('-DCGMRES_WARNINGS_AS_ERRORS=ON')
        if sanitizers:
            build_options.append('-DCGMRES_ENABLE_SANITIZERS=ON')
        print('CMake options:', *build_options)
        build_api.build_cpp(
            generator, build_dir, build_options, config=config, parallel=parallel
        )

    def build_python_interface(self, generator: str='Auto', vectorize: bool=True,
                               remove_build_dir: bool=False, config: str='Release',
                               parallel: Optional[int]=None, warnings_as_errors: bool=False,
                               sanitizers: bool=False) -> None:
        """ Builds Python interfaces. 

            Args: 
                generator: An optional variable for Windows user to choose the
                    generator. If 'MSYS', then 'MSYS Makefiles' is used. If 
                    'MinGW', then 'MinGW Makefiles' is used. The default value 
                    is 'Auto' and the generator is selected automatically. If 
                    sh.exe exists in your PATH, MSYS is choosed, and otherwise 
                    MinGW is used. If different value from 'MSYS' and 'MinGW', 
                    generator is selected automatically.
                vectorize: If True, vectorization ('-march=native' compile option) is enabled.
                    Default is True.
                remove_build_dir: If true, the existing build directory is 
                    removed and if False, the build directory is not removed.
                    Need to be set True is you change CMake configuration, e.g., 
                    if you change the generator. The default value is False.
                warnings_as_errors: Treat compiler warnings as build errors.
                    The default value is False.
                sanitizers: Enable AddressSanitizer and UndefinedBehaviorSanitizer.
                    Requires GCC or Clang. The default value is False.
        """
        if remove_build_dir:
            build_api.remove_build_directory(self.get_ocp_dir())
        build_dir = self.get_ocp_build_dir()
        os.makedirs(build_dir, exist_ok=True)
        if vectorize:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=ON', '-DBUILD_MAIN=OFF', '-DBUILD_PYTHON_INTERFACE=ON', '-DPython_EXECUTABLE='+sys.executable]
        else:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=OFF', '-DBUILD_MAIN=OFF', '-DBUILD_PYTHON_INTERFACE=ON', '-DPython_EXECUTABLE='+sys.executable]
        if warnings_as_errors:
            build_options.append('-DCGMRES_WARNINGS_AS_ERRORS=ON')
        if sanitizers:
            build_options.append('-DCGMRES_ENABLE_SANITIZERS=ON')
        print('CMake options:', *build_options)
        build_api.build_cpp(
            generator, build_dir, build_options, config=config, parallel=parallel
        )

    def get_executable_path(self, config: str='Release') -> Path:
        """Return the generated simulation executable for any CMake generator."""
        return build_api.find_executable(self.get_ocp_build_dir(), self.__ocp_name, config)

    def install_python_interface(
        self, install_prefix: Optional[Pathish] = None
    ) -> Path:
        """Installs generated bindings into the running Python environment.

        When ``install_prefix`` is omitted, the active interpreter's
        site-packages directory is used. Activate a virtual environment before
        starting Jupyter to install the bindings into that environment.
        """
        return install_python_interface(
            self.get_ocp_dir(), self.get_ocp_name(), install_prefix
        )

    def run_simulation(self) -> None:
        """ Run numerical simulation. Call after build() succeeded.
        """
        shutil.rmtree(self.get_ocp_log_dir(), ignore_errors=True)
        os.makedirs(self.get_ocp_log_dir(), exist_ok=True)
        executable = self.get_executable_path()
        # Generated simulations write to ../log. Always run from the CMake
        # build root, including when a multi-config generator places the
        # executable in a configuration subdirectory such as Release/.
        build_dir = Path(self.get_ocp_build_dir()).resolve()
        subprocess.run([str(executable)], cwd=build_dir, check=True)
        print('The log files are generated at ', self.get_ocp_log_dir())

def generate_docs() -> None:
    """ Generate docs. Doxygen and webbrowser are required.
    """
    subprocess.run(
        ['doxygen'], 
        cwd=os.path.join(os.getcwd(), 'doc'), 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE, 
        shell=True
    )

def open_docs() -> None:
    import webbrowser
    webbrowser.open('file:///'+str(os.path.join(os.getcwd(), 'doc', 'html', 'annotated.html')))
