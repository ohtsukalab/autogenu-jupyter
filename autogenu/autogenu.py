import subprocess
import platform
from enum import Enum, auto
from collections import namedtuple
import sympy
import os
import sys
import glob, shutil

from . import symutils


class ScalarVariable:
    def __init__(self, symbol: sympy.Symbol, name: str, value=0.0):
        self.symbol = symbol
        self.name = name 
        self.value = value

class ArrayVariable:
    def __init__(self, symbol, name: str, size: int, values=[]):
        assert size > 0
        self.symbol = symbol
        self.name = name 
        self.size = size
        self.values = values

class ControlInputBound:
    def __init__(self, uindex: int, umin, umax, dummy_weight):
        assert uindex >= 0
        assert umin < umax
        assert dummy_weight >= 0.0
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
    def __init__(self, ocp_name: str, nx: int, nu: int):
        assert nx > 0, 'nx must be positive integer!'
        assert nu > 0, 'nu must be positive integer!'
        self.__ocp_name = ocp_name
        self.__nx = nx
        self.__nu = nu
        self.__nc = 0
        self.__nh = 0
        self.__scalar_vars = []
        self.__array_vars = []
        self.__ubounds = []
        self.__symbolic_functions = None
        self.__nlp_type = None
        self.__horizon_params = None
        self.__solver_params = None
        self.__initialization_params = None
        self.__simulation_params = None

    def define_t(self):
        """ Returns symbolic scalar variable 't'.
        """
        return sympy.Symbol('t')

    def define_x(self):
        """ Returns symbolic vector variable 'x' whose size is nx.
        """
        return sympy.symbols('x[0:%d]' %(self.__nx))

    def define_u(self):
        """ Returns symbolic vector variable 'u' whose size is nu.
        """
        return sympy.symbols('u[0:%d]' %(self.__nu))

    def define_scalar_var(self, var_name: str):
        """ Returns symbolic variable whose name is var_name. 

            Args:
                var_name: Name of the scalar variable.
        """
        var_symbol = sympy.Symbol(var_name)
        self.__scalar_vars.append(ScalarVariable(var_symbol, var_name))
        return var_symbol

    def define_scalar_vars(self, *var_name_list):
        """ Returns symbolic variables whose names are given by 
            var_name_list. 

            Args:
                var_name_list: Names of the scalar variables.
        """
        var_symbols = []
        for var_name in var_name_list:
            assert isinstance(var_name, str), 'The input must be list of strings!'
            var_symbol = self.define_scalar_var(var_name)
            var_symbols.append(var_symbol)
        return var_symbols

    def define_array_var(self, var_name: str, size: int):
        """ Returns symbolic vector variable whose names is var_name. 

            Args:
                var_name: Name of the array variable.
                size: Size of the array variable.
        """
        assert size > 0, 'The second argument must be positive integer!'
        array_var = sympy.symbols(var_name+'[0:%d]' %(size))
        self.__array_vars.append(ArrayVariable(array_var, var_name, size, []))
        return array_var

    def set_FB_epsilon(self, FB_epsilon):
        """ Set reguralization term of the semi-smooth Fischer-Burumeister (FB) 
            method. Set the array whose size is dimension of the inequality 
            constraints considered by semi-smooth FB method.

            Args:
                FB epsilon: Array of the reguralization term. 
        """
        for eps in FB_epsilon:
            assert eps >= 0, "FB epsilon must be non-negative!"
        self.__FB_epsilon = FB_epsilon

    def set_scalar_var(self, name: str, value):
        """ Set the value of the scalar variable you defied. 

            Args:
                name: Name of the scalar variable.
                value: Value of the scalar variable.
        """
        for scalar_var in self.__scalar_vars:
            if name == scalar_var.name:
                scalar_var.value = value

    def set_scalar_vars(self, *name_and_value_list):
        """ Set the values of the scalar variables you defied. 

            Args:
                name_and_value_lis: A list composed of the name of 
                the scalar variable and value of the scalar variable.
        """
        for name_and_value in name_and_value_list:
            name = name_and_value[0]
            value = name_and_value[1]
            self.set_scalar_var(name, value)

    def set_array_var(self, name: str, values):
        """ Set the value of the array variable you defied. 

            Args:
                name: Name of the arrray variable.
                values: Values of the arry variable. This size is used as  
                        the size of the array variable.
        """
        for array_var in self.__array_vars:
            if name == array_var.name:
                assert array_var.size == len(values)
                array_var.values = values

    def set_functions(self, f, C, h, L, phi):
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
        assert len(f) > 0 
        assert len(f) == self.__nx, "Dimension of f must be nx!"
        self.__nc = len(C)
        self.__nh = len(h)
        x = sympy.symbols('x[0:%d]' %(self.__nx))
        u = sympy.symbols('u[0:%d]' %(self.__nu+self.__nc+self.__nh))
        lmd = sympy.symbols('lmd[0:%d]' %(self.__nx))
        hamiltonian = L + sum(lmd[i] * f[i] for i in range(self.__nx))
        hamiltonian += sum(u[self.__nu+i] * C[i] for i in range(self.__nc))
        nuc = self.__nu + self.__nc
        hamiltonian += sum(u[nuc+i] * h[i] for i in range(self.__nh))
        hx = symutils.diff_scalar_func(hamiltonian, x)
        hu = symutils.diff_scalar_func(hamiltonian, u)
        fb_eps = sympy.symbols('fb_eps[0:%d]' %(self.__nh))
        for i in range(self.__nh):
            hu[nuc+i] = sympy.sqrt(u[nuc+i]**2 + h[i]**2 + fb_eps[i]) - (u[nuc+i] - h[i])
        phix = symutils.diff_scalar_func(phi, x)
        self.__symbolic_functions = SymbolicFunctions(f, phix, hx, hu)

    def add_control_input_bounds(
        self, uindex: int, umin, umax, dummy_weight
        ):
        """ Adds the bax constraints on the control input that is condensed in 
            linear problem. 

            Args: 
                uindex: The index of the constrianed control input element. 
                umin: The minimum value of the constrianed control input. 
                umax: The minimum value of the constrianed control input. 
                dummy_weight: An weight to stabilize the numerical computation.
        """
        assert uindex >= 0
        assert uindex < self.__nu
        assert umin < umax
        assert dummy_weight >= 0, "dummy_weight must be non-negative!"
        find_same_index = False
        for ub in self.__ubounds:
            if ub.uindex == uindex:
                find_same_index = True
                ub.umin = umin
                ub.umax = umax
                ub.dummy_weight = dummy_weight
        if not find_same_index:
            self.__ubounds.append(ControlInputBound(uindex, umin, umax, dummy_weight))

    def set_nlp_type(self, nlp_type: NLPType):
        """ Sets solver types of the C/GMRES methods. 

            Args: 
                nlp_type: The solver type. Choose from 
                NLPType.SingleShooting and NLPType.MultipleShooting, 
        """
        self.__nlp_type = nlp_type

    def set_horizon_params(self, Tf, alpha=0.0):
        """ Sets parameters of the horizon of NMPC. If alpha > 0, then the 
            length of the horzion at time t is given by Tf * (1-exp(-alpha*t)). 
            If alpha is not positive, the it is given by Tf.

            Args: 
                Tf, alpha: Parameter about the length of the horizon of NMPC.
        """
        assert Tf > 0
        self.__horizon_params = HorizonParams(Tf, alpha)

    def set_solver_params(
            self, sampling_time, N: int, finite_difference_epsilon, zeta, kmax: int
        ):
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
        assert sampling_time > 0
        assert N > 0
        assert finite_difference_epsilon > 0
        assert zeta > 0
        assert kmax > 0
        self.__solver_params = SolverParams(sampling_time, N, finite_difference_epsilon, zeta, kmax)

    def set_initialization_params(
            self, solution_initial_guess, tolerance=1.0e-04, max_iterations: int=100
        ):
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
        assert len(solution_initial_guess) == self.__nu + self.__nc + self.__nh
        assert tolerance >= 0
        assert max_iterations >= 0
        self.__initialization_params = InitializationParams(solution_initial_guess, tolerance, max_iterations)

    def set_simulation_params(
            self, initial_time, initial_state, simulation_length
        ):
        """ Set parameters for numerical simulation. 

            Args: 
                initial_time: The time parameter at the beginning of the 
                    simulation. 
                initial_state: The state of the system at the beginning of the 
                    simulation. 
                simulation_length: The length of the numerical simulation. 
        """
        assert len(initial_state) == self.__nx, "The dimension of initial_state must be nx!"
        assert simulation_length > 0
        self.__simulation_params = SimulationParams(initial_time, initial_state, simulation_length)

    def generate_ocp_definition(self, simplification: bool=False, common_subexpression_elimination: bool=False):
        """ Generates the C++ source file in which the equations to solve the 
            optimal control problem are described. Before call this method, 
            set_functions() must be called.

            Args: 
                simplification: The flag for simplification. If True, the 
                    Symbolic functions are simplified. Default is False.
                common_subexpression_elimination: The flag for common subexpression elimination. If True, 
                    common subexpressions are eliminated. Default is False.
        """
        assert self.__symbolic_functions is not None, \
                "Symbolic functions are not set!. Before call this method, call set_functions()"
        if self.__nh > 0:
            assert len(self.__FB_epsilon) == self.__nh
        self.__make_ocp_dir()
        if simplification:
            symutils.simplify(self.__f)
            symutils.simplify(self.__hx)
            symutils.simplify(self.__hu)
            symutils.simplify(self.__phix)
        f_model_h = open('generated/'+str(self.__ocp_name)+'/ocp.hpp', 'w')
        f_model_h.write('// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). \n')
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
        print('\'ocp.hpp\', the definition of the OCP, is generated at', os.path.abspath('generated/'+str(self.__ocp_name)+'/ocp.hpp'))


    def generate_main(self):
        """ Generates main.cpp that defines NMPC solver, set parameters for the 
            solver, and run numerical simulation. Befire call this method,
            set_nlp_type(), set_horizon_params(), set_solver_params(), 
            set_initialization_params(), and set_simulation_params(),
            must be called!
        """
        assert self.__nlp_type is not None, "Solver type is not set! Before call this method, call set_nlp_type()"
        assert self.__horizon_params is not None, "Horizon params are not set! Before call this method, call set_horizon_params()"
        assert self.__solver_params is not None, "Solver params are not set! Before call this method, call set_solver_params()"
        assert self.__initialization_params is not None, "Initialization params are not set! Before call this method, call set_initialization_params()"
        assert self.__simulation_params is not None, "Simulation params are not set! Before call this method, call set_simulation_params()"
        """ Makes a directory where the C++ source files are generated.
        """
        f_main = open('generated/'+str(self.__ocp_name)+'/main.cpp', 'w')
        f_main.writelines([
""" 
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include "ocp.hpp"
#include "cgmres/zero_horizon_ocp_solver.hpp"
""" 
        ])
        if self.__nlp_type == NLPType.SingleShooting:
            f_main.write('#include "cgmres/single_shooting_cgmres_solver.hpp"')
        elif self.__nlp_type == NLPType.MultipleShooting:
            f_main.write('#include "cgmres/multiple_shooting_cgmres_solver.hpp"')
        else:
            return NotImplementedError()
        f_main.writelines([
""" 

#include "cgmres/simulator/simulator.hpp"
#include <string>

int main() {
""" 
        ])
        f_main.write(
            '  // Define the optimal control problem.\n'
            '  cgmres::OCP_'+str(self.__ocp_name)+' ocp;\n'
            '\n'
        )
        f_main.write(
            '  // Define the horizon.\n'
            '  const double Tf = '+str(self.__horizon_params.Tf)+';\n'
            '  const double alpha = '+str(self.__horizon_params.alpha)+';\n'
            '  cgmres::Horizon horizon(Tf, alpha);\n'
            '\n'
        )
        f_main.write(
            '  // Define the solver settings.\n'
            '  cgmres::SolverSettings settings;\n'
            '  settings.sampling_time = '+str(self.__solver_params.sampling_time)+'; // sampling period \n'
            '  settings.zeta = '+str(self.__solver_params.zeta)+';\n'
            '  settings.finite_difference_epsilon = '+str(self.__solver_params.finite_difference_epsilon)+';\n'
            '  // For initialization.\n'
            '  settings.max_iter = '+str(self.__initialization_params.max_iteraions)+';\n'
            '  settings.opterr_tol = '+str(self.__initialization_params.tolerance)+';\n'
            '\n'
        )
        f_main.write('  // Define the initial time and initial state.\n')
        f_main.write('  const double t0 = '+str(self.__simulation_params.initial_time)+';\n')
        f_main.write(
            '  cgmres::Vector<'+str(len(self.__simulation_params.initial_state))+'> x0;\n'
            +'  x0 << '
        )
        for i in range(len(self.__simulation_params.initial_state)-1):
            f_main.write(str(self.__simulation_params.initial_state[i])+', ')
        f_main.write(str(self.__simulation_params.initial_state[-1])+';\n')
        f_main.write('\n')
        # initial guess for the initialization of the solution
        f_main.write('  // Initialize the solution of the C/GMRES method.\n')
        f_main.write('  constexpr int kmax_init = '+str(min(self.__solver_params.kmax, self.__nu+self.__nc+self.__nh))+';\n')
        f_main.write(
            '  cgmres::ZeroHorizonOCPSolver<cgmres::OCP_'+self.__ocp_name+', kmax_init> '
            +'initializer(ocp, settings);\n'
        )
        f_main.write('  cgmres::Vector<'+str(self.__nu+self.__nc+self.__nh)+'> uc0;\n'+'  uc0 << ')
        for i in range(len(self.__initialization_params.solution_initial_guess)-1):
            f_main.write(str(self.__initialization_params.solution_initial_guess[i])+', ')
        f_main.write(str(self.__initialization_params.solution_initial_guess[-1])+';\n')
        f_main.write('  initializer.set_uc(uc0);\n')
        f_main.write('  initializer.solve(t0, x0);\n')
        f_main.write('\n')

        f_main.write('  // Define the C/GMRES solver.\n')
        f_main.write('  constexpr int N = '+str(self.__solver_params.N)+';\n')
        f_main.write('  constexpr int kmax = '+str(min(self.__solver_params.kmax, self.__solver_params.N*(self.__nu+self.__nc+self.__nh)))+';\n')
        if self.__nlp_type == NLPType.SingleShooting:
            f_main.write(
                '  cgmres::SingleShootingCGMRESSolver<cgmres::OCP_'+self.__ocp_name+', N, kmax> mpc(ocp, horizon, settings);\n'
                '  mpc.set_uc(initializer.ucopt());\n'
                '  mpc.init_dummy_mu();\n'
            )
        elif self.__nlp_type == NLPType.MultipleShooting:
            f_main.write(
                '  cgmres::MultipleShootingCGMRESSolver<cgmres::OCP_'+self.__ocp_name+', N, kmax> mpc(ocp, horizon, settings);\n'
                '  mpc.set_uc(initializer.ucopt());\n'
                '  mpc.init_x_lmd(t0, x0);\n'
                '  mpc.init_dummy_mu();\n'
            )
        else:
            return NotImplementedError()
        f_main.write('\n\n')
        f_main.write(
            '  // Perform a numerical simulation.\n'
            '  const double tsim = '+str(self.__simulation_params.simulation_length)+';\n'
            '  const double sampling_time = settings.sampling_time;\n'
            '  const std::string save_dir_name("../simulation_result");\n'
            '  cgmres::simulation(ocp, mpc, x0, t0, tsim, sampling_time, ' 
            +"save_dir_name"+', "' +self.__ocp_name +'");\n\n'
            '  return 0;\n'
            '}\n'
        )
        f_main.close()
        print('\'main.cpp\', the closed-loop simulation code, is generated at', os.path.abspath('generated/'+str(self.__ocp_name)+'/main.cpp'))

    def generate_python_bindings(self):
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/'+str(self.__ocp_name)+'/ocp.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
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
        if (v.size() != self.umin.size()) { 
          throw std::invalid_argument("[OCP]: 'umin.size()' must be "+std::to_string(self.umin.size()));
        } Map<VectorX>(self.umin.data(), self.umin.size()) = v; })
    .def_property("umax", 
      [](const OCP& self) { return Map<const VectorX>(self.umax.data(), self.umax.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != self.umax.size()) { 
          throw std::invalid_argument("[OCP]: 'umax.size()' must be "+std::to_string(self.umax.size()));
        } Map<VectorX>(self.umax.data(), self.umax.size()) = v; })
    .def_property("dummy_weight", 
      [](const OCP& self) { return Map<const VectorX>(self.dummy_weight.data(), self.dummy_weight.size()); },
      [](OCP& self, const VectorX& v) { 
        if (v.size() != self.dummy_weight.size()) {
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
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/'+str(self.__ocp_name)+'/zero_horizon_ocp_solver.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/zero_horizon_ocp_solver.hpp"
#include "cgmres/python/zero_horizon_ocp_solver.hpp"
#include "ocp.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

""" 
        ])
        f_pybind11.write('constexpr int kmax_init = '+str(min(self.__solver_params.kmax, self.__nc+self.__nu+self.__nh))+';\n')
        f_pybind11.write('DEFINE_PYBIND11_MODULE_ZERO_HORIZON_OCP_SOLVER(OCP_'+str(self.__ocp_name)+', kmax_init)\n')
        f_pybind11.writelines([
"""

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/'+str(self.__ocp_name)+'/single_shooting_cgmres_solver.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/single_shooting_cgmres_solver.hpp"
#include "cgmres/python/single_shooting_cgmres_solver.hpp"
#include "ocp.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

""" 
        ])
        f_pybind11.write('constexpr int N = '+str(self.__solver_params.N)+';\n')
        f_pybind11.write('constexpr int kmax = '+str(min(self.__solver_params.kmax, self.__solver_params.N*(self.__nc+self.__nu+self.__nh)))+';\n')
        f_pybind11.write('DEFINE_PYBIND11_MODULE_SINGLE_SHOOTING_CGMRES_SOLVER(OCP_'+str(self.__ocp_name)+', N, kmax)\n')
        f_pybind11.writelines([
"""

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/'+str(self.__ocp_name)+'/multiple_shooting_cgmres_solver.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/multiple_shooting_cgmres_solver.hpp"
#include "cgmres/python/multiple_shooting_cgmres_solver.hpp"
#include "ocp.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

""" 
        ])
        f_pybind11.write('constexpr int N = '+str(self.__solver_params.N)+';\n')
        f_pybind11.write('constexpr int kmax = '+str(min(self.__solver_params.kmax, self.__solver_params.N*(self.__nc+self.__nu+self.__nh)))+';\n')
        f_pybind11.write('DEFINE_PYBIND11_MODULE_MULTIPLE_SHOOTING_CGMRES_SOLVER(OCP_'+str(self.__ocp_name)+', N, kmax)\n')
        f_pybind11.writelines([
"""

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/common/horizon.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/types.hpp"
#include "cgmres/horizon.hpp"
#include "cgmres/python/horizon.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

DEFINE_PYBIND11_MODULE_HORIZON()

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/common/solver_settings.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/solver_settings.hpp"
#include "cgmres/python/solver_settings.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

DEFINE_PYBIND11_MODULE_SOLVER_SETTINGS()

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/common/timer.cpp', 'w')
        f_pybind11.writelines([
"""
// This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
// The autogenu-jupyter copyright holders make no ownership claim of its contents. 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "cgmres/timer.hpp"
#include "cgmres/python/timer.hpp"

#include <iostream>
#include <stdexcept>

namespace cgmres {
namespace python {

namespace py = pybind11;

DEFINE_PYBIND11_MODULE_TIMER()

} // namespace python
} // namespace cgmres
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/'+str(self.__ocp_name)+'/__init__.py', 'w')
        f_pybind11.writelines([
"""
# This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
# The autogenu-jupyter copyright holders make no ownership claim of its contents. 

from .ocp import *
from .zero_horizon_ocp_solver import *
from .single_shooting_cgmres_solver import *
from .multiple_shooting_cgmres_solver import *
""" 
        ])
        f_pybind11.close()
        f_pybind11 = open('generated/'+str(self.__ocp_name)+'/python/common/__init__.py', 'w')
        f_pybind11.writelines([
"""
# This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
# The autogenu-jupyter copyright holders make no ownership claim of its contents. 

from .horizon import *
from .solver_settings import *
from .timer import *
""" 
        ])
        f_pybind11.close()
        print('pybind11 source codes are generated at', os.path.abspath('generated/'+str(self.__ocp_name)+'/python'))


    def generate_cmake(self):
        """ Generates CMakeLists.txt in a directory where your .ipynb files 
            locates.
        """
        f_cmake = open('generated/'+str(self.__ocp_name)+'/CMakeLists.txt', 'w')
        f_cmake.writelines([
"""
# This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
# The autogenu-jupyter copyright holders make no ownership claim of its contents. 

cmake_minimum_required(VERSION 3.1)
""" 
        ])
        f_cmake.write('project('+str(self.__ocp_name)+' CXX)')
        f_cmake.writelines([
"""

set(CMAKE_CXX_STANDARD 17)

option(VECTORIZE "Enable -march=native" OFF)
option(BUILD_MAIN "Build C++ simulation" ON)
option(BUILD_PYTHON_INTERFACE "Build Python interface" OFF)

set(CGMRES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../../include)

if (BUILD_MAIN)
  add_executable(
      ${PROJECT_NAME}
      main.cpp
  )
  target_include_directories(
      ${PROJECT_NAME}
      PRIVATE
      ${CGMRES_INCLUDE_DIR}
  )
  if (VECTORIZE)
    target_compile_options(
      ${PROJECT_NAME}
      PRIVATE
      -march=native
    )
  endif()
endif()

if (BUILD_PYTHON_INTERFACE)
    add_subdirectory(python/common)
    add_subdirectory(python/${PROJECT_NAME})
endif()
"""
            ])
        f_cmake.close()
        f_cmake_python = open('generated/'+self.__ocp_name+'/python/'+self.__ocp_name+'/CMakeLists.txt', 'w')
        f_cmake_python.writelines([
"""
# This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
# The autogenu-jupyter copyright holders make no ownership claim of its contents. 

macro(pybind11_add_cgmres_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_include_directories(
    ${MODULE} 
    PRIVATE
    ${CGMRES_INCLUDE_DIR}
    ${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/eigen
    ${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/pybind11
    ${PROJECT_SOURCE_DIR}
  )
    if (VECTORIZE)
    target_compile_options(
        ${MODULE}
        PRIVATE
        -march=native
    )
    endif()
endmacro()

add_subdirectory(${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/pybind11)
pybind11_add_cgmres_module(ocp)
pybind11_add_cgmres_module(zero_horizon_ocp_solver)
pybind11_add_cgmres_module(single_shooting_cgmres_solver)
pybind11_add_cgmres_module(multiple_shooting_cgmres_solver)

set(CGMRES_PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
"""
            ])
        if platform.system() == 'Windows':
            f_cmake_python.write(
                'set(CGMRES_PYTHON_BINDINGS_LIBDIR $ENV{HOMEPATH}/.local/lib/python${CGMRES_PYTHON_VERSION}/site-packages/cgmres/'+self.__ocp_name+')')
        else:
            f_cmake_python.write(
                'set(CGMRES_PYTHON_BINDINGS_LIBDIR $ENV{HOME}/.local/lib/python${CGMRES_PYTHON_VERSION}/site-packages/cgmres/'+self.__ocp_name+')')
        f_cmake_python.writelines([
"""
file(GLOB PYTHON_BINDINGS_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/*.cpython*)
file(GLOB PYTHON_FILES_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
install(
  FILES ${PYTHON_FILES_${CURRENT_MODULE_DIR}} ${PYTHON_BINDINGS_${CURRENT_MODULE_DIR}} 
  DESTINATION ${CGMRES_PYTHON_BINDINGS_LIBDIR}/${CURRENT_MODULE_DIR}
)
"""
            ])
        f_cmake_python.close()
        f_cmake_python = open('generated/'+self.__ocp_name+'/python/common/CMakeLists.txt', 'w')
        f_cmake_python.writelines([
"""
# This file was automatically generated by autogenu-jupyter (https://github.com/mayataka/autogenu-jupyter). 
# The autogenu-jupyter copyright holders make no ownership claim of its contents. 

macro(pybind11_add_cgmres_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_include_directories(
    ${MODULE} 
    PRIVATE
    ${CGMRES_INCLUDE_DIR}
    ${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/eigen
    ${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/pybind11
    ${PROJECT_SOURCE_DIR}
  )
    if (VECTORIZE)
    target_compile_options(
        ${MODULE}
        PRIVATE
        -march=native
    )
    endif()
endmacro()

add_subdirectory(${CGMRES_INCLUDE_DIR}/cgmres/thirdparty/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/pybind11)
pybind11_add_cgmres_module(solver_settings)
pybind11_add_cgmres_module(horizon)
pybind11_add_cgmres_module(timer)

set(CGMRES_PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
"""
            ])
        if platform.system() == 'Windows':
            f_cmake_python.write(
                'set(CGMRES_PYTHON_BINDINGS_LIBDIR $ENV{HOMEPATH}/.local/lib/python${CGMRES_PYTHON_VERSION}/site-packages/cgmres/common)')
        else:
            f_cmake_python.write(
                'set(CGMRES_PYTHON_BINDINGS_LIBDIR $ENV{HOME}/.local/lib/python${CGMRES_PYTHON_VERSION}/site-packages/cgmres/common)')
        f_cmake_python.writelines([
"""
file(GLOB PYTHON_BINDINGS_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/*.cpython*)
file(GLOB PYTHON_FILES_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
install(
  FILES ${PYTHON_FILES_${CURRENT_MODULE_DIR}} ${PYTHON_BINDINGS_${CURRENT_MODULE_DIR}} 
  DESTINATION ${CGMRES_PYTHON_BINDINGS_LIBDIR}/${CURRENT_MODULE_DIR}
)
"""
            ])
        f_cmake_python.close()
        print('CMakeLists.txt are generated at', os.path.abspath('generated/'+str(self.__ocp_name)))

    def git_submodule_update(self):
        """ Updates git submodules
        """
        print('Update git submodules...')
        if platform.system() == 'Windows':
            subprocess.run(
                ['git', 'submodule', 'update', '--init', '--recursive'], 
                cwd='.', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
        else:
            subprocess.run(
                ['git', 'submodule', 'update', '--init', '--recursive'], 
                cwd='.', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
        print('Successfully updated git submodules\n')

    def build_main(self, generator: str='Auto', vectorize: bool=True, 
                   remove_build_dir: bool=False):
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
        """
        if remove_build_dir:
            self.__remove_build_dir()
        if vectorize:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=ON', '-DBUILD_MAIN=ON', '-DBUILD_PYTHON_INTERFACE=OFF']
        else:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=OFF', '-DBUILD_MAIN=ON', '-DBUILD_PYTHON_INTERFACE=OFF']
        print('CMake options:', *build_options)
        if platform.system() == 'Windows':
            subprocess.run(
                ['mkdir', 'build'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            if generator == 'MSYS':
                proc = subprocess.Popen(
                    ['cmake', '..', '-G', 'MSYS Makefiles', *build_options], 
                    cwd='generated/'+self.__ocp_name+'/build', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    shell=True
                )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            elif generator == 'MinGW':
                proc = subprocess.Popen(
                    ['cmake', '..', '-G', 'MinGW Makefiles', *build_options], 
                    cwd='generated/'+self.__ocp_name+'/build', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    shell=True
                )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            else:
                proc = subprocess.Popen(
                    ['where', 'sh.exe'], 
                    cwd='C:', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, 
                    shell=True
                )
                if proc.stderr.readline() == b'':
                    proc = subprocess.Popen(
                        ['cmake', '..', '-G', 'MSYS Makefiles', *build_options], 
                        cwd='generated/'+self.__ocp_name+'/build', 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT, 
                        shell=True
                    )
                else:
                    proc = subprocess.Popen(
                        ['cmake', '..', '-G', 'MinGW Makefiles', *build_options], 
                        cwd='generated/'+self.__ocp_name+'/build', 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT, 
                        shell=True
                    )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            proc = subprocess.Popen(
                ['cmake', '--build', '.'], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
            for line in iter(proc.stdout.readline,b''):
                print(line.rstrip().decode("utf8"))
            print('\n')
            
        else:
            subprocess.run(
                ['mkdir', 'build'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            proc = subprocess.Popen(
                ['cmake', '..', *build_options], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
            print('\n')
            proc = subprocess.Popen(
                ['cmake', '--build', '.'], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout = subprocess.PIPE, 
                stderr = subprocess.STDOUT
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
            print('\n')

    def build_python_interface(self, generator: str='Auto', vectorize: bool=True, 
                               remove_build_dir: bool=False):
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
        """
        if remove_build_dir:
            self.__remove_build_dir()
        if vectorize:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=ON', '-DBUILD_MAIN=OFF', '-DBUILD_PYTHON_INTERFACE=ON']
        else:
            build_options = ['-DCMAKE_BUILD_TYPE=Release', '-DVECTORIZE=OFF', '-DBUILD_MAIN=OFF', '-DBUILD_PYTHON_INTERFACE=ON']
        print('CMake options:', *build_options)
        if platform.system() == 'Windows':
            subprocess.run(
                ['mkdir', 'build'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            if generator == 'MSYS':
                proc = subprocess.Popen(
                    ['cmake', '..', '-G', 'MSYS Makefiles', *build_options], 
                    cwd='generated/'+self.__ocp_name+'/build', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    shell=True
                )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            elif generator == 'MinGW':
                proc = subprocess.Popen(
                    ['cmake', '..', '-G', 'MinGW Makefiles', *build_options], 
                    cwd='generated/'+self.__ocp_name+'/build', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    shell=True
                )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            else:
                proc = subprocess.Popen(
                    ['where', 'sh.exe'], 
                    cwd='C:', 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, 
                    shell=True
                )
                if proc.stderr.readline() == b'':
                    proc = subprocess.Popen(
                        ['cmake', '..', '-G', 'MSYS Makefiles', *build_options], 
                        cwd='generated/'+self.__ocp_name+'/build', 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT, 
                        shell=True
                    )
                else:
                    proc = subprocess.Popen(
                        ['cmake', '..', '-G', 'MinGW Makefiles', *build_options], 
                        cwd='generated/'+self.__ocp_name+'/build', 
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.STDOUT, 
                        shell=True
                    )
                for line in iter(proc.stdout.readline, b''):
                    print(line.rstrip().decode("utf8"))
                print('\n')
            proc = subprocess.Popen(
                ['cmake', '--build', '.', '--parallel', '8'], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
            for line in iter(proc.stdout.readline,b''):
                print(line.rstrip().decode("utf8"))
            print('\n')
        else:
            subprocess.run(
                ['mkdir', 'build'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            proc = subprocess.Popen(
                ['cmake', '..', *build_options], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
            print('\n')
            proc = subprocess.Popen(
                ['cmake', '--build', '.', '-j8'], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout = subprocess.PIPE, 
                stderr = subprocess.STDOUT
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
            print('\n')

    def install_python_interface(self, install_prefix=None):
        if install_prefix is None:
            python_version = 'python' + str(sys.version_info.major) + '.' + str(sys.version_info.minor)
            if platform.system() == 'Windows':
                install_prefix = os.path.join(os.path.abspath(os.environ['HOMEPATH']), '.local/lib', python_version, 'site-packages')
            else:
                install_prefix = os.path.join(os.path.abspath(os.environ['HOME']), '.local/lib', python_version, 'site-packages')
            install_destination = os.path.join(os.path.abspath(install_prefix), 'cgmres')
        else:
            install_destination = os.path.join(os.path.abspath(install_prefix), 'cgmres')
        pybind11_sharedlibs = glob.glob('generated/'+self.__ocp_name+'/build/python/'+self.__ocp_name+'/*.so') \
                                + glob.glob('generated/'+self.__ocp_name+'/build/python/'+self.__ocp_name+'/*.pyd')
        pybind11_sharedlibs_common = glob.glob('generated/'+self.__ocp_name+'/build/python/common/*.so') \
                                        + glob.glob('generated/'+self.__ocp_name+'/build/python/common/*.pyd')
        if platform.system() == 'Windows':
            os.makedirs(os.path.join(install_destination, self.__ocp_name), exist_ok=True)
            os.makedirs(os.path.join(install_destination, 'common'), exist_ok=True)
        else:
            subprocess.run(
                ['mkdir', '-p', str(os.path.join(install_destination, self.__ocp_name)), str(os.path.join(install_destination, 'common'))], 
                cwd='.',
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
            )
        for e in pybind11_sharedlibs:
            shutil.copy(e, str(os.path.join(install_destination, self.__ocp_name)))
        for e in pybind11_sharedlibs_common:
            shutil.copy(e, str(os.path.join(install_destination, 'common')))
        python_files = glob.glob('generated/'+self.__ocp_name+'/python/'+self.__ocp_name+'/*.py')
        python_files_common = glob.glob('generated/'+self.__ocp_name+'/python/common/*.py')
        for e in python_files:
            shutil.copy(e, str(os.path.join(install_destination, self.__ocp_name)))
        for e in python_files_common:
            shutil.copy(e, str(os.path.join(install_destination, 'common')))
        print('Python interfaces have been installed at ' + str(install_prefix))
        print('To use Python interfaces, run \n')
        print('    export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + '\n')
        print('in the terminal to recognize the PYTHONPATH temporary.')
        print('Or set the PATH in Ubuntu as\n')
        print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~./bashrc\n')
        print('or in Mac OSX as\n')
        print('    echo export PYTHONPATH=$PYTHONPATH:' + str(install_prefix) + ' >> ~./zshrc\n')

    def run_simulation(self):
        """ Run numerical simulation. Call after build() succeeded.
        """
        if platform.system() == 'Windows':
            subprocess.run(
                ['rmdir', '/q', '/s', 'simulation_result'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', 'simulation_result'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            proc = subprocess.Popen(
                [self.__ocp_name+'.exe'], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                shell=True
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))
        else:
            subprocess.run(
                ['rm', '-rf', 'simulation_result'], 
                cwd='generated/'+self.__ocp_name, 
                stdout = subprocess.PIPE, 
                stderr = subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', 'simulation_result'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            proc = subprocess.Popen(
                ['./'+self.__ocp_name], 
                cwd='generated/'+self.__ocp_name+'/build', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT
            )
            for line in iter(proc.stdout.readline, b''):
                print(line.rstrip().decode("utf8"))

    def __make_ocp_dir(self):
        """ Makes a directory where the C source files of OCP formulations are 
            generated.
        """
        if platform.system() == 'Windows':
            subprocess.run(
                ['mkdir', 'generated'], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', self.__ocp_name], 
                cwd='generated', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', 'python'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', self.__ocp_name], 
                cwd='generated/'+self.__ocp_name+'/python', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
            subprocess.run(
                ['mkdir', 'common'], 
                cwd='generated/'+self.__ocp_name+'/python', 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
        else:
            subprocess.run(
                ['mkdir', 'generated'], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            subprocess.run(
                ['mkdir', self.__ocp_name], 
                cwd='generated',
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            subprocess.run(
                ['mkdir', '-p', 'python/'+self.__ocp_name], 
                cwd='generated/'+self.__ocp_name,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )
            subprocess.run(
                ['mkdir', '-p', 'python/common'], 
                cwd='generated/'+self.__ocp_name,
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )

    def __remove_build_dir(self):
        """ Removes a build directory. This function is mainly for Windows 
            users with MSYS.
        """
        if platform.system() == 'Windows':
            subprocess.run(
                ['rmdir', '/q', '/s', 'build'], 
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True
            )
        else:
            subprocess.run(
                ['rm', '-r', 'build'],
                cwd='generated/'+self.__ocp_name, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE
            )