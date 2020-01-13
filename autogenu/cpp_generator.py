import linecache
import subprocess
import platform

import sympy

from autogenu import solver_parameters as slpr
from autogenu import initialization_parameters as inipr 
from autogenu import simulation_parameters as simpr
from autogenu import cpp_executor as cppexe


def make_model_dir(model_name):
    """ Makes a directory where the C++ source files are generated.

        Args: 
            model_name: A string representing the name of the simulation model.
    """
    if platform.system() == 'Windows':
        subprocess.run(
            ['mkdir', 'models'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
        subprocess.run(
            ['mkdir', model_name], 
            cwd='models', 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            shell=True
        )
    else:
        subprocess.run(
            ['mkdir', 'models'], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
        subprocess.run(
            ['mkdir', model_name], 
            cwd='models',
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )


def generate_cpp(fxu, phix, hx, hu, model_name, cse_flag=False):
    """ Generates the C++ source file in which the equations to solve the 
        optimal control problem are described.

        Args: 
            fxu: The state equation.
            phix: The partial derivative of the terminal cost with respect to 
                the state.
            hx: The partial derivative of the hamiltonina with respect to the 
                state.
            hu: The partial derivative of the hamiltonina with respect to the 
                control input.
            model_name: A string representing the name of the simulation model.
            cse_flag: The flag for common subexpression elimination. If True, 
                common subexpressions in fxu, phix, hx, and hu are eliminated. 
                Default is False.
    """
    f_model_cpp = open('models/'+str(model_name)+'/nmpc_model.cpp', 'w')
    cpp_template = 'autogenu/.cpp_templates/nmpc_model.cpp'
    if cse_flag:
        fxu = sympy.cse(fxu)
        phix = sympy.cse(phix)
        hx = sympy.cse(hx)
        hu = sympy.cse(hu)
    f_model_cpp.writelines(
        [linecache.getline(cpp_template, i) for i in range(0, 12)]
    )
    if cse_flag:
        for i in range(len(fxu[0])):
            cse_exp, cse_rhs = fxu[0][i]
            f_model_cpp.write(
                '  double '+sympy.ccode(cse_exp)+' = '+sympy.ccode(cse_rhs)+';\n'
            )
        f_model_cpp.writelines(
            ['  f[%d] = '%i+sympy.ccode(fxu[1][i])+';\n' for i in range(len(fxu[1]))]
        )
    else:
        f_model_cpp.writelines(
            ['  f[%d] = '%i+sympy.ccode(fxu[i])+';\n' for i in range(len(fxu))]
        )
    f_model_cpp.writelines(
        [linecache.getline(cpp_template, i) for i in range(13, 21)]
    )
    if cse_flag:
        for i in range(len(phix[0])):
            cse_exp, cse_rhs = phix[0][i]
            f_model_cpp.write(
                '  double '+sympy.ccode(cse_exp)+' = '+sympy.ccode(cse_rhs)+';\n'
            )
        f_model_cpp.writelines(
            ['  phix[%d] = '%i+sympy.ccode(phix[1][i])+';\n' for i in range(len(phix[1]))]
        )
    else:
        f_model_cpp.writelines(
            ['  phix[%d] = '%i+sympy.ccode(phix[i])+';\n' for i in range(len(phix))]
        )
    f_model_cpp.writelines(
        [linecache.getline(cpp_template, i) for i in range(22, 32)]
    )
    if cse_flag:
        for i in range(len(hx[0])):
            cse_exp, cse_rhs = hx[0][i]
            f_model_cpp.write(
                '  double '+sympy.ccode(cse_exp)+' = '+sympy.ccode(cse_rhs)+';\n'
            )
        f_model_cpp.writelines(
            ['  hx[%d] = '%i+sympy.ccode(hx[1][i])+';\n' for i in range(len(hx[1]))]
        )
    else:
        f_model_cpp.writelines(
            ['  hx[%d] = '%i+sympy.ccode(hx[i])+';\n' for i in range(len(hx))]
        )
    f_model_cpp.writelines(
        [linecache.getline(cpp_template, i) for i in range(33, 43)]
    )
    if cse_flag:
        for i in range(len(hu[0])):
            cse_exp, cse_rhs = hu[0][i]
            f_model_cpp.write(
                '  double '+sympy.ccode(cse_exp)+' = '+sympy.ccode(cse_rhs)+';\n'
            )
        f_model_cpp.writelines(
            ['  hu[%d] = '%i+sympy.ccode(hu[1][i])+';\n' for i in range(len(hu[1]))]
        )
    else:
        f_model_cpp.writelines(
            ['  hu[%d] = '%i+sympy.ccode(hu[i])+';\n' for i in range(len(hu))]
        )
    f_model_cpp.writelines(
        [linecache.getline(cpp_template, i) for i in range(44, 62)]
    )
    f_model_cpp.close()


def generate_hpp(
        dimx, dimu, dimc, scalar_parameters, array_parameters, model_name
    ):
    """ Generates the C++ source header in which the equations to solve the 
        optimal control problem are decleared and parameters are defined.

        Args: 
            dimx: Dimension of the state.
            dimu: Dimension of the control input.
            dimc: Dimenstion of the equality constraints.
            scalar_parameters: Parameters defined by users and used in the 
                state equation, constraints, and cost function.
            array_parameters: Parameters defined by users and used in the cost 
                function.
            model_name: A string representing the name of the simulation model.
    """
    f_model_hpp = open('models/'+str(model_name)+'/nmpc_model.hpp', 'w')
    hpp_template = 'autogenu/.cpp_templates/nmpc_model.hpp'
    f_model_hpp.writelines(
        [linecache.getline(hpp_template, i) for i in range(0, 18)]
    )
    f_model_hpp.write(
        '  static constexpr int dim_state_ = %d;\n' %dimx
    )
    f_model_hpp.write(
        '  static constexpr int dim_control_input_ = %d;\n' %dimu
    )
    f_model_hpp.write(
        '  static constexpr int dim_constraints_ = %d;\n' %dimc
    )
    f_model_hpp.write('\n')
    f_model_hpp.writelines([
            '  static constexpr double '+str(scalar_parameters[i][0])
            +' = '+str(scalar_parameters[i][1])+';\n' 
            for i in range(len(scalar_parameters))
    ])
    f_model_hpp.write('\n')
    f_model_hpp.writelines([
        '  double '+str(array_parameters[i][0])+'['
        +str(array_parameters[i][1])+']'+' = '
        +str(array_parameters[i][2])+';\n' 
        for i in range(len(array_parameters))
    ])
    f_model_hpp.writelines(
        [linecache.getline(hpp_template, i) for i in range(18, 68)]
    )
    f_model_hpp.close()


def generate_main(
        solver_index, solver_parameters, initialization_parameters, 
        simulation_parameters, model_name, saturation_list=None
    ):
    """ Makes a directory where the C++ source files are generated.

        Args: 
            solver_index: Representing the choice of the solvers.
            solver_parameters: Parameters for NMPC solver.
            initialization_parameters: Parameters for initialization of the 
                solution of NMPC. 
            simulation_parameters: Parameters about the numerical simulation.
            model_name: A string representing the name of the simulation model.
            saturation_list: Optional parameters with respect to the saturation
                of the control input and just valid when solver_index is 3.
    """
    # include header
    f_main = open('models/'+str(model_name)+'/main.cpp', 'w')
    f_main.write('#include "nmpc_model.hpp"\n')
    if solver_index == 1:
        f_main.write(
            '#include "continuation_gmres.hpp"\n'
        )
    elif solver_index == 2:
        f_main.write(
            '#include "multiple_shooting_cgmres.hpp"\n'
        )
    else:
        f_main.write(
            '#include "input_saturation_set.hpp"\n'
            '#include "ms_cgmres_with_input_saturation.hpp"\n'
        )
    f_main.write(
        '#include "cgmres_simulator.hpp"\n'
    )
    f_main.write('#include <string>\n')
    if platform.system() == 'Windows':
        f_main.write(
            '#include <direct.h>\n'
            '\n'
        )
    else:
        f_main.write(
            '#include <sys/stat.h>\n'
            '\n'
        )
    f_main.write(
        '\n'
        'int main()\n'
        '{\n'
        '    // Define the model in NMPC.\n'
        '    cgmres::NMPCModel nmpc_model;\n'
        '\n'
    )
    # define solver
    f_main.write('    // Define the solver.\n')
    if solver_index == 1:
        f_main.write(
            '    cgmres::ContinuationGMRES nmpc_solver('
            +str(solver_parameters.T_f)+', '
            +str(solver_parameters.alpha)+', '
            +str(solver_parameters.N)+', '
            +str(solver_parameters.finite_difference_increment)+', '
            +str(solver_parameters.zeta)+', '
            +str(solver_parameters.kmax)+');\n'
        )
    elif solver_index == 2:
        f_main.write(
            '    cgmres::MultipleShootingCGMRES nmpc_solver('
            +str(solver_parameters.T_f)+', '
            +str(solver_parameters.alpha)+', '
            +str(solver_parameters.N)+', '
            +str(solver_parameters.finite_difference_increment)+', '
            +str(solver_parameters.zeta)+', '
            +str(solver_parameters.kmax)+');\n'
        )
    else:
        f_main.write(
            '    cgmres::InputSaturationSet '
            'input_saturation_set;\n'
        )
        for i in range(len(saturation_list)):
            f_main.write(
                '    input_saturation_set.'
                'appendInputSaturation('
                +str(saturation_list[i][0])+', '
                +str(saturation_list[i][1])+', '
                +str(saturation_list[i][2])+', '
                +str(saturation_list[i][3])+', '
                +str(saturation_list[i][4])+');\n'
            )
        f_main.write(
            '    cgmres::MSCGMRESWithInputSaturation '
            'nmpc_solver(input_saturation_set, '
            +str(solver_parameters.T_f)+', '
            +str(solver_parameters.alpha)+', '
            +str(solver_parameters.N)+', '
            +str(solver_parameters.finite_difference_increment)+', '
            +str(solver_parameters.zeta)+', '
            +str(solver_parameters.kmax)+');\n'
        )
    f_main.write('\n')
    f_main.write('\n')
    # initial state
    f_main.write('    // Set the initial state.\n')
    f_main.write(
        '    double initial_state['
        +str(len(simulation_parameters.initial_state))+
        '] = {'
    )
    for i in range(len(simulation_parameters.initial_state)-1):
        f_main.write(str(simulation_parameters.initial_state[i])+', ')
    f_main.write(str(simulation_parameters.initial_state[-1])+'};\n')
    f_main.write('\n')
    # initial guess for the initialization of the solution
    f_main.write('    // Set the initial guess of the solution.\n')
    f_main.write(
        '    double initial_guess_solution['
        +str(len(initialization_parameters.initial_guess_solution))
        +'] = {'
    )
    for i in range(len(initialization_parameters.initial_guess_solution)-1):
        f_main.write(
            str(initialization_parameters.initial_guess_solution[i])+', '
        )
    f_main.write(
        str(initialization_parameters.initial_guess_solution[-1])+'};\n'
    )
    f_main.write('\n')
    f_main.write('\n')
    # initialization of the solution of the C/GMRES method
    f_main.write('    // Initialize the solution of the C/GMRES method.\n')
    f_main.write(
        '    nmpc_solver.setParametersForInitialization('
        +'initial_guess_solution, '
        +str(initialization_parameters.newton_residual_torelance)+', '
        +str(initialization_parameters.max_newton_iteration)+');\n'
    )
    f_main.write('\n')
    if (solver_index == 3 
        and initialization_parameters.initial_Lagrange_multiplier is not None):
        f_main.write(
            '    // Set the initial guess of the lagrange multiplier '
            'for the condensed constraints with respect to the saturation '
            'on the function of the control input .\n'
            )
        f_main.write(
            '    double initial_guess_lagrange_multiplier['
            +str(len(initialization_parameters.initial_Lagrange_multiplier))
            +'] = {'
        )
        for i in range(len(initialization_parameters.initial_Lagrange_multiplier)-1):
            f_main.write(
                str(initialization_parameters.initial_Lagrange_multiplier[i])+', '
            )
        f_main.write(
            str(initialization_parameters.initial_Lagrange_multiplier[-1])+'};\n'
        )
        f_main.write(
            '\n'+'    nmpc_solver.setInitialInputSaturationMultiplier'
            +'(initial_guess_lagrange_multiplier);'
            +'\n'
        )
    f_main.write('\n')
    f_main.write('\n')
    if platform.system() == 'Windows':
        f_main.write(
            '    // Makes a directory for saving simulation results.\n'
            '    std::string save_dir_name("../simulation_result");\n'
            '    int mkdir_err = mkdir(save_dir_name.c_str());\n'
            '\n'
        )
    else:
        f_main.write(
            '    // Makes a directory for saving simulation results.\n'
            '    std::string save_dir_name("../simulation_result");\n'
            '    int mkdir_err = mkdir(save_dir_name.c_str(), 0755);'
            '\n'
        )
    f_main.write('\n')
    f_main.write(
        '    // Perform a numerical simulation.\n'
        '    cgmres::simulation(nmpc_solver, initial_state, '
        +str(simulation_parameters.initial_time)+', '
        +str(simulation_parameters.simulation_time)+', '
        +str(simulation_parameters.sampling_time)+', '
        +"save_dir_name"+', "'+model_name+'");\n'
        '\n'
        '    return 0;\n'
        '}\n'
    )
    f_main.close()


def generate_cmake(solver_index, model_name):
    """ Generates CMakeLists.txt in a directory where your .ipynb files rocates.

        Args: 
            solver_index: Index representing the choice of the solvers.
            model_name: A string representing the name of the simulation model.
    """
    f_cmake = open('CMakeLists.txt', 'w')
    f_cmake.write(
        'cmake_minimum_required(VERSION 3.1)\n'
        'project(cgmres_simulator CXX)\n'
        '\n'
        'set(CMAKE_CXX_STANDARD 11)\n'
        'set(CMAKE_CXX_FLAGS "-O3")\n'
        '\n'
        '\n'
        'set(MODEL_DIR ${CMAKE_SOURCE_DIR}/models/'+str(model_name)+')\n'
        'set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/cgmres)\n'
        'set(SIMULATOR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/cgmres/simulator)\n'
        'set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)\n'
        'set(SIMULATOR_SRC_DIR ${CMAKE_SOURCE_DIR}/src/simulator)\n'
        '\n'
        'include_directories(${MODEL_DIR})\n'
        'include_directories(${INCLUDE_DIR})\n'
        'include_directories(${SIMULATOR_INCLUDE_DIR})\n'
        '\n'
        '\n'
        'add_subdirectory(${MODEL_DIR})\n'
        '\n'
    )
    if solver_index == 1:
        f_cmake.write(
            'add_library(\n'
            '    cgmres\n'
            '    STATIC\n'
            '    ${SRC_DIR}/continuation_gmres.cpp\n'
            '    ${SRC_DIR}/single_shooting_continuation.cpp\n'
            '    ${SRC_DIR}/single_shooting_ocp.cpp\n'
            '    ${SRC_DIR}/time_varying_smooth_horizon.cpp\n'
            '    ${SRC_DIR}/cgmres_initializer.cpp\n'
            '    ${SRC_DIR}/zero_horizon_ocp.cpp\n'
            '    ${SRC_DIR}/optimal_control_problem.cpp\n'
            '    ${SRC_DIR}/linear_algebra.cpp\n'
            ')\n'
            'target_include_directories(\n'
            '    cgmres\n'
            '    PRIVATE\n'
            '    ${MODEL_DIR}\n'
            '    ${INCLUDE_DIR}\n'
            ')\n'
            '\n'
        )
    elif solver_index == 2:
        f_cmake.write(
            'add_library(\n'
            '    multiple_shooting_cgmres\n'
            '    STATIC\n'
            '    ${SRC_DIR}/multiple_shooting_cgmres.cpp\n'
            '    ${SRC_DIR}/multiple_shooting_continuation.cpp\n'
            '    ${SRC_DIR}/multiple_shooting_ocp.cpp\n'
            '    ${SRC_DIR}/time_varying_smooth_horizon.cpp\n'
            '    ${SRC_DIR}/cgmres_initializer.cpp\n'
            '    ${SRC_DIR}/zero_horizon_ocp.cpp\n'
            '    ${SRC_DIR}/optimal_control_problem.cpp\n'
            '    ${SRC_DIR}/linear_algebra.cpp\n'
            ')\n'
            'target_include_directories(\n'
            '    multiple_shooting_cgmres\n'
            '    PRIVATE\n'
            '    ${MODEL_DIR}\n'
            '    ${INCLUDE_DIR}\n'
            ')\n'
            '\n'
        )
    elif solver_index == 3:
        f_cmake.write(
                'add_library(\n'
            '    ms_cgmres_with_input_saturation\n'
            '    STATIC\n'
            '    ${SRC_DIR}/ms_cgmres_with_input_saturation.cpp\n'
            '    ${SRC_DIR}/ms_continuation_with_input_saturation.cpp\n'
            '    ${SRC_DIR}/ms_ocp_with_input_saturation.cpp\n'
            '    ${SRC_DIR}/time_varying_smooth_horizon.cpp\n'
            '    ${SRC_DIR}/ms_cgmres_with_input_saturation_initializer.cpp\n'
            '    ${SRC_DIR}/zero_horizon_ocp_with_input_saturation.cpp\n'
            '    ${SRC_DIR}/input_saturation_functions.cpp\n'
            '    ${SRC_DIR}/input_saturation_set.cpp\n'
            '    ${SRC_DIR}/input_saturation.cpp\n'
            '    ${SRC_DIR}/optimal_control_problem.cpp\n'
            '    ${SRC_DIR}/linear_algebra.cpp\n'
            ')\n'
            'target_include_directories(\n'
            '    ms_cgmres_with_input_saturation\n'
            '    PRIVATE\n'
            '    ${MODEL_DIR}\n'
            '    ${INCLUDE_DIR}\n'
            ')\n'
            '\n'
        )
    f_cmake.write(
        'add_library(\n'
        '    cgmres_simulator\n'
        '    STATIC\n'
        '    ${SIMULATOR_SRC_DIR}/save_simulation_data.cpp\n'
        '    ${SIMULATOR_SRC_DIR}/numerical_integrator.cpp\n'
        ')\n'
        'target_include_directories(\n'
        '    cgmres_simulator\n'
        '    PRIVATE\n'
        '    ${MODEL_DIR}\n'
        '    ${INCLUDE_DIR}\n'
        '    ${SIMULATOR_INCLUDE_DIR}\n'
        ')\n'
    )
    if platform.system() == 'Windows':
        f_cmake.write(
            'add_executable(main ${MODEL_DIR}/main.cpp)\n'
            'target_link_libraries(main\n'
            '    PRIVATE\n'
            '    cgmres_simulator\n'
        )
        if solver_index == 1:
            f_cmake.write(
                '    cgmres\n'
            )
        elif solver_index == 2:
            f_cmake.write(
                '    multiple_shooting_cgmres\n'
            )
        elif solver_index == 3:
            f_cmake.write(
                '    ms_cgmres_with_input_saturation\n'
            )
        f_cmake.write(
            '    nmpcmodel\n'
            ')\n'
            'target_compile_options(\n'
            '    main\n'
            '    PRIVATE\n'
            '    -O3\n'
            ')\n'
            '\n'
        )
    else:
        f_cmake.write(
            'add_executable(a.out ${MODEL_DIR}/main.cpp)\n'
            'target_link_libraries(a.out\n'
            '    PRIVATE\n'
            '    cgmres_simulator\n'
        )
        if solver_index == 1:
            f_cmake.write(
                '    cgmres\n'
            )
        elif solver_index == 2:
            f_cmake.write(
                '    multiple_shooting_cgmres\n'
            )
        elif solver_index == 3:
            f_cmake.write(
                '    ms_cgmres_with_input_saturation\n'
            )
        f_cmake.write(
            '    nmpcmodel\n'
            ')\n'
            'target_compile_options(\n'
            '    a.out\n'
            '    PRIVATE\n'
            '    -O3\n'
            ')\n'
        )
    f_cmake.close()


def generate_cmake_for_model(model_name):
    """ Generates CMakeLists.txt in a directory where C++ source and header 
        files representing NMPC model are generated.

        Args: 
            model_name: A string representing the name of the simulation model.
    """
    f_cmake = open('models/'+str(model_name)+'/CMakeLists.txt', 'w')
    f_cmake.write(
        'add_library(\n'
        '    nmpcmodel \n'
        '    STATIC \n'
        '    nmpc_model.cpp \n'
        ')'
    )