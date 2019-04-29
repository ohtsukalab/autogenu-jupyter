from AutoGenU_modules import solver_params as solverparam
from AutoGenU_modules import simulation_params as simparam
from AutoGenU_modules import cpp_executor
import sympy
import linecache
import subprocess
import platform


def makeModelDir(model_name):
    if(platform.system() == 'Windows'):
        subprocess.run(['mkdir', 'models'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
        subprocess.run(['mkdir', model_name], cwd='models', stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    else:
        subprocess.run(['mkdir', 'models'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        subprocess.run(['mkdir', model_name], cwd='models', stdout = subprocess.PIPE, stderr = subprocess.PIPE)



def generateCpp(dimx, dimu, dimc, fxu, Cxu, phix, hx, hu, model_name):
    f_model_cpp = open('models/'+str(model_name)+'/nmpc_model.cpp', 'w')
    f_model_cpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.cpp', i) for i in range(0,11)])
    f_model_cpp.writelines(['    f[%d] = '%i + sympy.ccode(fxu[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.cpp', i) for i in range(12,22)])
    f_model_cpp.writelines(['    phix[%d] = '%i + sympy.ccode(phix[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.cpp', i) for i in range(23,36)])
    f_model_cpp.writelines(['    hx[%d] = '%i + sympy.ccode(hx[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.cpp', i) for i in range(37,50)])
    f_model_cpp.writelines(['    hu[%d] = '%i + sympy.ccode(hu[i]) + ';\n' for i in range(dimu)])
    f_model_cpp.write('}')
    f_model_cpp.close()


def generateHpp(dimx, dimu, dimc, scalar_params, array_params, model_name):
    f_model_hpp = open('models/'+str(model_name)+'/nmpc_model.hpp', 'w')
    f_model_hpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.hpp', i) for i in range(0,18)])
    f_model_hpp.write('    static constexpr int dim_state_ = %d;\n' %dimx)
    f_model_hpp.write('    static constexpr int dim_control_input_ = %d;\n' %dimu)
    f_model_hpp.write('    static constexpr int dim_constraints_ = %d;\n' %dimc)
    f_model_hpp.write('\n')
    f_model_hpp.writelines(['    static constexpr double ' + str(scalar_params[i][0]) + ' = ' + str(scalar_params[i][1]) + ';\n' for i in range(len(scalar_params))])
    f_model_hpp.write('\n\n')
    f_model_hpp.writelines(['    double ' + str(array_params[i][0]) + '[' + str(array_params[i][1]) + ']' + ' = ' + str(array_params[i][2]) + ';\n' for i in range(len(array_params))])
    f_model_hpp.writelines([linecache.getline('AutoGenU_modules/.cpp_templates/nmpc_model.hpp', i) for i in range(18,51)])
    f_model_hpp.close()


def generateMain(solver_index, solver_params, initialization_params, simulation_params, model_name, saturation_list=None):
    # include header
    f_main = open('models/'+str(model_name)+'/main.cpp', 'w')
    f_main.write('#include "nmpc_model.hpp"\n')
    if(solver_index == 1):
        f_main.write('#include "continuation_gmres.hpp"\n')
        f_main.write('#include "cgmres_simulator.hpp"\n')
    elif(solver_index == 2):
        f_main.write('#include "multiple_shooting_cgmres.hpp"\n')
        f_main.write('#include "multiple_shooting_cgmres_simulator.hpp"\n')
    else:
        f_main.write('#include "control_input_saturation_sequence.hpp"\n')
        f_main.write('#include "multiple_shooting_cgmres_with_saturation.hpp"\n')
        f_main.write('#include "multiple_shooting_cgmres_with_saturation_simulator.hpp"\n')
    f_main.write('\n')

    f_main.write('int main()\n')
    f_main.write('{\n')

    # define NMPCModel 
    f_main.write('    // Define the model in NMPC.\n')
    f_main.write('    NMPCModel nmpc_model;\n')
    f_main.write('\n')

    # define solver
    f_main.write('    // Define the solver.\n')
    if(solver_index == 1):
        f_main.write('    ContinuationGMRES nmpc_solver('+str(solver_params.T_f)+', '+str(solver_params.alpha)+', '+str(solver_params.horizon_division_num)+', '+str(solver_params.difference_increment)+', '+str(solver_params.zeta)+', '+str(solver_params.max_dim_krylov)+');\n')
    elif(solver_index == 2):
        f_main.write('    MultipleShootingCGMRES nmpc_solver('+str(solver_params.T_f)+', '+str(solver_params.alpha)+', '+str(solver_params.horizon_division_num)+', '+str(solver_params.difference_increment)+', '+str(solver_params.zeta)+', '+str(solver_params.max_dim_krylov)+');\n')
    else:
        f_main.write('    ControlInputSaturationSequence control_input_saturation_seq;\n')
        for i in range(len(saturation_list)):
            f_main.write('    control_input_saturation_seq.appendControlInputSaturation('+str(saturation_list[i][0])+', '+str(saturation_list[i][1])+', '+str(saturation_list[i][2])+', '+str(saturation_list[i][3])+', '+str(saturation_list[i][4])+');\n')
        f_main.write('    MultipleShootingCGMRESWithSaturation nmpc_solver(control_input_saturation_seq, '+str(solver_params.T_f)+', '+str(solver_params.alpha)+', '+str(solver_params.horizon_division_num)+', '+str(solver_params.difference_increment)+', '+str(solver_params.zeta)+', '+str(solver_params.max_dim_krylov)+');\n')
    f_main.write('\n')
    f_main.write('\n')

    # initial state
    f_main.write('    // Set the initial state.\n')
    f_main.write('    double initial_state['+str(len(simulation_params.initial_state))+'] = {')
    for i in range(len(simulation_params.initial_state)-1):
        f_main.write(str(simulation_params.initial_state[i]) + ', ')
    f_main.write(str(simulation_params.initial_state[-1]) + '};\n')    
    f_main.write('\n')

    # initial guess solution
    f_main.write('    // Set the initial guess of the solution.\n')
    f_main.write('    double initial_guess_solution['+str(len(initialization_params.initial_guess_solution))+'] = {')
    for i in range(len(initialization_params.initial_guess_solution)-1):
        f_main.write(str(initialization_params.initial_guess_solution[i]) + ', ')
    f_main.write(str(initialization_params.initial_guess_solution[-1]) + '};\n')    
    f_main.write('\n')
    if solver_index == 3 and initialization_params.initial_guess_lagrange_multiplier != None:
        f_main.write('    // Set the initial guess of the lagrange multiplier for the condensed constraints with respect to the saturation on the function of the control input .\n')
        f_main.write('    double initial_guess_lagrange_multiplier['+str(len(initialization_params.initial_guess_lagrange_multiplier))+'] = {')
        for i in range(len(initialization_params.initial_guess_lagrange_multiplier)-1):
            f_main.write(str(initialization_params.initial_guess_lagrange_multiplier[i]) + ', ')
        f_main.write(str(initialization_params.initial_guess_lagrange_multiplier[-1]) + '};\n')    
    f_main.write('\n')
    f_main.write('\n')

    # initialization of the solution of the C/GMRES method
    f_main.write('    // Initialize the solution of the C/GMRES method.\n')
    if solver_index == 3 and initialization_params.initial_guess_lagrange_multiplier != None:
        f_main.write('    nmpc_solver.initSolution('+str(simulation_params.initial_time)+', initial_state, initial_guess_solution, initial_guess_lagrange_multiplier, ' +str(initialization_params.convergence_radius)+', '+str(initialization_params.maximum_itr_newton)+');\n')
    else:
        f_main.write('    nmpc_solver.initSolution('+str(simulation_params.initial_time)+', initial_state, initial_guess_solution, '+str(initialization_params.convergence_radius)+', '+str(initialization_params.maximum_itr_newton)+');\n')
    f_main.write('\n')

    f_main.write('    // Perform a numerical simulation.\n')
    f_main.write('    nmpcsim::simulation(nmpc_solver, initial_state, '+str(simulation_params.initial_time)+', '+str(simulation_params.simulation_time)+', '+str(simulation_params.sampling_time)+', "'+model_name+'");\n')
    f_main.write('\n')

    f_main.write('    return 0;\n')
    f_main.write('}\n')
    f_main.close()


def generateCMake(solver_index, model_name):
    f_cmake = open('CMakeLists.txt', 'w')
    f_cmake.write('cmake_minimum_required(VERSION 3.1)\n')
    f_cmake.write('project(cgmres_simulator CXX)\n')
    f_cmake.write('\n')
    f_cmake.write('set(CMAKE_CXX_STANDARD 11)\n')
    f_cmake.write('set(CMAKE_CXX_FLAGS "-O3")\n')
    f_cmake.write('\n')
    f_cmake.write('set(SOLVER_DIR ${CMAKE_SOURCE_DIR}/src/solver)\n')
    f_cmake.write('set(SIMULATOR_DIR ${CMAKE_SOURCE_DIR}/src/simulator)\n')
    f_cmake.write('set(MODEL_DIR ${CMAKE_SOURCE_DIR}/'+'models/'+str(model_name)+')\n')
    f_cmake.write('\n')
    f_cmake.write('include_directories(${SOLVER_DIR})\n')
    f_cmake.write('include_directories(${SIMULATOR_DIR})\n')
    f_cmake.write('include_directories(${MODEL_DIR})\n')
    f_cmake.write('\n')
    f_cmake.write('\n')
    f_cmake.write('add_subdirectory(${MODEL_DIR})\n')
    f_cmake.write('\n')

    if(solver_index == 1):
        f_cmake.write('add_library(\n')
        f_cmake.write('    cgmres\n')
        f_cmake.write('    STATIC\n')
        f_cmake.write('    ${SOLVER_DIR}/continuation_gmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/init_cgmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/matrixfree_gmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/linear_funcs.cpp\n')
        f_cmake.write(')\n')
        f_cmake.write('\n')
        f_cmake.write('add_library(\n')
        f_cmake.write('    cgmres_simulator\n')
        f_cmake.write('    STATIC\n')
        if(platform.system() == 'Windows'):
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data_for_windows.cpp\n')
        else:
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/numerical_integrator.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/cgmres_simulator.cpp\n')
        f_cmake.write(')\n')
    elif(solver_index == 2):
        f_cmake.write('add_library(\n')
        f_cmake.write('    multiple_shooting_cgmres\n')
        f_cmake.write('    STATIC\n')
        f_cmake.write('    ${SOLVER_DIR}/multiple_shooting_cgmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/init_cgmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/matrixfree_gmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/linear_funcs.cpp\n')
        f_cmake.write(')\n')
        f_cmake.write('\n')
        f_cmake.write('add_library(\n')
        f_cmake.write('    multiple_shooting_cgmres_simulator\n')
        f_cmake.write('    STATIC\n')
        if(platform.system() == 'Windows'):
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data_for_windows.cpp\n')
        else:
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/numerical_integrator.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/multiple_shooting_cgmres_simulator.cpp\n')
        f_cmake.write(')\n')
    elif(solver_index == 3):
        f_cmake.write('add_library(\n')
        f_cmake.write('    multiple_shooting_cgmres_with_saturation\n')
        f_cmake.write('    STATIC\n')
        f_cmake.write('    ${SOLVER_DIR}/multiple_shooting_cgmres_with_saturation.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/init_cgmres_with_saturation.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/control_input_saturation.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/control_input_saturation_sequence.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/matrixfree_gmres.cpp\n')
        f_cmake.write('    ${SOLVER_DIR}/linear_funcs.cpp\n')
        f_cmake.write(')\n')
        f_cmake.write('\n')
        f_cmake.write('add_library(\n')
        f_cmake.write('    multiple_shooting_cgmres_with_saturation_simulator\n')
        f_cmake.write('    STATIC\n')
        if(platform.system() == 'Windows'):
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data_for_windows.cpp\n')
        else:
            f_cmake.write('    ${SIMULATOR_DIR}/save_simulation_data.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/numerical_integrator.cpp\n')
        f_cmake.write('    ${SIMULATOR_DIR}/multiple_shooting_cgmres_with_saturation_simulator.cpp\n')
        f_cmake.write(')\n')

    f_cmake.write('\n')
    f_cmake.write('\n')
    if(platform.system() == 'Windows'):
        f_cmake.write('add_executable(a.exe ${MODEL_DIR}/main.cpp)\n')
        if(solver_index == 1):
            f_cmake.write('target_link_libraries(a.exe cgmres)\n')
            f_cmake.write('target_link_libraries(a.exe cgmres_simulator)\n')
        elif(solver_index == 2):
            f_cmake.write('target_link_libraries(a.exe multiple_shooting_cgmres)\n')
            f_cmake.write('target_link_libraries(a.exe multiple_shooting_cgmres_simulator)\n')
        elif(solver_index == 3):
            f_cmake.write('target_link_libraries(a.exe multiple_shooting_cgmres_with_saturation)\n')
            f_cmake.write('target_link_libraries(a.exe multiple_shooting_cgmres_with_saturation_simulator)\n')
        f_cmake.write('target_link_libraries(a.exe nmpcmodel)\n')
    else:
        f_cmake.write('add_executable(a.out ${MODEL_DIR}/main.cpp)\n')
        if(solver_index == 1):
            f_cmake.write('target_link_libraries(a.out cgmres)\n')
            f_cmake.write('target_link_libraries(a.out cgmres_simulator)\n')
        elif(solver_index == 2):
            f_cmake.write('target_link_libraries(a.out multiple_shooting_cgmres)\n')
            f_cmake.write('target_link_libraries(a.out multiple_shooting_cgmres_simulator)\n')
        elif(solver_index == 3):
            f_cmake.write('target_link_libraries(a.out multiple_shooting_cgmres_with_saturation)\n')
            f_cmake.write('target_link_libraries(a.out multiple_shooting_cgmres_with_saturation_simulator)\n')
        f_cmake.write('target_link_libraries(a.out nmpcmodel)\n')


    f_cmake.close()


def generateCMakeForModel(model_name):
    f_cmake = open('models/'+str(model_name)+'/CMakeLists.txt', 'w')
    f_cmake.write('add_library(\n')
    f_cmake.write('    nmpcmodel \n')
    f_cmake.write('    STATIC \n')
    f_cmake.write('    nmpc_model.cpp \n')
    f_cmake.write(')')
