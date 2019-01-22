import sympy
import sys
import linecache

def generateCpp(dimx, dimu, dimc, fxu, Cxu, phix, hx, hu):
    f_model_cpp = open('nmpc_model.cpp', 'w')
    f_model_cpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.cpp', i) for i in range(0,11)])
    f_model_cpp.writelines(['    f[%d] = '%i + sympy.ccode(fxu[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.cpp', i) for i in range(12,22)])
    f_model_cpp.writelines(['    phix[%d] = '%i + sympy.ccode(phix[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.cpp', i) for i in range(23,36)])
    f_model_cpp.writelines(['    hx[%d] = '%i + sympy.ccode(hx[i]) + ';\n' for i in range(dimx)])
    f_model_cpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.cpp', i) for i in range(37,50)])
    f_model_cpp.writelines(['    hu[%d] = '%i + sympy.ccode(hu[i]) + ';\n' for i in range(dimu)])
    f_model_cpp.write('}')
    f_model_cpp.close()


def generateHpp(dimx, dimu, dimc, scalar_params, array_params):
    f_model_hpp = open('nmpc_model.hpp', 'w')
    f_model_hpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.hpp', i) for i in range(0,18)])

    f_model_hpp.write('    static constexpr int dim_state_ = %d;\n' %dimx)
    f_model_hpp.write('    static constexpr int dim_control_input_ = %d;\n' %dimu)
    f_model_hpp.write('    static constexpr int dim_constraints_ = %d;\n' %dimc)
    f_model_hpp.write('\n')

    f_model_hpp.writelines(['    static constexpr double ' + str(scalar_params[i][0]) + ' = ' + str(scalar_params[i][1]) + ';\n' for i in range(len(scalar_params))])
    f_model_hpp.write('\n\n')

    f_model_hpp.writelines(['    double ' + str(array_params[i][0]) + '[' + str(array_params[i][1]) + ']' + ' = ' + array_params[i][2] + ';\n' for i in range(len(array_params))])

    f_model_hpp.writelines([linecache.getline('templates_for_codeGen/nmpc_model.hpp', i) for i in range(38,71)])

    f_model_hpp.close()
