import sympy


def diff_scalar_func(scalar_func, var):
    """ Calculate partial derivative of a function with respect to a scalar or
        a vector. 

        Args:
            scalar_func: A symbolic scalar function.
            var: A symbolic scalar or a symbolic vector.

        Returns: 
            Partial derivative of scalar_func with respect to var. If var is a 
            vector, Returns Jacobian.
    """
    return [sympy.diff(scalar_func, var[i]) for i in range(len(var))]


def simplify(func):
    """ Simplifies a scalar-valued or vector-valued function.

        Args:
            func: A symbolic functions.
    """
    if type(func) == list:
        for i in range(len(func)):
            func[i] = sympy.simplify(sympy.nsimplify(func[i]))
    else:
        func = sympy.simplify(sympy.nsimplify(func))

def write_symfunc(writable_file, function, output_value_name: str, common_subexpression_elimination: bool):
    """ Write input symbolic function onto writable_file. The function's 
        return value name must be set. common_subexpression_elimination is optional.

        Args: 
            writable_file: A writable file, i.e., a file streaming that is 
                already opened as writing mode.
            function: A symbolic function wrote onto the writable_file.
            output_value_name: The name of the output value.
            common_subexpression_elimination: If true, common subexpression elimination is used. If 
                False, it is not used.
    """
    if common_subexpression_elimination:
        func_cse = sympy.cse(function)
        for i in range(len(func_cse[0])):
            cse_exp, cse_rhs = func_cse[0][i]
            writable_file.write(
                '    const double '+sympy.ccode(cse_exp)
                +' = '+sympy.ccode(cse_rhs)+';\n'
            )
        for i in range(len(func_cse[1])):
            writable_file.write(
                '    '+output_value_name+'[%d] = '%i
                +sympy.ccode(func_cse[1][i])+';\n'
            )
    else:
        writable_file.writelines(
            ['    '+output_value_name+'[%d] = '%i
            +sympy.ccode(function[i])+';\n' for i in range(len(function))]
        )