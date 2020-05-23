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