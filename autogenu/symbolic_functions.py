import sympy


def dot_product(vec1, vec2):
    """ Calculate dot product of two symbolic vectors or vector valued 
        functions. 

        Args:
            vec1, vec2: Symbolic vectors or vector valued functions.

        Returns: 
            Dot product of two symbolic vectors or vector valued functions.
    """
    return sum(vec1[i]*vec2[i] for i in range(len(vec1)))


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


def simplify_vector_func(vector_func):
    """ Simplifies a vector-valued function.

        Args:
            vector_func: A list of SymPy scalar functions.
    """
    for i in range(len(vector_func)):
        vector_func[i] = sympy.simplify(sympy.nsimplify(vector_func[i]))