import sympy


def dotProduct(vec1, vec2):
    return sum(vec1[i]*vec2[i] for i in range(len(vec1)))

def diffScalarFunc(scalar_func, var):
    return [sympy.diff(scalar_func, var[i]) for i in range(len(var))]