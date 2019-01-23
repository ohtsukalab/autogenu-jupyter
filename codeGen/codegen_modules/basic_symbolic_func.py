import sympy


def dotProduct(vec1, vec2):
    return sum(vec1[i]*vec2[i] for i in range(len(vec1)))

def diffScalarFunc(scalar_func, var):
    return [sympy.diff(scalar_func, var[i]) for i in range(len(var))]

def generateBasicVars(dimx, dimu, dimc):
    t = sympy.Symbol('t')
    x = sympy.symbols(f'x[0:{dimx}]')
    u = sympy.symbols(f'u[0:{dimu+dimc}]')
    lmd = sympy.symbols(f'lmd[0:{dimx}]')

    return t, x, u, lmd