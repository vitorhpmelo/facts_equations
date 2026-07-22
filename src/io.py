import sympy as sym 
from sympy.printing.julia import julia_code

def print_julia(expr, mapping=None):
    """Return the Julia code corresponding to a SymPy expression."""
    if mapping is not None:
        expr = expr.xreplace(mapping)
    return julia_code(expr)


def print_latex(expr, mapping=None):
    """Return the LaTeX code corresponding to a SymPy expression."""
    if mapping is not None:
        expr = expr.xreplace(mapping)
    return sym.latex(expr)