import sympy as sym

def calc_currents_sh(cIsh):
    """
    Compute the real and imaginary components of a shunt current.

    Parameters
    ----------
    cIsh : sympy.Expr
        Complex current injected by the shunt device.

    Returns
    -------
    tuple
        (I_re, I_im)
    """
    I_re = sym.simplify(sym.re(cIsh))
    I_im = sym.simplify(sym.im(cIsh))

    return I_re, I_im



def calc_power_sh(cV, cIsh):
    """
    Compute the active and reactive power injection of a shunt device.

    Parameters
    ----------
    cV : sympy.Expr
        Complex voltage at the connected bus.
    cIsh : sympy.Expr
        Complex current injected by the shunt device.

    Returns
    -------
    tuple
        (P, Q)
    """
    S = sym.simplify(cV * sym.conjugate(cIsh))

    P = sym.simplify(sym.re(S))
    Q = sym.simplify(sym.im(S))

    return P, Q