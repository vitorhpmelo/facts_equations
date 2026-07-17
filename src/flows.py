import sympy as sym
def calc_currents_flows(cIkm, cImk):
    """
    Compute the real and imaginary components of the branch currents.

    Parameters
    ----------
    cIkm : sympy.Expr
        Complex current flowing from bus k to bus m.
    cImk : sympy.Expr
        Complex current flowing from bus m to bus k.

    Returns
    -------
    tuple
        (Ikm_re, Ikm_im, Imk_re, Imk_im)
    """
    Ikm_re = sym.simplify(sym.re(cIkm))
    Ikm_im = sym.simplify(sym.im(cIkm))

    Imk_re = sym.simplify(sym.re(cImk))
    Imk_im = sym.simplify(sym.im(cImk))

    return Ikm_re, Ikm_im, Imk_re, Imk_im


def calc_power_flows(cVk, cVm, cIkm, cImk):
    """
    Compute the active and reactive power flows at both terminals of a branch.

    Parameters
    ----------
    cVk : sympy.Expr
        Complex voltage at bus k.
    cVm : sympy.Expr
        Complex voltage at bus m.
    cIkm : sympy.Expr
        Complex current flowing from bus k to bus m.
    cImk : sympy.Expr
        Complex current flowing from bus m to bus k.

    Returns
    -------
    tuple
        (Pkm, Qkm, Pmk, Qmk)
    """
    Skm = sym.simplify(cVk * sym.conjugate(cIkm))
    Smk = sym.simplify(cVm * sym.conjugate(cImk))

    Pkm = sym.simplify(sym.re(Skm))
    Qkm = sym.simplify(sym.im(Skm))

    Pmk = sym.simplify(sym.re(Smk))
    Qmk = sym.simplify(sym.im(Smk))

    return Pkm, Qkm, Pmk, Qmk