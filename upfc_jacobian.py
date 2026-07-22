#%%
import sympy as sym
from sympy import Symbol
from src.flows import *
from src.io import *




gt_se,bt_se,gt_sh,bt_sh= sym.symbols("g_{tf}^{se} b_{tf}^{se} g_{tf}^{sh} b_{tf}^{sh}", real=True)
# %%
# Define the voltage magnitudes at buses k and m.
Vk, Vm, Vsh, Vse = sym.symbols(r"V_k V_m V_{sh} V_{se}", positive=True)

# Define the voltage phase angles at buses k and m.
thk, thm, thsh, thse = sym.symbols(r"\theta_k \theta_m \theta_{sh} \theta_{se}", real=True)




cVk = Vk * sym.exp(sym.I * thk)
cVm = Vm * sym.exp(sym.I * thm)

cVse = Vse * sym.exp(sym.I * thse)
cVsh = Vsh * sym.exp(sym.I * thsh)

yse = gt_se + sym.I*bt_se
ysh = gt_sh + sym.I*bt_sh

V=[cVk,cVm,cVse,cVsh]
Y = [[yse+ysh,-yse,-yse,-ysh],[-yse,yse,yse,0]]

cI = sym.Matrix(Y) * sym.Matrix(V)
cIkm = sym.sympify(cI[0])
cImk = sym.sympify(cI[1])
# %%
Ikm_re, Ikm_im, Imk_re, Imk_im = calc_currents_flows(cIkm, cImk)

# Compute the power flows.
Pkm, Qkm, Pmk, Qmk = calc_power_flows(cVk, cVm, cIkm, cImk)


cIse = -cImk
cIsh = sym.simplify((cIkm+cImk))
#%%
Sse = cVse * sym.conjugate(cIse)
Ssh = cVsh * sym.conjugate(cIsh)

Pse= sym.simplify(sym.re(Sse))
Qse= sym.simplify(sym.im(Sse))
Psh= sym.simplify(sym.re(Ssh))
Qsh= sym.simplify(sym.im(Ssh))
Pb= sym.simplify(Pse+Psh)


#%%

dPkm_dthk=sym.simplify(sym.diff(Pkm,thk))
dPkm_dthm=sym.simplify(sym.diff(Pkm,thm))
dPkm_dtse=sym.simplify(sym.diff(Pkm,thse))
dPkm_dtsh=sym.simplify(sym.diff(Pkm,thsh))
dPmk_dthk=sym.simplify(sym.diff(Pmk,thk))
dPmk_dthm=sym.simplify(sym.diff(Pmk,thm))
dPmk_dtse=sym.simplify(sym.diff(Pmk,thse))
dPmk_dtsh=sym.simplify(sym.diff(Pmk,thsh))


dQkm_dthk=sym.simplify(sym.diff(Qkm,thk))
dQkm_dthm=sym.simplify(sym.diff(Qkm,thm))
dQkm_dtse=sym.simplify(sym.diff(Qkm,thse))
dQkm_dtsh=sym.simplify(sym.diff(Qkm,thsh))
dQmk_dthk=sym.simplify(sym.diff(Qmk,thk))
dQmk_dthm=sym.simplify(sym.diff(Qmk,thm))
dQmk_dtse=sym.simplify(sym.diff(Qmk,thse))
dQmk_dtsh=sym.simplify(sym.diff(Qmk,thsh))



dIkm_re_dthk=sym.simplify(sym.diff(Ikm_re,thk))
dIkm_re_dthm=sym.simplify(sym.diff(Ikm_re,thm))
dIkm_re_dtse=sym.simplify(sym.diff(Ikm_re,thse))
dIkm_re_dtsh=sym.simplify(sym.diff(Ikm_re,thsh))
dImk_re_dthk=sym.simplify(sym.diff(Imk_re,thk))
dImk_re_dthm=sym.simplify(sym.diff(Imk_re,thm))
dImk_re_dtse=sym.simplify(sym.diff(Imk_re,thse))
dImk_re_dtsh=sym.simplify(sym.diff(Imk_re,thsh))


dIkm_im_dthk=sym.simplify(sym.diff(Ikm_im,thk))
dIkm_im_dthm=sym.simplify(sym.diff(Ikm_im,thm))
dIkm_im_dtse=sym.simplify(sym.diff(Ikm_im,thse))
dIkm_im_dtsh=sym.simplify(sym.diff(Ikm_im,thsh))
dImk_im_dthk=sym.simplify(sym.diff(Imk_im,thk))
dImk_im_dthm=sym.simplify(sym.diff(Imk_im,thm))
dImk_im_dtse=sym.simplify(sym.diff(Imk_im,thse))
dImk_im_dtsh=sym.simplify(sym.diff(Imk_im,thsh))


dPb_dthk=sym.simplify(sym.diff(Pb,thk))
dPb_dthm=sym.simplify(sym.diff(Pb,thm))
dPb_dtse=sym.simplify(sym.diff(Pb,thse))
dPb_dtsh=sym.simplify(sym.diff(Pb,thsh))




mapping = {
    Vk:   sym.Symbol("upfc.from_bus.vm"),
    Vm:   sym.Symbol("upfc.virtual_bus.vm"),
    Vse:  sym.Symbol("upfc.v_se.vm"),
    Vsh:  sym.Symbol("upfc.v_sh.vm"),

    thk:  sym.Symbol("upfc.from_bus.va"),
    thm:  sym.Symbol("upfc.virtual_bus.va"),
    thse: sym.Symbol("upfc.v_se.va"),
    thsh: sym.Symbol("upfc.v_sh.va"),

    gt_se: sym.Symbol("upfc.tf_se.G[1,1]"),
    bt_se: sym.Symbol("upfc.tf_se.B[1,1]"),
    gt_sh: sym.Symbol("upfc.tf_sh.G[1,1]"),
    bt_sh: sym.Symbol("upfc.tf_sh.B[1,1]"),
}



#%%

print("Julia code")
print("=" * 80)

# ============================================================================
# Active power flow Pkm
# ============================================================================
print("\nActive power Pkm")
print("θk")
print(print_julia(dPkm_dthk, mapping))

print("θm")
print(print_julia(dPkm_dthm, mapping))

print("θse")
print(print_julia(dPkm_dtse, mapping))

print("θsh")
print(print_julia(dPkm_dtsh, mapping))


# ============================================================================
# Active power flow Pmk
# ============================================================================
print("\nActive power Pmk")
print("θk")
print(print_julia(dPmk_dthk, mapping))

print("θm")
print(print_julia(dPmk_dthm, mapping))

print("θse")
print(print_julia(dPmk_dtse, mapping))

print("θsh")
print(print_julia(dPmk_dtsh, mapping))


# ============================================================================
# Reactive power flow Qkm
# ============================================================================
print("\nReactive power Qkm")
print("θk")
print(print_julia(dQkm_dthk, mapping))

print("θm")
print(print_julia(dQkm_dthm, mapping))

print("θse")
print(print_julia(dQkm_dtse, mapping))

print("θsh")
print(print_julia(dQkm_dtsh, mapping))


# ============================================================================
# Reactive power flow Qmk
# ============================================================================
print("\nReactive power Qmk")
print("θk")
print(print_julia(dQmk_dthk, mapping))

print("θm")
print(print_julia(dQmk_dthm, mapping))

print("θse")
print(print_julia(dQmk_dtse, mapping))

print("θsh")
print(print_julia(dQmk_dtsh, mapping))


# ============================================================================
# Current Ikm (real)
# ============================================================================
print("\nCurrent Ikm (real)")
print("θk")
print(print_julia(dIkm_re_dthk, mapping))

print("θm")
print(print_julia(dIkm_re_dthm, mapping))

print("θse")
print(print_julia(dIkm_re_dtse, mapping))

print("θsh")
print(print_julia(dIkm_re_dtsh, mapping))


# ============================================================================
# Current Imk (real)
# ============================================================================
print("\nCurrent Imk (real)")
print("θk")
print(print_julia(dImk_re_dthk, mapping))

print("θm")
print(print_julia(dImk_re_dthm, mapping))

print("θse")
print(print_julia(dImk_re_dtse, mapping))

print("θsh")
print(print_julia(dImk_re_dtsh, mapping))


# ============================================================================
# Current Ikm (imaginary)
# ============================================================================
print("\nCurrent Ikm (imaginary)")
print("θk")
print(print_julia(dIkm_im_dthk, mapping))

print("θm")
print(print_julia(dIkm_im_dthm, mapping))

print("θse")
print(print_julia(dIkm_im_dtse, mapping))

print("θsh")
print(print_julia(dIkm_im_dtsh, mapping))


# ============================================================================
# Current Imk (imaginary)
# ============================================================================
print("\nCurrent Imk (imaginary)")
print("θk")
print(print_julia(dImk_im_dthk, mapping))

print("θm")
print(print_julia(dImk_im_dthm, mapping))

print("θse")
print(print_julia(dImk_im_dtse, mapping))

print("θsh")
print(print_julia(dImk_im_dtsh, mapping))


# ============================================================================
# Power balance Pb
# ============================================================================
print("\nPower balance Pb")
print("θk")
print(print_julia(dPb_dthk, mapping))

print("θm")
print(print_julia(dPb_dthm, mapping))

print("θse")
print(print_julia(dPb_dtse, mapping))

print("θsh")
print(print_julia(dPb_dtsh, mapping))