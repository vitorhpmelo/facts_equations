#%%
import sympy as sym
from sympy import Symbol
from sympy.printing.julia import julia_code

from src.flows import *




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
# %%
#cIse= -cImk

#cIkm = cIsh + (-cImk)
#cIsh = cIkm + cImk

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
print("Currents re:")
print("direction from")
print(julia_code(Ikm_re))
print("direction to")
print(julia_code(Imk_re))


print("Currents im:")
print("direction from")
print(julia_code(Ikm_im))
print("direction to")
print(julia_code(Imk_im))

print("Active Powers:")
print("direction from")
print(julia_code(Pkm))
print("direction to")
print(julia_code(Pmk))
print("Reactive Powers:")
print("direction from")
print(julia_code(Qkm))
print("direction to")
print(julia_code(Qmk))
# %%

print("Power Balance:")
print(julia_code(Pb))

# %%


Pkm.subs({Vk: 1.022609222019097, Vm: 1.0275274934266296, Vse: 0.010403779466191827, Vsh: 1, thk: -0.23907669146542013, thm: -0.23987960384008997, thse: -2.7045447412782946, thsh: -0.23691105365751117, gt_se: 0.99009900990099, bt_se: -9.900990099009901, gt_sh: 0.99009900990099, bt_sh: -9.900990099009901})

print("Power Balance:")
Pb.subs({Vk: 1.022609222019097, Vm: 1.0275274934266296, Vse: 0.010403779466191827, Vsh: 1, thk: -0.23907669146542013, thm: -0.23987960384008997, thse: -2.7045447412782946, thsh: -0.23691105365751117, gt_se: 0.99009900990099, bt_se: -9.900990099009901, gt_sh: 0.99009900990099, bt_sh: -9.900990099009901})    
#%%
mapping = {
    Vk: Symbol("upfc.from_bus.vm"),
    Vm: Symbol("upfc.virtual_bus.vm"),
    Vse: Symbol("upfc.v_se.vm"),
    Vsh: Symbol("upfc.v_sh.vm"),

    thk: Symbol("upfc.from_bus.va"),
    thm: Symbol("upfc.virtual_bus.va"),
    thse: Symbol("upfc.v_se.va"),
    thsh: Symbol("upfc.v_sh.va"),

    gt_se: Symbol("upfc.tf_se.G[1,1]"),
    bt_se: Symbol("upfc.tf_se.B[1,1]"),
    gt_sh: Symbol("upfc.tf_sh.G[1,1]"),
    bt_sh: Symbol("upfc.tf_sh.B[1,1]"),
}
#%%


expr = Pkm.xreplace(mapping)

exprPb = Pb.xreplace(mapping)

print("Pkm in julia code:")
print(julia_code(expr))

print("Pb in julia code:")
print(julia_code(exprPb))
# %%
expr.subs({Vk: 1.022609222019097,
            Vm: 1.022609222019097, 
            Vse: 0.010403779466191827, 
            Vsh: 1,
            thk: -0.23907669146542013,
            thm: -0.23987960384008997, 
            thse: -2.7045447412782946, 
            thsh: -0.23691105365751117, 
            gt_se: 0.99009900990099,
            bt_se: -9.900990099009901, 
            gt_sh: 0.99009900990099,
            bt_sh: -9.900990099009901})
# %%
