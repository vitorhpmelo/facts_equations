#%%
import sympy as sym
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

cIse = -cImk
cIsh = cIkm - cIse

Sse = cVse * sym.conjugate(cImk)
Ssh = cVsh * sym.conjugate(cIsh)

Pse= sym.simplify(sym.re(Sse))
Qse= sym.simplify(sym.im(Sse))
Psh= sym.simplify(sym.re(Ssh))
Qsh= sym.simplify(sym.im(Ssh))
Pb= sym.simplify(Pse+Psh)


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
