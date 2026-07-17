#%%
import sympy as sym
from src.flows import *

# Define the symbolic state variable representing the TCSC series reactance.
x_tcsc = sym.Symbol("x_{tcsc}", real=True)

# Define the voltage magnitudes at buses k and m.
Vk, Vm = sym.symbols(r"V_k V_m", positive=True)

# Define the voltage phase angles at buses k and m.
thk, thm = sym.symbols(r"\theta_k \theta_m", real=True)

#%%
# Construct the complex bus voltages in polar form:
#   V = |V|e^{jθ}
cVk = Vk * sym.exp(sym.I * thk)
cVm = Vm * sym.exp(sym.I * thm)

#%%
# Compute the branch currents using Ohm's law:
#   I = (V_k - V_m)/(jx_tcsc)
cIkm = (cVk - cVm) / (sym.I * x_tcsc)
cImk = (cVm - cVk) / (sym.I * x_tcsc)


# Compute the current components.
Ikm_re, Ikm_im, Imk_re, Imk_im = calc_currents_flows(cIkm, cImk)

# Compute the power flows.
Pkm, Qkm, Pmk, Qmk = calc_power_flows(cVk, cVm, cIkm, cImk)
# %%
