#%%
import sympy as sym
from src.sh import *

# Define the symbolic state variable representing the SVC series reactance.
rt,xt,b_svc = sym.symbols("r_t x_t b_{svc}", real=True)

xeq = xt - 1/b_svc

gk = rt/((rt**2 )+ (xeq**2))
bk = -xeq/((rt**2 )+ (xeq**2))
yk = gk+sym.I*bk

# Define the voltage magnitudes at bus k.
Vk = sym.symbols(r"V_k", positive=True)

# Define the voltage phase angles at bus k.
thk = sym.symbols(r"\theta_k", real=True)

# Construct the complex bus voltage in polar form:
#   V = |V|e^{jθ}
cVk = Vk * sym.exp(sym.I * thk)

cIsh = yk*cVk

I_re, I_im = calc_currents_sh(cIsh)

# Power measurements
P, Q = calc_power_sh(cVk, cIsh)
#%%