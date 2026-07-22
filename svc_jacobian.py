#%%
import sympy as sym
from src.sh import *
from src.io import *

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

dI_re_dVk=sym.simplify(sym.diff(I_re,Vk))

dI_re_dthk=sym.simplify(sym.diff(I_re,thk))

dI_re_dbsvc=sym.simplify(sym.diff(I_re,b_svc))

dI_im_dVk=sym.simplify(sym.diff(I_im,Vk))

dI_im_dthk=sym.simplify(sym.diff(I_im,thk))

dI_im_dbsvc=sym.simplify(sym.diff(I_im,b_svc))

dP_dVk=sym.simplify(sym.diff(P,Vk))

dP_dthk=sym.simplify(sym.diff(P,thk))

dP_dbsvc=sym.simplify(sym.diff(P,b_svc))

dQ_dVk=sym.simplify(sym.diff(Q,Vk))

dQ_dthk=sym.simplify(sym.diff(Q,thk))

dQ_dbsvc=sym.simplify(sym.diff(Q,b_svc))
# %%
#LATEX equations
print("Latex equations")
print("Active power")
print("Vk")
print(dP_dVk)
print(r"tk")
print(dP_dthk)
print("bsvc")
print(dP_dbsvc)
print("Reactive power")
print(dQ_dVk)
print(r"tk$")
print(dQ_dthk)
print("bsvc")
print(dQ_dbsvc)
print("Current re")
print("Vk")
print(dI_re_dVk)
print(r"tk")
print(dI_re_dthk)
print("bsvc")
print(dI_re_dbsvc)
print("Current im")
print(dI_im_dVk)
print(r"tk")
print(dI_im_dthk)
print("bsvc")
print(dI_im_dbsvc)
#%% Julia code 
mapping = {
    Vk: sym.Symbol("svc.bus.vm"),

    thk: sym.Symbol("svc.bus.va"),
    b_svc: sym.Symbol("svc.bsvc"),

    rt: sym.Symbol("svc.bsvc*svc.rtf"),
    xt: sym.Symbol("svc.bsvc*svc.xtf"),

}

# %%
print("Julia code")
print("P")
print("Vk")
print_julia(dP_dVk,mapping)
print(r"tk")
print_julia(dP_dthk,mapping)  
print("bsvc")
print_julia(dP_dbsvc,mapping)   
print("Reactive power")
print_julia(dQ_dVk,mapping)
print(r"tk")
print_julia(dQ_dthk,mapping)
print("bsvc")
print_julia(dQ_dbsvc,mapping)
print("Current re")
print("Vk")
print_julia(dI_re_dVk,mapping)
print(r"tk")
print_julia(dI_re_dthk,mapping)
print("bsvc")
print_julia(dI_re_dbsvc,mapping)
print("Current im")
print("Vk")
print_julia(dI_im_dVk,mapping)
print(r"tk")
print_julia(dI_im_dthk,mapping)
print("bsvc")
print_julia(dI_im_dbsvc,mapping)

# %%
