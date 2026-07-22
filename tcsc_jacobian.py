#%%
import sympy as sym
from src.flows import *
from src.io import print_julia,print_latex

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
mapping = {
    Vk: sym.Symbol("tcsc.from_bus.vm"),
    Vm: sym.Symbol("tcsc.virtual_bus.vm"),


    thk: sym.Symbol("tcsc.from_bus.va"),
    thm: sym.Symbol("tcsc.virtual_bus.va"),

    x_tcsc: sym.Symbol("tcsc.x_tcsc")

}

#%%
# ============================================================================
# Jacobian entries
# ============================================================================

# ---------------------------------------------------------------------------
# Active power flow Pkm
# ---------------------------------------------------------------------------
dPkm_dVk     = sym.simplify(sym.diff(Pkm, Vk))
dPkm_dVm     = sym.simplify(sym.diff(Pkm, Vm))
dPkm_dthk    = sym.simplify(sym.diff(Pkm, thk))
dPkm_dthm    = sym.simplify(sym.diff(Pkm, thm))
dPkm_dx_tcsc = sym.simplify(sym.diff(Pkm, x_tcsc))

# ---------------------------------------------------------------------------
# Reactive power flow Qkm
# ---------------------------------------------------------------------------
dQkm_dVk     = sym.simplify(sym.diff(Qkm, Vk))
dQkm_dVm     = sym.simplify(sym.diff(Qkm, Vm))
dQkm_dthk    = sym.simplify(sym.diff(Qkm, thk))
dQkm_dthm    = sym.simplify(sym.diff(Qkm, thm))
dQkm_dx_tcsc = sym.simplify(sym.diff(Qkm, x_tcsc))

# ---------------------------------------------------------------------------
# Active power flow Pmk
# ---------------------------------------------------------------------------
dPmk_dVk     = sym.simplify(sym.diff(Pmk, Vk))
dPmk_dVm     = sym.simplify(sym.diff(Pmk, Vm))
dPmk_dthk    = sym.simplify(sym.diff(Pmk, thk))
dPmk_dthm    = sym.simplify(sym.diff(Pmk, thm))
dPmk_dx_tcsc = sym.simplify(sym.diff(Pmk, x_tcsc))

# ---------------------------------------------------------------------------
# Reactive power flow Qmk
# ---------------------------------------------------------------------------
dQmk_dVk     = sym.simplify(sym.diff(Qmk, Vk))
dQmk_dVm     = sym.simplify(sym.diff(Qmk, Vm))
dQmk_dthk    = sym.simplify(sym.diff(Qmk, thk))
dQmk_dthm    = sym.simplify(sym.diff(Qmk, thm))
dQmk_dx_tcsc = sym.simplify(sym.diff(Qmk, x_tcsc))

# ---------------------------------------------------------------------------
# Real current Ikm
# ---------------------------------------------------------------------------
dIkm_re_dVk     = sym.simplify(sym.diff(Ikm_re, Vk))
dIkm_re_dVm     = sym.simplify(sym.diff(Ikm_re, Vm))
dIkm_re_dthk    = sym.simplify(sym.diff(Ikm_re, thk))
dIkm_re_dthm    = sym.simplify(sym.diff(Ikm_re, thm))
dIkm_re_dx_tcsc = sym.simplify(sym.diff(Ikm_re, x_tcsc))

# ---------------------------------------------------------------------------
# Imaginary current Ikm
# ---------------------------------------------------------------------------
dIkm_im_dVk     = sym.simplify(sym.diff(Ikm_im, Vk))
dIkm_im_dVm     = sym.simplify(sym.diff(Ikm_im, Vm))
dIkm_im_dthk    = sym.simplify(sym.diff(Ikm_im, thk))
dIkm_im_dthm    = sym.simplify(sym.diff(Ikm_im, thm))
dIkm_im_dx_tcsc = sym.simplify(sym.diff(Ikm_im, x_tcsc))

# ---------------------------------------------------------------------------
# Real current Imk
# ---------------------------------------------------------------------------
dImk_re_dVk     = sym.simplify(sym.diff(Imk_re, Vk))
dImk_re_dVm     = sym.simplify(sym.diff(Imk_re, Vm))
dImk_re_dthk    = sym.simplify(sym.diff(Imk_re, thk))
dImk_re_dthm    = sym.simplify(sym.diff(Imk_re, thm))
dImk_re_dx_tcsc = sym.simplify(sym.diff(Imk_re, x_tcsc))

# ---------------------------------------------------------------------------
# Imaginary current Imk
# ---------------------------------------------------------------------------
dImk_im_dVk     = sym.simplify(sym.diff(Imk_im, Vk))
dImk_im_dVm     = sym.simplify(sym.diff(Imk_im, Vm))
dImk_im_dthk    = sym.simplify(sym.diff(Imk_im, thk))
dImk_im_dthm    = sym.simplify(sym.diff(Imk_im, thm))
dImk_im_dx_tcsc = sym.simplify(sym.diff(Imk_im, x_tcsc))

print("Latex code")
print("="*80)

print("Julia code")
print("="*80)

# ============================================================================
# Active power flow Pkm
# ============================================================================
print("Active power Pkm")
print("Vk")
print(print_latex(dPkm_dVk))

print("Vm")
print(print_latex(dPkm_dVm))

print("θk")
print(print_latex(dPkm_dthk))

print("θm")
print(print_latex(dPkm_dthm))

print("x_tcsc")
print(print_latex(dPkm_dx_tcsc))


# ============================================================================
# Reactive power flow Qkm
# ============================================================================
print("\nReactive power Qkm")
print("Vk")
print(print_latex(dQkm_dVk))

print("Vm")
print(print_latex(dQkm_dVm))

print("θk")
print(print_latex(dQkm_dthk))

print("θm")
print(print_latex(dQkm_dthm))

print("x_tcsc")
print(print_latex(dQkm_dx_tcsc))


# ============================================================================
# Active power flow Pmk
# ============================================================================
print("\nActive power Pmk")
print("Vk")
print(print_latex(dPmk_dVk))

print("Vm")
print(print_latex(dPmk_dVm))

print("θk")
print(print_latex(dPmk_dthk))

print("θm")
print(print_latex(dPmk_dthm))

print("x_tcsc")
print(print_latex(dPmk_dx_tcsc))


# ============================================================================
# Reactive power flow Qmk
# ============================================================================
print("\nReactive power Qmk")
print("Vk")
print(print_latex(dQmk_dVk))

print("Vm")
print(print_latex(dQmk_dVm))

print("θk")
print(print_latex(dQmk_dthk))

print("θm")
print(print_latex(dQmk_dthm))

print("x_tcsc")
print(print_latex(dQmk_dx_tcsc))


# ============================================================================
# Current Ikm (real)
# ============================================================================
print("\nCurrent Ikm (real)")
print("Vk")
print(print_latex(dIkm_re_dVk))

print("Vm")
print(print_latex(dIkm_re_dVm))

print("θk")
print(print_latex(dIkm_re_dthk))

print("θm")
print(print_latex(dIkm_re_dthm))

print("x_tcsc")
print(print_latex(dIkm_re_dx_tcsc))


# ============================================================================
# Current Ikm (imaginary)
# ============================================================================
print("\nCurrent Ikm (imaginary)")
print("Vk")
print(print_latex(dIkm_im_dVk))

print("Vm")
print(print_latex(dIkm_im_dVm))

print("θk")
print(print_latex(dIkm_im_dthk))

print("θm")
print(print_latex(dIkm_im_dthm))

print("x_tcsc")
print(print_latex(dIkm_im_dx_tcsc))


# ============================================================================
# Current Imk (real)
# ============================================================================
print("\nCurrent Imk (real)")
print("Vk")
print(print_latex(dImk_re_dVk))

print("Vm")
print(print_latex(dImk_re_dVm))

print("θk")
print(print_latex(dImk_re_dthk))

print("θm")
print(print_latex(dImk_re_dthm))

print("x_tcsc")
print(print_latex(dImk_re_dx_tcsc))


# ============================================================================
# Current Imk (imaginary)
# ============================================================================
print("\nCurrent Imk (imaginary)")
print("Vk")
print(print_latex(dImk_im_dVk))

print("Vm")
print(print_latex(dImk_im_dVm))

print("θk")
print(print_latex(dImk_im_dthk))

print("θm")
print(print_latex(dImk_im_dthm))

print("x_tcsc")
print(print_latex(dImk_im_dx_tcsc))



#%%
print("Julia code")
print("="*80)

# ============================================================================
# Active power flow Pkm
# ============================================================================
print("Active power Pkm")
print("Vk")
print(print_julia(dPkm_dVk,mapping))

print("Vm")
print(print_julia(dPkm_dVm,mapping))

print("θk")
print(print_julia(dPkm_dthk,mapping))

print("θm")
print(print_julia(dPkm_dthm,mapping))

print("x_tcsc")
print(print_julia(dPkm_dx_tcsc,mapping))


# ============================================================================
# Reactive power flow Qkm
# ============================================================================
print("\nReactive power Qkm")
print("Vk")
print(print_julia(dQkm_dVk,mapping))

print("Vm")
print(print_julia(dQkm_dVm,mapping))

print("θk")
print(print_julia(dQkm_dthk,mapping))

print("θm")
print(print_julia(dQkm_dthm,mapping))

print("x_tcsc")
print(print_julia(dQkm_dx_tcsc,mapping))


# ============================================================================
# Active power flow Pmk
# ============================================================================
print("\nActive power Pmk")
print("Vk")
print(print_julia(dPmk_dVk,mapping))

print("Vm")
print(print_julia(dPmk_dVm,mapping))

print("θk")
print(print_julia(dPmk_dthk,mapping))

print("θm")
print(print_julia(dPmk_dthm,mapping))

print("x_tcsc")
print(print_julia(dPmk_dx_tcsc,mapping))


# ============================================================================
# Reactive power flow Qmk
# ============================================================================
print("\nReactive power Qmk")
print("Vk")
print(print_julia(dQmk_dVk,mapping))

print("Vm")
print(print_julia(dQmk_dVm,mapping))

print("θk")
print(print_julia(dQmk_dthk,mapping))

print("θm")
print(print_julia(dQmk_dthm,mapping))

print("x_tcsc")
print(print_julia(dQmk_dx_tcsc,mapping))


# ============================================================================
# Current Ikm (real)
# ============================================================================
print("\nCurrent Ikm (real)")
print("Vk")
print(print_julia(dIkm_re_dVk,mapping))

print("Vm")
print(print_julia(dIkm_re_dVm,mapping))

print("θk")
print(print_julia(dIkm_re_dthk,mapping))

print("θm")
print(print_julia(dIkm_re_dthm,mapping))

print("x_tcsc")
print(print_julia(dIkm_re_dx_tcsc,mapping))


# ============================================================================
# Current Ikm (imaginary)
# ============================================================================
print("\nCurrent Ikm (imaginary)")
print("Vk")
print(print_julia(dIkm_im_dVk,mapping))

print("Vm")
print(print_julia(dIkm_im_dVm,mapping))

print("θk")
print(print_julia(dIkm_im_dthk,mapping))

print("θm")
print(print_julia(dIkm_im_dthm,mapping))

print("x_tcsc")
print(print_julia(dIkm_im_dx_tcsc,mapping))


# ============================================================================
# Current Imk (real)
# ============================================================================
print("\nCurrent Imk (real)")
print("Vk")
print(print_julia(dImk_re_dVk,mapping))

print("Vm")
print(print_julia(dImk_re_dVm,mapping))

print("θk")
print(print_julia(dImk_re_dthk,mapping))

print("θm")
print(print_julia(dImk_re_dthm,mapping))

print("x_tcsc")
print(print_julia(dImk_re_dx_tcsc,mapping))


# ============================================================================
# Current Imk (imaginary)
# ============================================================================
print("\nCurrent Imk (imaginary)")
print("Vk")
print(print_julia(dImk_im_dVk,mapping))

print("Vm")
print(print_julia(dImk_im_dVm,mapping))

print("θk")
print(print_julia(dImk_im_dthk,mapping))

print("θm")
print(print_julia(dImk_im_dthm,mapping))

print("x_tcsc")
print(print_julia(dImk_im_dx_tcsc,mapping))
# %%
