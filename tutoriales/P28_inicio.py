"""
P28 — Transporte de agua en la membrana y desbalance hídrico (punto de partida)
===============================================================================
Guía: docs/proyectos/P28_agua.md

El transporte de agua ya está en rfbmodel1D (membrane_residual): λ(x,t) con difusión D(λ)
y arrastre electro-osmótico ξ(λ)·i/F, más leyes de sorción en las dos caras.
Su aporte: ciclar, medir cuánta agua cruza y evaluar una estrategia de reposición.

Qué hace este script (ya funciona tal cual):
  - Corre UN ciclo (carga + descarga, 30 min cada tramo, ±500 A/m2).
  - En cada paso de tiempo separa el flujo de agua en difusión y arrastre
    (flujos_agua_membrana en iiq3843_utils.py) y registra el flujo en la cara del cátodo.
  - Verifica el balance de agua de la membrana y acumula el agua neta que gana el catolito.
Lo que usted programa (TODO): n ciclos, la estrategia de reposición y el barrido de corriente.

    python tutoriales/P28_inicio.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from iiq3843_utils import (configurar, construir_celda, integrar, cambiar_operacion,
                           flujos_agua_membrana, area_celda)

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)
J, T_TRAMO = 500.0, 1800.0
V_W = 1.809e-5     # volumen molar del agua líquida [m3/mol] (MembraneParams.V_w)

configurar(T_C=32.0, NC=100, NM=60)
cell, y, yp = construir_celda(0.5, "galvanostatico", +J)
A = area_celda(cell)

c_f, dx_m = cell.membrane.c_f, cell._dx_m
reg = {k: [] for k in ("t", "j", "N_cat", "N_an", "N_dif_medio", "N_eod_medio",
                       "N_cara1", "N_caraN", "lam0", "W_int")}
t_ini = 0.0
for j in (+J, -J):                                  # carga y luego descarga
    if j < 0:
        y, yp = cambiar_operacion(cell, y, "galvanostatico", j)
    sol = integrar(cell, y, yp, T_TRAMO)            # pasos del solver: NO use n_salidas aquí
    for k in range(sol.t.size):
        f = flujos_agua_membrana(cell, sol.y[:, k], sol.t[k])
        lam = cell.split_state(sol.y[:, k])["membrane"][:, 0]
        N = f["N_dif"] + f["N_eod"]                 # flujo total en caras interiores
        reg["t"].append(t_ini + sol.t[k]); reg["j"].append(j)
        reg["N_cat"].append(f["N_cat"]); reg["N_an"].append(f["N_an"])
        reg["N_dif_medio"].append(f["N_dif"].mean()); reg["N_eod_medio"].append(f["N_eod"].mean())
        reg["N_cara1"].append(N[0]); reg["N_caraN"].append(N[-1])
        reg["lam0"].append(lam[0]); reg["W_int"].append(c_f * np.sum(lam[1:-1] * dx_m[1:-1]))
    t_ini += sol.t[-1]
    y = sol.y[:, -1]
r = {k: np.array(v) for k, v in reg.items()}

# Balance de agua de los volúmenes interiores de la membrana:
#   d/dt Σ c_f λ_k Δx_k = N(cara 1) − N(última cara)
# Se usan caras interiores porque el flujo en los bordes (leyes de sorción) es muy rígido.
for signo, nombre in ((+J, "carga"), (-J, "descarga")):
    m = r["j"] == signo
    dW   = r["W_int"][m][-1] - r["W_int"][m][0]
    int_ = np.trapezoid(r["N_cara1"][m] - r["N_caraN"][m], r["t"][m])
    print(f"{nombre}: ΔW interior = {dW:.4e} mol/m2, ∫(N_1 − N_N)dt = {int_:.4e} mol/m2, "
          f"error = {abs(dW - int_) / abs(dW) * 100:.2f} %")

# Agua que gana el catolito = −(agua que entra a la membrana por la cara del cátodo):
#   ∫N_cat dt = ∫N(cara 1) dt + c_f Δx_0 [λ_0(fin) − λ_0(inicio)]   (balance del volumen 0)
in_mem = np.trapezoid(r["N_cara1"], r["t"]) + c_f * dx_m[0] * (r["lam0"][-1] - r["lam0"][0])
n_cat = -in_mem * A                                  # [mol]
print(f"Agua neta ganada por el catolito en 1 ciclo: {n_cat:.3e} mol = {n_cat*V_W*1e6:.3f} mL "
      f"(área {A*1e4:.0f} cm2)")

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
tm = r["t"] / 60
ax[0].plot(tm, r["N_dif_medio"] * 1e3, label="difusión (media en caras)")
ax[0].plot(tm, r["N_eod_medio"] * 1e3, label="arrastre (media en caras)")
ax[0].plot(tm, r["N_cat"] * 1e3, "k--", label="cara cátodo")
ax[0].set(xlabel="t [min]", ylabel="flujo de agua [mmol/m2/s], +x hacia el ánodo"); ax[0].legend()
ax[1].plot(tm, r["W_int"]); ax[1].set(xlabel="t [min]", ylabel="agua en la membrana (interior) [mol/m2]")
fig.tight_layout(); fig.savefig(SALIDA / "P28_ciclo_agua.png", dpi=110)

# ---------------------------------------------------------------------------------------
# TODO (P28):
#  1. Repita el ciclo n veces (copie el lazo de arriba dentro de un for) y acumule n_cat.
#  2. src/estrategias_agua.py: con V_cat(t) = V0 + V_W·(agua ganada), aplique reposición
#     cuando |V − V0|/V0 > θ y cuente el volumen movido. El modelo mantiene V_tank constante
#     (su derivada es 0 en tank_residual): la estrategia se evalúa a posteriori.
#  3. Barra J = 250, 500, 1000 A/m2 y compare el agua por ciclo.
# ---------------------------------------------------------------------------------------
