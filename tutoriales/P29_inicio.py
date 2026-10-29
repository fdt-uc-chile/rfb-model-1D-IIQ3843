"""
P29 — Balance térmico y riesgo de precipitación de V2O5 (punto de partida)
==========================================================================
Guía: docs/proyectos/P29_termico.md

rfbmodel1D es isotérmico: la temperatura es un parámetro (CellParams.T) con el que se
calculan viscosidad, velocidad, presión de saturación y cinética al construir la celda.
Por eso el alcance mínimo acopla en UNA dirección: el modelo entrega el calor generado
Q_gen(t) a T fija y usted integra aparte un balance de energía concentrado del electrolito.

Qué hace este script (ya funciona tal cual):
  - Corre descarga y carga (30 min cada una, ±500 A/m2) y calcula
        Q_gen = i · A · (V − OCV)   [W]   (positivo en carga y en descarga)
  - Escala a un stack de N_CELDAS celdas y llama a balance_energia(), hoy vacío.
Lo que usted programa (TODO): balance_energia() y T_amb(t).

    python tutoriales/P29_inicio.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from iiq3843_utils import (configurar, construir_celda, integrar, cambiar_operacion,
                           series_temporales, area_celda)

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)
J, T_TRAMO = 500.0, 1800.0
N_CELDAS = 40            # celdas del stack (supuesto: declárelo en el informe)

# ---- Calor generado por la celda ---------------------------------------------------
configurar(T_C=32.0, NC=100, NM=60)
cell, y, yp = construir_celda(0.75, "galvanostatico", -J)
A = area_celda(cell)
t_all, Q_all, t_ini = [], [], 0.0
for j in (-J, +J):
    if j > 0:
        y, yp = cambiar_operacion(cell, y, "galvanostatico", j)
    sol = integrar(cell, y, yp, T_TRAMO)
    ts = series_temporales(cell, sol)
    Q = ts["i_io"] * A * (ts["V"] - ts["OCV"])      # [W] por celda
    t_all.append(t_ini + ts["t"]); Q_all.append(Q)
    print(f"j = {j:+.0f} A/m2: V {ts['V'][0]:.4f}→{ts['V'][-1]:.4f} V, OCV {ts['OCV'][0]:.4f}→{ts['OCV'][-1]:.4f} V, "
          f"Q_gen {Q[1]:.4f}→{Q[-1]:.4f} W/celda")
    t_ini += sol.t[-1]; y = sol.y[:, -1]
t_q = np.concatenate(t_all); Q_q = np.concatenate(Q_all) * N_CELDAS
print(f"Energía disipada en el ciclo: {np.trapezoid(Q_q, t_q)/3.6e6*1e3:.2f} Wh (stack de {N_CELDAS} celdas)")


# ---- Balance de energía concentrado ------------------------------------------------
def T_amb(t):
    """TODO (P29): temperatura ambiente [K] en el instante t [s] (datos SINCA o sinusoide)."""
    return 273.15 + 25.0


def balance_energia(t, T, m_cp, hA):
    """TODO (P29): dT/dt [K/s] del electrolito.
        m_cp dT/dt = Q_gen(t) − hA (T − T_amb(t)) − ε σ A_ext (T⁴ − T_cielo⁴)
    Q_gen(t) se interpola de (t_q, Q_q) y hoy es solo el calor irreversible; agregue el reversible
    Q_rev = i·A·T·ΔS_r/F (CellParams.dS_r), ver la guía. Hoy devuelve 0: la temperatura no cambia.
    """
    Q = np.interp(t, t_q, Q_q)
    return [0.0]


m_cp = 0.5e-3 * N_CELDAS * 1350 * 3000      # [J/K]: tanque de 0,5 L por celda (supuesto)
hA   = 2.0                                  # [W/K] (supuesto: barrer en el proyecto)
sol_T = solve_ivp(balance_energia, (t_q[0], t_q[-1]), [T_amb(0.0)], args=(m_cp, hA), max_step=30.0)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(t_q / 60, Q_q); ax[0].set(xlabel="t [min]", ylabel="Q_gen stack [W]")
ax[1].plot(sol_T.t / 60, sol_T.y[0] - 273.15); ax[1].axhline(40, color="r", ls="--", label="40 °C")
ax[1].set(xlabel="t [min]", ylabel="T electrolito [°C]"); ax[1].legend()
fig.tight_layout(); fig.savefig(SALIDA / "P29_calor.png", dpi=110)
print("Listo: resultados/P29_calor.png. Ahora implemente balance_energia() y T_amb().")
