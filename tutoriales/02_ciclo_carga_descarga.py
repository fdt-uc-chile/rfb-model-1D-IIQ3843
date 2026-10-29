"""
Tutorial 2 — Encadenar carga y descarga (ciclos)
================================================
El modelo no tiene un "ciclador": cada tramo es una llamada a `integrar`. Para cambiar la
corriente se usa `cambiar_operacion`, que vuelve a hacer consistente el estado (y, yp).

Un ciclo de 30 min de carga + 30 min de descarga a ±500 A/m2 desde SOC = 0,5 (unos 20 s).

    python tutoriales/02_ciclo_carga_descarga.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from iiq3843_utils import (configurar, construir_celda, integrar, cambiar_operacion,
                           series_temporales)

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)
J, T_TRAMO, N_CICLOS = 500.0, 1800.0, 1

configurar(T_C=32.0, NC=100, NM=60)
cell, y, yp = construir_celda(0.5, "galvanostatico", +J)   # parte cargando

t_tot, V_tot, soc_tot, t_ini = [], [], [], 0.0
for ciclo in range(N_CICLOS):
    for j in (+J, -J):                                     # carga, luego descarga
        if not (ciclo == 0 and j > 0):
            y, yp = cambiar_operacion(cell, y, "galvanostatico", j)
        sol = integrar(cell, y, yp, T_TRAMO)
        ts  = series_temporales(cell, sol)
        t_tot.append(t_ini + ts["t"]); V_tot.append(ts["V"]); soc_tot.append(ts["soc"])
        t_ini += sol.t[-1]
        y = sol.y[:, -1]                                   # el tramo siguiente parte de aquí
        print(f"ciclo {ciclo+1}, j = {j:+.0f} A/m2: V {ts['V'][0]:.4f} → {ts['V'][-1]:.4f} V, "
              f"SOC {ts['soc'][0]:.4f} → {ts['soc'][-1]:.4f}")

t = np.concatenate(t_tot) / 60
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(t, np.concatenate(V_tot)); ax[0].set(xlabel="t [min]", ylabel="V celda [V]")
ax[1].plot(t, np.concatenate(soc_tot) * 100); ax[1].set(xlabel="t [min]", ylabel="SOC tanque [%]")
fig.tight_layout(); fig.savefig(SALIDA / "tutorial2_ciclo.png", dpi=110)
plt.show()
