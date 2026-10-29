"""
Tutorial 1 — Corrida base y perfiles en todos los subdominios
=============================================================
Descarga galvanostática de 10 min a -500 A/m2 desde SOC = 0,75, con mallas reducidas
(100 volúmenes en el cátodo y 60 en la membrana). Tarda unos 10 s en un computador portátil.

Ejecutar desde la raíz del repositorio:
    python tutoriales/01_corrida_base_perfiles.py

Salidas (carpeta resultados/):
    tutorial1_perfiles.png   nueve paneles: tanque, cátodo, membrana y ánodo
    tutorial1.npz            t, y de la simulación
"""
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from iiq3843_utils import (configurar, construir_celda, integrar, perfiles,
                           series_temporales, area_celda)
from rfbmodel1D.post_processing import calculate_current_density

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)

# 1. Punto de operación ---------------------------------------------------------
SOC0   = 0.75      # estado de carga inicial [-]
J_APPL = -500.0    # densidad de corriente aplicada [A/m2] (negativa = descarga)
T_FIN  = 600.0     # tiempo simulado [s]

# 2. Parámetros: SIEMPRE antes de construir la celda ---------------------------
configurar(T_C=32.0, Q_ml_min=45.0, NC=100, NM=60)

# 3. Celda + estado inicial consistente ----------------------------------------
cell, y0, yp0 = construir_celda(SOC0, "galvanostatico", J_APPL)
print("incógnitas:", y0.size)

# 4. Integración ---------------------------------------------------------------
t0 = time.time()
sol = integrar(cell, y0, yp0, T_FIN)
print(f"éxito={sol.success}  pasos={sol.t.size}  tiempo={time.time() - t0:.1f} s")
np.savez_compressed(SALIDA / "tutorial1.npz", t=sol.t, y=sol.y)

# 5. Post-proceso --------------------------------------------------------------
y_fin = sol.y[:, -1]
p   = perfiles(cell, y_fin, t=sol.t[-1])          # todos los perfiles al tiempo final
ts  = series_temporales(cell, sol)                # V, OCV, i_io, SOC vs. t
cur = calculate_current_density(cell, y_fin)      # corrientes en las caras
L_ca = cell.cathode.thickness
x_ca = p["x_ca"] * 1e3                            # [mm]
x_m  = (L_ca + p["x_m"]) * 1e3                    # [mm], membrana a continuación del cátodo

fig, ax = plt.subplots(3, 3, figsize=(15, 12))
a = ax[0, 0]; a.plot(ts["t"], ts["soc"] * 100)
a.set(title="Tanque: SOC", xlabel="t [s]", ylabel="SOC [%]")
a = ax[0, 1]; a.plot(x_ca, p["c_V4"], label="V(IV)"); a.plot(x_ca, p["c_V5"], label="V(V)")
a.set(title="Cátodo: vanadio", xlabel="x [mm]", ylabel="c [mol/m3]"); a.legend()
a = ax[0, 2]; a.plot(x_ca, p["c_H"], label="H+"); a.plot(x_ca, p["c_HSO4"], label="HSO4-")
a.set(title="Cátodo: ácido", xlabel="x [mm]", ylabel="c [mol/m3]"); a.legend()
a = ax[1, 0]; a.plot(x_ca, p["phi_el_ca"], label="Φ_el"); a.plot(x_ca, p["phi_io_ca"], label="Φ_io")
a.set(title="Cátodo: potenciales", xlabel="x [mm]", ylabel="V"); a.legend()
a = ax[1, 1]; a.plot(x_ca, p["eta_ca"] * 1e3)
a.set(title="Cátodo: sobrepotencial", xlabel="x [mm]", ylabel="η [mV]")
a = ax[1, 2]
a.plot(cur["x_ca"] * 1e3, cur["i_el_ca"], label="i_el cátodo")
a.plot(cur["x_ca"] * 1e3, cur["i_io_ca"], label="i_io cátodo")
a.plot(cur["x_mem"] * 1e3, cur["i_io_mem"], label="i_io membrana")
a.set(title="Densidades de corriente", xlabel="x [mm]", ylabel="A/m2"); a.legend()
a = ax[2, 0]; a.plot(p["x_m"] * 1e6, p["lam"])
a.set(title="Membrana: contenido de agua", xlabel="x en la membrana [µm]", ylabel="λ")
a = ax[2, 1]; a.plot(x_ca, p["phi_io_ca"], label="cátodo"); a.plot(x_m, p["phi_io_m"], label="membrana")
a.set(title="Φ_io en cátodo + membrana", xlabel="x [mm]", ylabel="Φ_io [V]"); a.legend()
a = ax[2, 2]; a.plot(ts["t"], ts["V"], label="V celda"); a.plot(ts["t"], ts["OCV"], "--", label="OCV")
a.set(title="Voltaje de celda", xlabel="t [s]", ylabel="V"); a.legend()
fig.tight_layout(); fig.savefig(SALIDA / "tutorial1_perfiles.png", dpi=110)

# 6. Controles rápidos ---------------------------------------------------------
an = p["anodo"]
print(f"V final = {ts['V'][-1]:.4f} V   OCV = {ts['OCV'][-1]:.4f} V")
print(f"i_io membrana = {p['i_io']:.2f} A/m2   ∫ j dx en el cátodo = {np.sum(p['j_ca'] * p['dx_ca']):.2f} A/m2")
print(f"λ membrana: {p['lam'].min():.2f} – {p['lam'].max():.2f}")
print(f"ánodo: η_act = {an['eta_acl']*1e3:.3f} mV, η_protón = {an['eta_proton']*1e3:.3f} mV, "
      f"s_ACL = {an['s_acl']:.4f}, c_H2 = {an['c_H2_acl']:.2f} mol/m3")
print(f"área = {area_celda(cell)*1e4:.0f} cm2;  figura en {SALIDA / 'tutorial1_perfiles.png'}")
plt.show()
