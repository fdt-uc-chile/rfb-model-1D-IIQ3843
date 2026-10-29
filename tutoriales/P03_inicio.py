"""
P03 — Resistencias de transporte en la ACL y la AGDL del ánodo (punto de partida)
=================================================================================
Guía: docs/proyectos/P03_anodo.md

El ánodo de rfbmodel1D está concentrado (un volumen para la ACL y otro para la AGDL),
pero ya incluye dos pérdidas que su proyecto debe separar y luego reproducir con un modelo 1D:
  η_act     = eta_acl, incógnita algebraica: j0 · 2Fη/(RT) + i_io = 0 (cinética lineal)
  η_protón  = cell.eta_proton, fórmula analítica de electrodo poroso con cinética lineal
              (update_BC_membrane_anode en rfbmodel1D/cells/DAE/dae_residual_full_cell.py)

Qué hace este script (ya funciona tal cual):
  Paso 1. Barre j0 (proxy de la carga de Pt) y el espesor de ACL con el modelo existente.
  Paso 2. Reproduce η_protón con la fórmula del código: es la referencia que su modelo 1D
          estacionario (src/anodo_1d.py) debe recuperar en el límite de cinética lineal.
  Paso 3. Muestra el enganche para una cinética no lineal (subclase de Cell, opcional).

    python tutoriales/P03_inicio.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from iiq3843_utils import configurar, construir_celda, integrar, perfiles, voltaje
from rfbmodel1D.cells.cell import Cell
from rfbmodel1D.parameters.params import SimulationParams, Vars

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)
J = -500.0


def eta_proton_analitico(cell, i_geo, eta_acl, c_H2):
    """Copia comentada de la fórmula de update_BC_membrane_anode (no la modifique: es la referencia).

    Electrodo poroso de espesor L con cinética lineal y conducción protónica σ_eff:
        c = sqrt(j0 n f / (σ_eff L)),  η_p = i / (σ_eff c) · [coth(cL) − 1/(cL)]
    """
    F, R, T = cell.F, cell.R, cell.T
    j0   = cell.anode.j_0 * 2                     # el código usa 2·j_0 aquí
    n, f = 2, F / (R * T)
    L    = cell.anode.thickness
    sig  = cell.membrane.kappa_io[-1] * cell.anode.ionomer_volume_fraction**1.5   # Bruggeman
    c    = (j0 * n * f / (sig * L))**0.5
    BV   = cell.BV_lumped(eta_acl, c_H2, "anode")  # = −i_io  [A/m2]
    return BV / (sig * c) * (1 / np.tanh(c * L) - 1 / (c * L))


def correr(clase=Cell):
    configurar(T_C=32.0, NC=100, NM=60)
    cell, y0, yp0 = construir_celda(0.75, "galvanostatico", J, clase=clase)
    sol = integrar(cell, y0, yp0, 300.0)
    y = sol.y[:, -1]
    p = perfiles(cell, y)
    an = p["anodo"]
    return cell, y, p, an


# ---- Paso 1: barrido de j0 y de espesor de ACL -------------------------------------
j0_base, L_base = SimulationParams.AnodeParams.j_0, SimulationParams.AnodeParams.thickness
filas = []
for fac_j0 in (1.0, 0.3, 0.1, 0.03, 0.01):
    for fac_L in (1.0, 3.0):
        SimulationParams.AnodeParams.j_0 = j0_base * fac_j0
        SimulationParams.AnodeParams.thickness = L_base * fac_L
        cell, y, p, an = correr()
        ref = eta_proton_analitico(cell, -p["i_io"], an["eta_acl"], an["c_H2_acl"])
        filas.append((fac_j0, fac_L, an["eta_acl"], an["eta_proton"], voltaje(cell, y)))
        print(f"j0 x{fac_j0:<5} L_ACL x{fac_L}: η_act = {an['eta_acl']*1e3:7.3f} mV  "
              f"η_protón = {an['eta_proton']*1e3:6.3f} mV (fórmula: {ref*1e3:6.3f})  V = {filas[-1][4]:.4f} V")
SimulationParams.AnodeParams.j_0, SimulationParams.AnodeParams.thickness = j0_base, L_base

f = np.array(filas)
fig, ax = plt.subplots(figsize=(6, 4))
for fac_L, mk in ((1.0, "o-"), (3.0, "s--")):
    s = f[f[:, 1] == fac_L]
    ax.loglog(s[:, 0], (s[:, 2] + s[:, 3]) * 1e3, mk, label=f"η_act + η_protón, L_ACL x{fac_L:g}")
ax.set(xlabel="j0 / j0_base (proxy de carga de Pt)", ylabel="pérdida anódica [mV]")
ax.legend(); fig.tight_layout(); fig.savefig(SALIDA / "P03_barrido_j0.png", dpi=110)
np.savetxt(SALIDA / "P03_barrido_j0.csv", f, delimiter=",",
           header="fac_j0,fac_L,eta_act_V,eta_proton_V,V_celda", comments="")


# ---- Paso 3 (opcional): cinética no lineal en el ánodo concentrado --------------------
class CeldaBVAnodo(Cell):
    """TODO (P03, ampliado): Butler–Volmer completo con dependencia de c_H2.

    El residuo llama a self.BV_lumped(eta, c_H2, "anode") en la ecuación de η y en η_protón.
    Hoy devuelve la cinética lineal original; reemplace el return por su expresión, p. ej.
        j0 * (c/c_ref * exp(αa F η/RT) − exp(−αc F η/RT))
    y verifique que para η pequeño coincide con la original (cuidado con n y los α).
    """
    def BV_lumped(self, eta, c, electrode):
        return super().BV_lumped(eta, c, electrode)


cell, y, p, an = correr(CeldaBVAnodo)
print(f"Enganche BV (sin cambios): η_act = {an['eta_acl']*1e3:.3f} mV (igual a la base).")
print("Listo: resultados/P03_barrido_j0.png y .csv. Siguiente paso: src/anodo_1d.py (ver la guía).")
