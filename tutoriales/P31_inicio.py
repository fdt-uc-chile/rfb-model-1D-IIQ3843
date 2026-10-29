"""
P31 — Distribución de corriente terciaria en el cátodo (punto de partida)
=========================================================================
Guía: docs/proyectos/P31_corriente_terciaria.md

Hallazgo clave: rfbmodel1D YA incluye la limitación por transferencia de masa. En
rfbmodel1D/electrochem/butler_volmer.py, cathode_surface_concentrations() calcula las
concentraciones en la superficie de las fibras con un coeficiente β0 = k_m · u_ip^0,4
(correlación de You et al., 2009). Su proyecto la expone, la reemplaza por k_m = D/δ y compara.

Qué hace este script (ya funciona tal cual):
  - CeldaTerciaria: subclase de Cell cuya cinética catódica llama a un método propio
    concentraciones_superficie(); con modo "codigo" reproduce exactamente la base.
  - Compara el modelo del código con el límite secundario (β0 → ∞, c_s = c_bulk).
  - Verifica ∫ j dx = i.
Lo que usted programa (TODO): el modo "delta" (k_m = D/δ) y los barridos.

    python tutoriales/P31_inicio.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from iiq3843_utils import configurar, construir_celda, integrar, perfiles, voltaje
from rfbmodel1D.cells.cell import Cell

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)


class CeldaTerciaria(Cell):
    modo  = "codigo"     # "codigo" (You et al.), "secundaria" (β0 → ∞) o "delta" (D/δ, TODO)
    delta = 50e-6        # espesor de capa límite [m] para el modo "delta"

    def coeficiente_beta0(self):
        """Coeficiente de transferencia de masa β0 [m/s] usado en las concentraciones de superficie."""
        if self.modo == "codigo":                      # idéntico a cathode_surface_concentrations
            k_m = 7 * self.cathode.D_V4_eff / self.d_f * (self.d_f / (self.cathode.mu_ca / self.cathode.rho_ca))**0.4
            return k_m * self.cathode.u_ip**0.4
        if self.modo == "secundaria":
            return 1e12                                # c_s = c_bulk: sin limitación de masa
        if self.modo == "delta":
            raise NotImplementedError("TODO (P31): devuelva D/δ (elija y justifique D)")
        raise ValueError(self.modo)

    def concentraciones_superficie(self, eta, c_V4, c_V5, c_H):
        """Copia de cathode_surface_concentrations con β0 reemplazable."""
        f, k_ca, a = self.F / (self.R * self.T), self.cathode.reaction_rate_constant, self.alpha_c
        beta_0 = self.coeficiente_beta0()
        eps = 1e-4
        c4 = np.sqrt(np.maximum(c_V4, 0)**2 + eps**2)
        c5 = np.sqrt(np.maximum(c_V5, 0)**2 + eps**2)
        cH = np.sqrt(np.maximum(c_H, 0)**2 + eps**2)
        P3 = (k_ca / beta_0) * (c5 / c4)**a * cH**(2 - 2 * a) * np.exp((1 - a) * f * eta)
        P4 = (k_ca / beta_0) * (c4 / c5)**(1 - a) * cH**(2 - 2 * a) * np.exp(-a * f * eta)
        cs4 = (P4 * c5 + (1 + P4) * c4) / (1 + P3 + P4)
        cs5 = (P3 * c4 + (1 + P3) * c5) / (1 + P3 + P4)
        return cs4, cs5, c4, c5

    def butler_volmer(self, eta, c, electrode):
        """Mismo Butler–Volmer volumétrico [A/m3] de electrochem/butler_volmer.py."""
        if electrode != "cathode":
            return super().butler_volmer(eta, c, electrode)
        F, f, a = self.F, self.F / (self.R * self.T), self.alpha_c
        NC = self.mesh.cathode_mesh.n_x
        cs4, cs5, c4, c5 = self.concentraciones_superficie(eta, c[:NC], c[NC:2*NC], c[2*NC:3*NC])
        cH = c[2*NC:3*NC]
        i0 = F * self.cathode.reaction_rate_constant * c5**(1 - a) * c4**a * cH**(2 - 2 * a)
        return self.cathode.specific_surface_area * i0 * (
            (cs4 / c4) * np.exp((1 - a) * f * eta) - (cs5 / c5) * np.exp(-a * f * eta))


J = -500.0
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for modo in ("codigo", "secundaria"):
    CeldaTerciaria.modo = modo
    configurar(T_C=32.0, NC=100, NM=60)
    cell, y0, yp0 = construir_celda(0.75, "galvanostatico", J, clase=CeldaTerciaria)
    sol = integrar(cell, y0, yp0, 300.0)
    y = sol.y[:, -1]
    p = perfiles(cell, y)
    beta0 = cell.coeficiente_beta0()
    j_lim = cell.cathode.specific_surface_area * cell.F * min(beta0, 1e6) * p["c_V5"]   # [A/m3], descarga
    print(f"modo {modo:10s}: V = {voltaje(cell, y):.5f} V  η(colector) = {p['eta_ca'][0]*1e3:7.3f} mV  "
          f"∫j dx = {np.sum(p['j_ca']*p['dx_ca']):8.2f} A/m2  β0 = {beta0:.3e} m/s  "
          f"max |j|/j_lim = {np.max(np.abs(p['j_ca']) / j_lim):.3f}")
    ax[0].plot(p["x_ca"] * 1e3, p["j_ca"] / 1e3, label=modo)
    ax[1].plot(p["x_ca"] * 1e3, p["eta_ca"] * 1e3, label=modo)
ax[0].set(xlabel="x [mm]", ylabel="j [kA/m3]"); ax[0].legend()
ax[1].set(xlabel="x [mm]", ylabel="η [mV]"); ax[1].legend()
fig.tight_layout(); fig.savefig(SALIDA / "P31_distribucion.png", dpi=110)
print("Listo: resultados/P31_distribucion.png. Ahora implemente el modo 'delta' y los barridos.")
