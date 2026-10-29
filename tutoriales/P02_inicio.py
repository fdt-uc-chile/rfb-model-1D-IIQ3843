"""
P02 — Pseudo-2D a lo largo del canal (punto de partida)
=======================================================
Guía: docs/proyectos/P02_pseudo2D.md

Idea: el canal se divide en Ny segmentos en serie. Cada segmento es el modelo 1D completo
(a través del espesor) alimentado con una composición de entrada fija; su salida es la
entrada del segmento siguiente. En modo potenciostático todos los segmentos comparten el
voltaje de celda (colectores conductores), así que no hace falta iterar sobre V.

Qué hace este script (ya funciona tal cual):
  - correr_segmento(): corre UN segmento hasta estado cuasi-estacionario y devuelve su
    corriente, la composición de salida mezclada y los términos del balance de especie.
  - Verifica Faraday para Ny = 1.
Lo que usted programa (TODO al final): el lazo sobre los Ny segmentos.

    python tutoriales/P02_inicio.py
"""
from pathlib import Path
import numpy as np

from iiq3843_utils import (configurar, construir_celda, integrar, cambiar_operacion,
                           perfiles, area_celda)

SALIDA = Path("resultados"); SALIDA.mkdir(exist_ok=True)


def correr_segmento(c_in, V, Ny, Q_ml_min=45.0, t_final=200.0, NC=100, NM=60):
    """Un segmento j de Ny, con entrada c_in = [c_V4, c_V5, c_H, c_HSO4] [mol/m3], voltaje V [V]
    y caudal total de catolito Q_ml_min [ml/min].

    Truco: el 'tanque' del modelo se hace enorme (1 m3) para que su composición no cambie;
    así actúa como la corriente de entrada al segmento. La altura del fieltro se divide por Ny
    (el área del segmento es A/Ny) y se conserva la velocidad en el plano u_ip de la celda
    completa, calculada en Cell.__init__ por solve_velocity().

    Devuelve un dict con:
      i        densidad de corriente del segmento [A/m2] (<0 descarga)
      A_seg    área del segmento [m2]
      c_fieltro concentración media en el fieltro = salida del caudal que lo atraviesa
      c_mezcla composición de salida del segmento (fieltro + canal, mezcladas)
      balance  términos del balance de V(V) [mol/s]
    """
    configurar(T_C=32.0, Q_ml_min=Q_ml_min, NC=NC, NM=NM)
    soc = c_in[1] / (c_in[0] + c_in[1])
    cell, y0, yp0 = construir_celda(soc, "potenciostatico", V, V_tank=1.0,
                                    c_V=c_in[0] + c_in[1], c_H=c_in[2])
    cell.cathode.height = cell.cathode.height / Ny          # segmento de altura H/Ny
    y0, yp0 = cambiar_operacion(cell, y0, "potenciostatico", V)   # re-consistencia
    sol = integrar(cell, y0, yp0, t_final)
    y = sol.y[:, -1]
    p = perfiles(cell, y)
    w = p["dx_ca"] / p["dx_ca"].sum()
    c_fieltro = np.array([w @ p["c_V4"], w @ p["c_V5"], w @ p["c_H"], w @ p["c_HSO4"]])
    A_seg = area_celda(cell)
    Q, Q_ca = cell.cathode.Q_total, cell.cathode.Q_ca
    Q_ch = Q - Q_ca                                          # catolito que no entra al fieltro
    # Flujo difusivo en la cara x = 0 (fieltro/colector-canal), igual que update_BC_tank_cathode:
    D = np.array([cell.cathode.D_V4_eff, cell.cathode.D_V5_eff,
                  cell.cathode.D_H_plus_eff, cell.cathode.D_HSO4_minus_eff])
    c0 = np.array([p["c_V4"][0], p["c_V5"][0], p["c_H"][0], p["c_HSO4"][0]])
    N0 = -D * (c0 - c_in) / (p["dx_ca"][0] / 2)              # [mol/m2/s], >0 hacia el fieltro
    # Salida mezclada: el fieltro devuelve Q_ca a c_fieltro; el canal entrega Q_ch a c_in,
    # menos lo que cedió al fieltro por difusión en x = 0.
    c_mezcla = (Q_ca * c_fieltro + Q_ch * c_in - A_seg * N0) / Q
    consumo_V5 = Q * (c_in[1] - c_mezcla[1])                 # [mol/s] V(V) que desaparece
    faraday_V5 = -p["i_io"] * A_seg / cell.F                 # [mol/s] según la corriente
    return {"i": p["i_io"], "A_seg": A_seg, "c_fieltro": c_fieltro, "c_mezcla": c_mezcla,
            "Q": Q, "Q_ca": Q_ca, "fraccion_sin_penetrar": Q_ch / Q,
            "balance": {"convectivo_fieltro": Q_ca * (c_in[1] - c_fieltro[1]),
                        "difusivo_x0": A_seg * N0[1], "total": consumo_V5, "faraday": faraday_V5}}


if __name__ == "__main__":
    c_in = np.array([500.0, 500.0, 4963.0, 4963.0 * 3 / 5])   # SOC 0,5 a la entrada del canal
    V = 1.05                                                  # voltaje de celda [V]
    r = correr_segmento(c_in, V, Ny=1)
    b = r["balance"]
    err = abs(b["total"] - b["faraday"]) / abs(b["faraday"])
    print(f"Ny = 1: i = {r['i']:.1f} A/m2, fracción que no penetra el fieltro = {r['fraccion_sin_penetrar']:.3f}")
    print(f"  consumo V(V): fieltro {b['convectivo_fieltro']:.3e} + cara x=0 {b['difusivo_x0']:.3e} "
          f"= {b['total']:.3e} mol/s;  Faraday {b['faraday']:.3e} mol/s;  error {err*100:.2f} %")
    print(f"  c_V5: entrada {c_in[1]:.1f}, fieltro {r['c_fieltro'][1]:.1f}, salida mezclada {r['c_mezcla'][1]:.1f} mol/m3")

    # ---------------------------------------------------------------------------------
    # TODO (P02): lazo pseudo-2D.
    #   Ny = 4
    #   c = c_in.copy(); I_total = 0; perfiles_y = []
    #   for j in range(Ny):
    #       r = correr_segmento(c, V, Ny)
    #       I_total += r["i"] * r["A_seg"]
    #       perfiles_y.append((j, r["i"], r["c_mezcla"][1]))
    #       c = r["c_mezcla"]                      # la salida del segmento j entra al j+1
    #   Compare I_total con r["i"]*A de Ny = 1 y repita para Ny = 1, 2, 4, 8 (convergencia).
    # ---------------------------------------------------------------------------------
