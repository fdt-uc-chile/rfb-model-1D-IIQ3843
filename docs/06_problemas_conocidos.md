# 6. Problemas conocidos y correcciones

Revisado el 29-10-2026 sobre rfbmodel-1D v1 (Zenodo 10.5281/zenodo.22274872).

## 6.1 Corregido en esta copia

| Problema | Efecto en v1 | Corrección | Commit |
| --- | --- | --- | --- |
| La matriz `Div` del potencial electrónico del cátodo usaba anchos de celda en las filas interiores (`submesh.build_diffusion_non_uniform_operator`) y distancias entre centros en las filas de borde (`Cell.apply_electronic_*_BC`). La primera y la última cara no conservaban carga. | La corriente en la membrana y ∫j dx excedían la corriente aplicada en 1,5 % a 500 A/m² (−507,4 contra −500 A/m²). El inventario de V(V) se desviaba de Faraday en la misma proporción. | Distancias entre centros en todas las filas. Ahora la corriente cierra con error menor a 10⁻⁵ en NC = 100 y 250. El voltaje cambia 0,5 mV. | "Corrige conservación de carga…" |
| `build_diffusion_operator` (malla uniforme) no tenía el factor 1/Δx. | Solo afecta si se usa `MeshParams.Uniform_C = True` (no es el caso por defecto). | Se multiplica por n_x. | ídem |
| `scripts/Validation/run_simulation.py` armaba `atol` sin la incógnita de acoplamiento. | Vector de tolerancias una entrada más corto que el estado. | Se agrega la entrada. | ídem |

Para volver al comportamiento exacto de v1 (por ejemplo, para reproducir una figura del
artículo): `git revert <commit>` o compare con el primer commit del historial.

## 6.2 Limitaciones del modelo que debe declarar en su informe

| Limitación | Dónde | Consecuencia para los proyectos |
| --- | --- | --- |
| Isotérmico: T es un parámetro fijo al construir la celda | `CellParams.T`, usada en `Cell.__init__` | P29 acopla el balance de energía en una dirección |
| Volumen del tanque constante (dV/dt = 0) | `tank_residual` | P28 contabiliza el desbalance de agua a posteriori |
| Sin cruce de vanadio por la membrana (flujo cero en la cara) | condiciones Neumann de c_V4, c_V5 en `cathode_residual` | P01 no puede cerrar un balance de vanadio membrana–electrolito |
| Ánodo concentrado (un volumen por capa), HOR lineal | `anode_residual`, `BV_lumped` | P03 construye un modelo 1D del ánodo aparte |
| Composición uniforme en la dirección del flujo (el fieltro es un reactor de mezcla en el plano) | término u_ip/H (c_T − c) | P02 lo resuelve con segmentos en serie |

## 6.3 Advertencias de uso

- **Opción `CellParams.vanadium_uptake = True`.** Falla al inicializar (`Algebraic solve
  failed`): usa ΔΦ_io de la membrana como denominador y en el estado inicial ese valor es 0.
  Además toma la corriente del valor de la condición de Φ_el, que solo es una corriente en modo
  galvanostático. No la use; P01 implementa su propia corrección.
- **`compute_OCP` fija `F_y = 10.02`** e ignora `CathodeParams.F_y` (25,0). Si su proyecto
  depende del OCP, redefina `compute_OCP` en una subclase y declare qué valor usa.
- **Estado compartido.** El diccionario de condiciones de borde de `params.py` y los atributos de
  `SimulationParams` son objetos de módulo que todas las celdas comparten. Una segunda celda en
  la misma sesión hereda la configuración de la primera. Llame a `configurar(...)` y fije el modo
  de control antes de cada `construir_celda`, o use un proceso por caso.
- **Salidas interpoladas (`t_eval`).** La interpolación de `solve_dae` entre pasos no respeta las
  relaciones algebraicas del DAE. En el balance de agua de la membrana introduce 7,5 % de error,
  contra 0,01 % usando los pasos propios del solver. Para integrar flujos en el tiempo, use
  `sol.t` sin `n_salidas`.
- **Flujos en los bordes de la membrana.** Las leyes de sorción son muy rígidas: un error pequeño
  en λ₀ produce un error grande en el flujo de borde. Para balances, use la primera cara
  interior más el cambio de inventario del primer volumen (así lo hace `tutoriales/P28_inicio.py`).
- **`RuntimeWarning: invalid value encountered in scalar power`** durante la inicialización:
  `fsolve` prueba estados no físicos al buscar las condiciones iniciales. Es inofensiva si la
  corrida termina con `éxito=True`.
- **`results/` y `notebooks/`** contienen los ajustes y figuras del artículo, generados con v1
  (antes de la corrección 6.1). Las diferencias son del orden de 0,5 mV en el voltaje.

## 6.4 Pendiente

- La tesis de J. Barros (citada en las guías como material de apoyo) no está en este repositorio.
- El artículo asociado (`HVRFB_paper.pdf`) se distribuye por los materiales del curso, no aquí.
