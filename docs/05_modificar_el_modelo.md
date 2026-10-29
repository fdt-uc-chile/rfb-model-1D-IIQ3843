# 5. Cómo modificar el modelo

Primero decida qué tipo de cambio necesita; el trabajo crece de A a D.

| Tipo | Ejemplo | Cómo | ¿Cambia el vector de estado? |
| --- | --- | --- | --- |
| A. Ley constitutiva | ξ(λ), κ_m(λ), D(λ), OCP, cinética | Subclase de `Cell` que redefine el método | No |
| B. Término fuente o sumidero | reacción secundaria, pérdida por cruce | Editar `*_residual` en `dae_residual_full_cell.py` | No |
| C. Condición de borde o de interfaz | otra ley de sorción, cruce de vanadio | Editar `update_BC_*` y el flujo de borde en `*_residual` | No |
| D. Incógnita nueva | temperatura, una especie más | Todo lo anterior + diseño del estado, escalas, inicialización, tolerancias | **Sí** |
| E. Capa externa | segmentos en serie, ciclado, balance de energía aparte | Script que llama a `Cell` varias veces | No |

Antes de editar, cree una rama: `git checkout -b P0X-apellido`. Así siempre puede volver a la base.

## 5.1 Tipo A: subclase de `Cell` (recomendado)

Redefina el método que el residuo llama mediante `self` (lista en
[03_arquitectura.md §3.4](03_arquitectura.md#34-mapa-fenómeno--código)). Ejemplo: cambiar el
coeficiente de arrastre electro-osmótico.

```python
import numpy as np
from rfbmodel1D.cells.cell import Cell
from iiq3843_utils import configurar, construir_celda, integrar, perfiles

class CeldaArrastre(Cell):
    def EOD_coefficient(self, lam):
        lam = np.clip(np.asarray(lam, dtype=float), 0, 22)
        return 1.0 * lam / 22          # ξ(λ) nuevo, sin dependencia de T

configurar(NC=100, NM=60)
cell, y0, yp0 = construir_celda(0.75, "galvanostatico", -500.0, clase=CeldaArrastre)
sol = integrar(cell, y0, yp0, 600.0)
```

Reglas para que funcione:

- **Mantenga la firma y la forma de lo que devuelve.** Las leyes se llaman con arreglos en el
  interior y con escalares en las interfaces: use `np.asarray` y operaciones de NumPy, nunca un
  `if` de Python sobre un arreglo.
- **Mantenga la ley suave.** El método de Newton del integrador necesita un jacobiano continuo;
  por eso el código usa mezclas `tanh` y `logaddexp` en vez de `if` o `max`.
- **Llame a `super()`** si solo quiere corregir la ley original (así queda claro qué cambió).

## 5.2 Tipo B, ejemplo resuelto: pérdida de V(V) por cruce de hidrógeno

En una celda H₂–V, el H₂ que cruza la membrana reduce químicamente al V(V) en el catolito:
VO₂⁺ + H⁺ + ½ H₂ → VO²⁺ + H₂O. Se modela como reacción homogénea de primer orden
r_x = k_x c_V5 [mol m⁻³ s⁻¹] en el fieltro; la estequiometría da las fuentes nuevas:

```math
S_{V4} = -\frac{j_n}{F} + r_x,\qquad S_{V5} = \frac{j_n}{F} - r_x,\qquad S_{H} = \frac{2j_n}{F} + b_{dis} - r_x
```

**Paso 1 — parámetro**, en `cathodeparams` de `rfbmodel1D/parameters/params.py`, en cero para que
la base no cambie:

```python
    F_y:float        = 25.0                  # Fitting parameter for cathode OCV (activity coefficient)
    k_x: float       = 0.0                   # V(V) reduction by crossover H2, first order [1/s]
```

**Paso 2 — copiarlo al objeto del subdominio**, al final de `Cathode.__init__` en
`rfbmodel1D/cells/subdomains/cathode.py`:

```python
        # Chemical side reaction
        self.k_x = CathodeParams.k_x # V(V) reduction by crossover H2 [1/s]
```

**Paso 3 — tasa y fuentes** en `cathode_residual`, justo después de calcular `b_dis`
(línea 218), y agregar `r_x` a tres argumentos `source_term=`:

```python
    b_dis = -k_d * ((c_H - c_HSO4) / (c_H + c_HSO4 + 1e-12) - Beta)
    # V(V) reduction by crossover H2: VO2+ + H+ + 1/2 H2 -> VO2+ + H2O
    r_x = self.cathode.k_x * c_V5 # [mol/m3/s]
    ...
        source_term=-BV / F + r_x,              # V(IV)
    ...
        source_term=BV / F - r_x,)              # V(V)
    ...
        source_term=2 * BV / F + b_dis - r_x,)  # H+
```

`species_transport` ya multiplica `source_term` por `dx` y lo resta: el término entra
correctamente al balance integrado sin otros cambios. El vector de estado, las escalas, las
tolerancias y la inicialización no cambian porque no se agregó ninguna incógnita.

**Paso 4 — correr**: en el Tutorial 1, `SimulationParams.CathodeParams.k_x = 1e-3` antes de
`construir_celda`. Corra cada caso en un proceso de Python nuevo.

**Paso 5 — verificar**: inventario total de V(V) (tanque + poros del fieltro) tras 600 s a
−500 A/m² desde SOC 0,75, comparado con Faraday (resultados del 29-10-2026 con esta copia):

| Caso | max \|F(0, y₀, y′₀)\| | Δn_V5 modelo [mol] | Δn_V5 Faraday [mol] | SOC tanque a 600 s |
| --- | --- | --- | --- | --- |
| k_x = 0 (base) | 8,4 × 10⁻⁴ | −7,7732 × 10⁻³ | −7,7732 × 10⁻³ | 73,52 % |
| k_x = 10⁻³ s⁻¹ | 8,4 × 10⁻⁴ | −1,2647 × 10⁻² | −7,7732 × 10⁻³ | 72,59 % |

La reacción secundaria consumió 4,87 × 10⁻³ mol de V(V). La estimación a mano es
k_x · c̄_V5 · V_poros · t = 10⁻³ × 750 × 1,08 × 10⁻⁵ m³ × 600 s ≈ 4,9 × 10⁻³ mol: la
modificación hace lo que debe.

## 5.3 Tipo C: condiciones de borde e interfaz

Cambie el flujo que aporta la cara, no el valor por defecto de `params.py`. En la cara
cátodo/membrana, el lado cátodo está en `cathode_residual` (mediante
`cell.bc[...]["cathode"]["right"]`, que fija `update_BC_cathode_membrane`) y el lado membrana en
`membrane_residual`. Lo que sale de un subdominio debe entrar al otro con el mismo valor, o no se
conservan masa ni carga.

## 5.4 Tipo D: lista de chequeo para una incógnita nueva

El diseño del vector de estado está escrito en varios lugares. Si olvida uno, obtiene un error
de dimensiones o, peor, un índice silenciosamente equivocado.

1. `params.py` → `state_layout_params`: aumente `n_vars_<subdominio>`.
2. `*_residual`: lea la columna nueva de `view`, escriba su ecuación en `res_view[:, k]` y
   verifique que todas las columnas quedan llenas.
3. `dae_residual_full_cell`: agregue a la escala del bloque (`scale_cathode`, `scale_membrane`
   o `scale_anode`) una magnitud característica para la ecuación nueva.
4. `cells/methods/initialize_dae_state.py`:
    1. agregue un valor inicial en `y_initial_cathode`, `y_membrane` o `y_anode`;
    2. actualice los desplazamientos escritos a mano (`5`, `6 * Nc`, `2 * Nm`, `7`) en
       `solve_algebraic`, `solve_differential` y `solve_warm_state`;
    3. agregue el índice a `idx_alg` si su ecuación no tiene y′, o a `idx_diff` si la tiene.
5. `tutoriales/iiq3843_utils.py` → `tolerancias`: una entrada más por nodo, en orden.
6. `Cell.split_state`, `post_processing/profiles.py` y `iiq3843_utils.perfiles`: exponga la
   variable por su nombre.
7. Toda función de eventos que indexe `y` directamente (por ejemplo, el corte por λ en
   `scripts/Validation/run_simulation.py`, que supone 2 variables por nodo de membrana).

## 5.5 Cómo comprobar cualquier modificación

- **Interruptor en cero.** Con el parámetro nuevo en cero, los resultados deben coincidir con la
  base dentro de la tolerancia del integrador.
- **Arranque consistente.** Tras `initialize_dae_state`, `max(abs(F))` debe ser del orden de
  10⁻³ o menor. Un valor grande indica una ecuación mal escalada o un valor inicial lejano.
- **Conservación.** Inventarios contra Faraday (como en 5.2) y corrientes en toda la celda (como
  en el Tutorial 1). `scripts/Validation/check_conservation.py` hace ambas cosas para varias
  corrientes y SOC (corrida larga: 2 h simuladas por caso).
- **Malla.** Repita con NC y NM dobles, o use `scripts/Validation/grid_convergence_cathode.py`
  y `grid_convergence_membrane.py`.
- **Sentido físico.** Revise signos en un caso que pueda razonar: una pérdida debe bajar el SOC
  más rápido en descarga y frenar su subida en carga.
