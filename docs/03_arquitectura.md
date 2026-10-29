# 3. Arquitectura del código

Esta página responde a una pregunta: **¿en qué archivo y en qué función está implementado cada
fenómeno?** Los diagramas usan notación UML (clases y secuencia) en Mermaid, que GitHub dibuja
directamente.

## 3.1 Clases y archivos

`Cell` es el objeto central. Contiene un objeto por subdominio (propiedades), la malla y las
condiciones de borde. Casi toda la física está en **funciones de módulo** que reciben la celda
como primer argumento (`self`). `Cell` las expone como métodos de una línea; esa es la razón por
la que se pueden reemplazar con una subclase (sección 3.5).

```mermaid
classDiagram
    direction LR
    class SimulationParams {
        +cellparams CellParams
        +anodeparams AnodeParams
        +membraneparams MembraneParams
        +cathodeparams CathodeParams
        +agdlparams AGDL_Params
        +meshparams MeshParams
        +bcs BCParams
        +state_layout_params StateLayoutParams
    }
    class Cell {
        +Anode anode
        +Cathode cathode
        +Membrane membrane
        +AGDL agdl
        +fv_mesh mesh
        +dict bc
        +str bc_mode
        +float T, F, R
        +float i_io
        +float eta_proton
        +modify_BC(var, subdominio, lado, tipo, valor)
        +initialize_dae_state(c_tanque)
        +dae_residual_full_cell(t, y, yp)
        +split_state(y)
        +butler_volmer(eta, c, electrodo)
        +BV_lumped(eta, c, electrodo)
        +compute_OCP(c, electrodo)
        +compute_overpotential(...)
        +update_kappa_io_m(lam, Phi_io)
        +update_kappa_io_cathode(y_cat)
        +EOD_coefficient(lam)
        +membrane_water_diffusivity(lam, modelo)
        +return_lambda_eq(a)
        +capillary_pressure(s_r, dominio)
    }
    class Cathode {
        +thickness, width, height
        +porosity, specific_surface_area
        +reaction_rate_constant
        +D_V4_eff ... D_HSO4_minus_eff
        +mu_V4 ... mu_SO4_2minus
        +u_ip, u_tp, Q_ca
        +array kappa_io
        +array eta
    }
    class Membrane {
        +thickness, c_f, V_m, V_w
        +k_sorp_vw, k_desorp_vw
        +k_sorp_lw, k_desorp_lw, k_ca_lw
        +diffusivity_model
        +array kappa_io
    }
    class Anode {
        +thickness, porosity
        +ionomer_volume_fraction
        +j_0, permeability, s_im
        +k_cond, k_evap
    }
    class AGDL {
        +thickness, porosity
        +permeability, D_H2, D_vw
    }
    class fv_mesh {
        +submesh anode_mesh
        +submesh membrane_mesh
        +submesh cathode_mesh
    }
    class submesh {
        +int n_x
        +array x_nodes
        +array x_grid
        +array dx
        +array D
    }
    Cell ..> SimulationParams : lee al construirse
    Cell *-- Cathode
    Cell *-- Membrane
    Cell *-- Anode
    Cell *-- AGDL
    Cell *-- fv_mesh
    fv_mesh *-- submesh : 3
```

| Archivo | Contenido |
| --- | --- |
| `rfbmodel1D/cells/cell.py` | Clase `Cell`: construcción, condiciones de borde (`modify_BC`, `apply_*_BC`), operador de migración, `split_state` |
| `rfbmodel1D/cells/DAE/dae_residual_full_cell.py` | **Las ecuaciones de gobierno**: un residuo por subdominio y las actualizaciones de interfaz |
| `rfbmodel1D/cells/subdomains/` | `Cathode`, `Membrane`, `Anode`, `AGDL`: contenedores de propiedades |
| `rfbmodel1D/cells/methods/` | Leyes constitutivas y utilidades: κ, D(λ), ξ(λ), isotermas λ(a), p_sat, propiedades efectivas, estado inicial |
| `rfbmodel1D/electrochem/` | Butler–Volmer, OCP (Nernst), sobrepotencial, descomposición de pérdidas |
| `rfbmodel1D/meshes/grid_creation/` | `submesh` (malla de volúmenes finitos en [0, 1]) y `fv_mesh` |
| `rfbmodel1D/parameters/params.py` | Todos los parámetros, tamaños de malla, diseño del vector de estado y condiciones de borde por defecto |
| `rfbmodel1D/post_processing/` | Extracción de perfiles, densidades de corriente, gráficos simples |
| `tutoriales/iiq3843_utils.py` | (Curso) construir, integrar, extraer; no cambia la física |

## 3.2 Secuencia de una simulación

```mermaid
sequenceDiagram
    participant S as Script del estudiante
    participant C as Cell
    participant I as initialize_dae_state.py
    participant D as solve_dae (BDF)
    participant R as dae_residual_full_cell.py
    S->>C: Cell(params=SimulationParams)
    Note over C: D_eff, movilidades, viscosidad,<br/>velocidad en el fieltro (solve_velocity)
    S->>C: modify_BC(Vars.phi_el, "cathode", "left", tipo, valor)
    S->>C: initialize_dae_state(c_tanque)
    C->>I: estado uniforme igual al tanque
    I->>R: fsolve sobre las incógnitas algebraicas
    I->>R: fsolve sobre las derivadas y′
    I-->>S: y0, yp0 consistentes
    S->>D: solve_dae(cell.dae_residual_full_cell, (0, t_f), y0, yp0)
    loop cada paso de tiempo e iteración de Newton
        D->>R: F(t, y, y′)
        R-->>D: residuo escalado
    end
    D-->>S: sol.t, sol.y
```

## 3.3 Qué pasa dentro de una evaluación del residuo

```mermaid
flowchart TD
    A["dae_residual_full_cell(t, y, y′)"] --> B["1 · Desempaquetar bloques<br/>tanque | cátodo | membrana | ánodo | acoplamiento"]
    B --> C["2 · Propiedades<br/>update_vanadium<br/>self.update_kappa_io_cathode — κ(c) del catolito<br/>self.update_kappa_io_m — κ(λ) de la membrana"]
    C --> D["3 · Condiciones de interfaz (sobrescriben cell.bc)<br/>update_BC_tank_cathode — flujo en x = 0<br/>update_BC_cathode_membrane — i_io, Φ_io de cara, a_w<br/>update_BC_membrane_anode — Φ_io y actividad de agua del lado ánodo"]
    D --> E1["<b>tank_residual</b>"]
    D --> E2["<b>cathode_residual</b><br/>self.compute_overpotential<br/>→ self.compute_OCP<br/>self.butler_volmer<br/>→ cathode_surface_concentrations<br/>species_transport"]
    D --> E3["<b>membrane_residual</b><br/>self.membrane_water_diffusivity<br/>self.EOD_coefficient<br/>self.return_lambda_cath_eq<br/>self.return_lambda_eq"]
    D --> E4["<b>anode_residual</b><br/>self.BV_lumped<br/>self.capillary_pressure<br/>self.water_saturation_pressure"]
    D --> E5["<b>coupling_residual</b>"]
    E1 --> F["4 · Escalar cada bloque y concatenar → F"]
    E2 --> F
    E3 --> F
    E4 --> F
    E5 --> F
```

Tres consecuencias para quien modifica el modelo:

1. **Las condiciones de interfaz se recalculan en cada llamada.** Cambiar un valor por defecto
   en `params.py` no cambia una condición de interfaz: hay que cambiar la función `update_BC_*`.
2. **Las propiedades son efectos secundarios.** `cell.cathode.eta`, `cell.cathode.kappa_io`,
   `cell.membrane.kappa_io`, `cell.i_io` y `cell.eta_proton` guardan la última evaluación.
3. **Los residuos están escalados** por magnitudes características (100 A/m², 500 mol/m³ de
   vanadio, 3000 mol/m³ de H⁺). Una ecuación nueva necesita su escala.

## 3.4 Mapa fenómeno → código

Las líneas corresponden a `dae_residual_full_cell.py` en esta copia; si cambian, busque por nombre.

| Fenómeno | Archivo | Función (línea) | ¿Se reemplaza con subclase? |
| --- | --- | --- | --- |
| Balance del tanque | `cells/DAE/dae_residual_full_cell.py` | `tank_residual` (295) | No (función de módulo) |
| Nernst–Planck de V(IV), V(V), H⁺, HSO₄⁻ en el fieltro | ídem | `cathode_residual` → `species_transport` (127, 159) | No |
| Disociación HSO₄⁻ ⇌ H⁺ + SO₄²⁻ | ídem | `b_dis` en `cathode_residual` (218) | No |
| Carga en el sólido (Φ_el) y en el electrolito (Φ_io) del cátodo | ídem | `cathode_residual` (258–292) | No |
| Cinética catódica VO²⁺/VO₂⁺ | `electrochem/butler_volmer.py` | `butler_volmer`, `cathode_surface_concentrations` | **Sí**: `Cell.butler_volmer` |
| Sobrepotencial y OCP | `electrochem/compute_overpotential.py`, `compute_OCP.py` | `compute_overpotential`, `compute_OCP` | **Sí** |
| κ del catolito | `cells/methods/update_kappa_io_m.py` | `update_kappa_io_cathode` | **Sí** |
| κ(λ) de la membrana | ídem | `update_kappa_io_m` | **Sí** |
| Agua en la membrana: difusión + arrastre | `cells/DAE/dae_residual_full_cell.py` | `membrane_residual` (328–428) | No |
| D(λ), ξ(λ), isotermas λ(a) | `cells/methods/membrane_water_diffusivity.py`, `EOD_coefficient.py`, `return_lambda_eq.py` | homónimas | **Sí** |
| Ohm en la membrana | `cells/DAE/dae_residual_full_cell.py` | `membrane_residual` (429–465) | No |
| Ánodo: vapor, H₂, saturación, η (concentrado) | ídem | `anode_residual` (468) | No |
| Cinética HOR (lineal) | `electrochem/butler_volmer.py` | `BV_lumped` | **Sí** |
| Pérdida protónica en la ACL (analítica) | `cells/DAE/dae_residual_full_cell.py` | `update_BC_membrane_anode` (734) | No |
| Interfaz cátodo/membrana: i_io, c_H* | ídem | `update_BC_cathode_membrane` (689), `coupling_residual` (596) | No |
| Velocidad en el fieltro, reparto fieltro/canal | `cells/methods/initialize_properties.py` | `solve_velocity` | **Sí** (se llama en `__init__`) |
| Malla y operador de difusión | `meshes/grid_creation/submesh.py` | `submesh` | — |

## 3.5 El patrón de extensión del curso: subclase de `Cell`

El residuo llama a muchas leyes a través de `self` (por ejemplo `self.update_kappa_io_m(...)`).
Si usted define una subclase que redefine ese método, el modelo usa su versión **sin tocar
`rfbmodel1D/`**:

```python
from rfbmodel1D.cells.cell import Cell

class MiCelda(Cell):
    def update_kappa_io_m(self, lambda_w, Phi_io_membrane):
        super().update_kappa_io_m(lambda_w, Phi_io_membrane)   # ley original
        self.membrane.kappa_io = self.membrane.kappa_io * 0.9    # su modificación

cell, y0, yp0 = construir_celda(0.75, "galvanostatico", -500.0, clase=MiCelda)
```

Ventajas: su cambio queda en un archivo propio, el caso base sigue disponible para comparar y
no hay conflictos al actualizar el repositorio. Para ecuaciones que viven en funciones de módulo
(columna "No" de la tabla) hay que editar `dae_residual_full_cell.py`; ver
[05_modificar_el_modelo.md](05_modificar_el_modelo.md).

## 3.6 Condiciones de borde: `cell.bc`

`cell.bc[variable][subdominio][lado]` es un `boundary_condition(type, value)`, con `variable`
un miembro de `Vars` (`Vars.c_v4`, `Vars.phi_io`, `Vars.a`, ...), `subdominio` en
`"cathode"`, `"membrane"`, `"anode"` y `lado` en `"left"`, `"right"`.

| Condición | Quién la fija | Cuándo |
| --- | --- | --- |
| Φ_el en el colector catódico (corriente o voltaje) | usted, con `modify_BC` (o `construir_celda`) | una vez |
| c_i en x = 0 del cátodo (= tanque) | `update_BC_tank_cathode` | cada evaluación |
| c_H y Φ_io en la cara cátodo/membrana | `update_BC_cathode_membrane` | cada evaluación |
| Actividad de agua a_w del lado cátodo de la membrana | `update_BC_cathode_membrane` (vía SOC medio) | cada evaluación |
| Φ_io, a_w y s del lado ánodo de la membrana | `update_BC_membrane_anode` | cada evaluación |
| Humedad y c_H2 en el canal de H₂ | `params.py` (`bcs[Vars.a]["anode"]["right"]`, `bcs[Vars.c_h2]["anode"]["right"]`) | fijas |

## 3.7 Dónde trabaja cada proyecto

```mermaid
flowchart LR
    subgraph EXT["Fuera del residuo: scripts que llaman a Cell"]
        P02["P02 · Ny segmentos en serie<br/>tutoriales/P02_inicio.py"]
        P28["P28 · ciclado y contabilidad de agua<br/>flujos_agua_membrana"]
        P29["P29 · balance de energía con Q_gen = i·A·(V − OCV)"]
    end
    subgraph SUB["Dentro del residuo: subclase de Cell"]
        P01["P01 · update_kappa_io_m<br/>κ con vanadio (Donnan)"]
        P31["P31 · butler_volmer<br/>β0 = D/δ"]
        P03["P03 · BV_lumped (opcional)<br/>+ modelo 1D del ánodo aparte"]
    end
    EXT --> M["rfbmodel1D (sin cambios)"]
    SUB --> M
```

Ningún proyecto del alcance mínimo necesita editar `rfbmodel1D/`. Las ampliaciones que agregan
una incógnita nueva (temperatura como estado en P29, c_s como incógnita en P31) sí lo requieren;
siga la lista de la sección 5.4 de [05_modificar_el_modelo.md](05_modificar_el_modelo.md).
