# 4. Ecuaciones de gobierno, discretización e implementación

Toda la física está en una función, `dae_residual_full_cell` en
`rfbmodel1D/cells/DAE/dae_residual_full_cell.py`. Recibe (t, y, y′) y devuelve el residuo F; el
integrador ajusta y e y′ hasta que F = 0. Cambiar el modelo siempre significa cambiar lo que
devuelve F: usted nunca escribe un lazo de pasos de tiempo.

## 4.1 Incógnitas diferenciales y algebraicas

Una ecuación es **diferencial** si contiene y′ y **algebraica** si no. El integrador trata ambas,
pero la inicialización (`cells/methods/initialize_dae_state.py`) necesita saber cuál es cuál.

| Subdominio | Diferenciales (tienen ∂/∂t) | Algebraicas |
| --- | --- | --- |
| Tanque | c_V4, c_V5, c_H, c_HSO4, V_tanque | — |
| Cátodo | c_V4, c_V5, c_H, c_HSO4 | Φ_el, Φ_io (conservación de carga) |
| Membrana | λ | Φ_io |
| Ánodo (ACL + AGDL) | c_vw, c_H2, s (en ambas capas) | η de la ACL |
| Acoplamiento | — | c_H en la cara cátodo/membrana |

## 4.2 Discretización por volúmenes finitos

El cátodo y la membrana se dividen en volúmenes de control en el intervalo [0, 1], que luego se
escala por el espesor. `submesh` guarda las caras `x_nodes` (N + 1 valores), los centros `x_grid`
(N valores) y los anchos `dx`. `Cell.__init__` los pasa a metros: `_dx_ca` (anchos Δx_k) y
`_dx_face_ca` (distancia entre centros vecinos, Δx_{k+1/2}). La malla del cátodo es no uniforme,
con agrupamiento tanh (α = 1) que refina ambos extremos; la de la membrana es uniforme.

Cada ecuación de la forma ∂(acumulación)/∂t + ∂N/∂x = S se integra sobre el volumen k:

```math
\varepsilon\,\frac{dc_k}{dt}\,\Delta x_k \;+\; \left(N_{k+1/2} - N_{k-1/2}\right) \;-\; S_k\,\Delta x_k \;=\; 0
```

Por eso el residuo tiene unidades de flujo (mol m⁻² s⁻¹ o A m⁻²) y toda fuente volumétrica
aparece multiplicada por `dx`. Las aproximaciones en las caras son:

| Término | Aproximación en la cara | Dónde |
| --- | --- | --- |
| Difusión | N = −D (c_{k+1} − c_k)/Δx_{k+1/2} | `N_e_diff`, `N_w_diff` en `species_transport` |
| Convección a través del espesor | v · c_cara; central si Pe = vΔx/D < 2, si no aguas arriba | `schemes`, `c_e_all`, `c_w_all` (inactiva por defecto: `u_tp = 0`) |
| Migración | −zμF · ½(c_k + c_{k+1}) · (Φ_{k+1} − Φ_k)/Δx_{k+1/2} | `Cell.return_migration_elements`, `apply_migration_operator` |
| Corriente iónica, cátodo | −κ_cara (Φ_{k+1} − Φ_k)/Δx_{k+1/2}, κ_cara media aritmética ponderada por distancia | `kappa_faces` en `cathode_residual` |
| Corriente iónica, membrana | ídem con κ_cara media armónica | `kappa_faces` en `membrane_residual` |
| Corriente electrónica, cátodo | −σ_eff (Φ_{k+1} − Φ_k)/Δx_{k+1/2}, mediante la matriz `Div` | `self._Div_ca` (de `submesh.D`) |
| Agua en la membrana | difusión con D(λ̄_cara) + arrastre, con ξ evaluado aguas arriba según el signo de i_io | `N_e`, `N_w` en `membrane_residual` |

Las condiciones de borde entran como el flujo por la primera o la última cara:

- **Dirichlet** (valor c_b dado): la cara está a medio volumen del centro, N = −D (c_0 − c_b)/(Δx_0/2).
  Ver `Cell.apply_diffusion_left_BC_vec`.
- **Neumann** (flujo dado): el flujo se suma al vector fuente `b`. Ver `Cell.get_boundary_flux_left`.
- **Galvanostática** (solo potenciales): el flujo en la cara es la corriente impuesta.

## 4.3 Cátodo (fieltro de carbono): `cathode_residual`, línea 127

Para i ∈ {V(IV), V(V), H⁺, HSO₄⁻}, balance de Nernst–Planck con un término de intercambio que
representa el flujo en el plano del fieltro:

```math
\varepsilon\frac{\partial c_i}{\partial t} + \frac{\partial N_i}{\partial x} = S_i + \frac{u_{ip}}{H}\left(c_i^{T} - c_i\right),
\qquad
N_i = -D_i^{\mathrm{eff}}\frac{\partial c_i}{\partial x} - z_i\,\mu_i F c_i\frac{\partial \Phi_{io}}{\partial x} + v\,c_i
```

con D_eff = D ε^1,5 (Bruggeman) y μ = D_eff/(RT), calculadas una vez en
`cells/methods/initialize_properties.py`. La velocidad u_ip sale de igualar la caída de presión
en el fieltro y en el canal (`solve_velocity`), que también fija qué fracción del catolito no
entra al fieltro. La función interna `species_transport` devuelve exactamente la forma integrada:

```python
return (porosity * du_c * dx               # acumulación
        + N_flux                           # difusión + convección, este − oeste
        - z * mu * F * mig                 # migración, este − oeste
        - b                                # flujos de Neumann en los bordes
        - source_term * dx                 # fuentes de reacción S_i
        - source_term_convection_ip * dx)  # intercambio en el plano con el tanque
```

Las fuentes S_i se pasan en el argumento `source_term=` (líneas 220–256):

| Especie | z | S_i [mol m⁻³ s⁻¹] |
| --- | --- | --- |
| V(IV), VO²⁺ | +2 | −j_n/F |
| V(V), VO₂⁺ | +1 | +j_n/F |
| H⁺ | +1 | +2 j_n/F + b_dis |
| HSO₄⁻ | −1 | −b_dis, con b_dis = −k_d [(c_H − c_HSO4)/(c_H + c_HSO4) − β] |

Conservación de carga en el sólido (Φ_el, línea 269) y en el electrolito (Φ_io, líneas 290–292):

```math
\frac{\partial}{\partial x}\left(\sigma_{\mathrm{eff}}\frac{\partial \Phi_{el}}{\partial x}\right) = j_n,
\qquad
\frac{\partial}{\partial x}\left(\kappa\,\frac{\partial \Phi_{io}}{\partial x}\right) = -j_n
```

con σ_eff = σ (1 − ε)^1,5 y κ calculada con la composición local en `update_kappa_io_cathode`
(F² Σ z² μ c, con SO₄²⁻ por electroneutralidad). La corriente volumétrica j_n es Butler–Volmer
con corrección por transferencia de masa (`electrochem/butler_volmer.py`):

```math
j_n = aFk\,c_{V5}^{1-\alpha}c_{V4}^{\alpha}c_{H}^{2-2\alpha}\left[\frac{c^{s}_{V4}}{c_{V4}}e^{(1-\alpha)f\eta} - \frac{c^{s}_{V5}}{c_{V5}}e^{-\alpha f\eta}\right],
\qquad \eta = \Phi_{el} - \Phi_{io} - E_{ca}
```

```math
E_{ca} = 1 + \frac{\Delta S_r}{F}\left(T - 298{,}15\right) + \frac{RT}{F}\ln\frac{F_y\,c_{V5}\,c_{H}^{2}}{c_{V4}\,a_w}
```

con f = F/(RT). Las concentraciones de superficie c^s siguen a You et al. (2009), con un
coeficiente de transferencia de masa β₀ ∝ u_ip^0,4 (`cathode_surface_concentrations`).
j_n > 0 es oxidación (carga) y j_n < 0 reducción (descarga).

## 4.4 Tanque: `tank_residual`, línea 295

Volumen bien mezclado que intercambia especies con el fieltro por el flujo en el plano y por la
cara x = 0:

```math
V_T\,\frac{dc_i^{T}}{dt} = -A_{tp}\,N_i\big|_{x=0} - A_{ip}\,u_{ip}\left(c_i^{T} - \bar{c}_i\right)
```

c̄_i es el promedio en volumen del fieltro, A_tp = alto × ancho y A_ip = ancho × espesor. El
volumen V_T es constante (su residuo es dV_T/dt = 0). `update_BC_tank_cathode` (línea 634)
calcula N_i en x = 0 y fija la composición del tanque como condición de Dirichlet del cátodo.

## 4.5 Membrana (Nafion): `membrane_residual`, línea 328

El contenido de agua λ (mol H₂O por mol SO₃⁻) se transporta por difusión y arrastre
electro-osmótico; Φ_io obedece Ohm sin fuentes:

```math
\frac{\partial \lambda}{\partial t} + \frac{\partial N_\lambda}{\partial x} = 0,
\qquad
N_\lambda = -D_\lambda(\lambda)\frac{\partial \lambda}{\partial x} + \frac{\xi(\lambda)\,i_{io}}{F\,c_f},
\qquad
\frac{\partial}{\partial x}\left(\kappa_m(\lambda)\frac{\partial \Phi_{io}}{\partial x}\right) = 0
```

| Ley | Expresión | Archivo |
| --- | --- | --- |
| Arrastre ξ(λ) | 2,5 λ/22 · exp[4000/R (1/303,15 − 1/T)] | `cells/methods/EOD_coefficient.py` |
| Conductividad κ_m(λ) | 50 (f_v − 0,06)^1,5 · exp[15000/R (1/303,15 − 1/T)] S/m, con f_v = λV_w/(λV_w + V_m) | `cells/methods/update_kappa_io_m.py` |
| Difusividad D_λ(λ) | según `membraneparams.diffusivity_model` (por defecto `"hwang"`) | `cells/methods/membrane_water_diffusivity.py` |
| Isotermas λ_eq(a) | vapor/líquido (ánodo) y catolito (cátodo) | `cells/methods/return_lambda_eq.py` |

Las condiciones de borde del agua son leyes de transferencia de masa, suavizadas con mezclas
tanh para que el jacobiano sea continuo. Lado cátodo: N_λ = k_ca f_v (λ_eq,ca − λ_0), diez veces
más rápida en desorción. Lado ánodo: N_λ = (k_vw + k_lw + k_extra)(λ_N − λ_eq,an), con
constantes distintas de sorción y desorción para vapor y líquido. Multiplique N_λ por c_f para
obtener mol H₂O m⁻² s⁻¹.

## 4.6 Capa catalítica y difusora del ánodo: `anode_residual`, línea 468

El ánodo es delgado (ACL 19 µm, AGDL 171 µm): cada capa es un volumen concentrado. Los flujos
ACL–AGDL–canal usan medias distancias.

| Incógnita | Balance (acumulación = …) | Código |
| --- | --- | --- |
| c_vw, ACL | desorción desde la membrana − vapor hacia la AGDL + condensación/evaporación | `res_view[:, 0]` |
| c_H2, ACL | − H₂ hacia la AGDL − consumo por HOR j/(2F) | `res_view[:, 1]` |
| s, ACL | líquido desde la membrana − flujo de Darcy hacia la AGDL − condensación/evaporación | `res_view[:, 2]` |
| η, ACL (algebraica) | j₀ · 2Fη/(RT) + i_io = 0 (HOR linealizada) | `res_view[:, 3]` |
| c_vw, AGDL | desde la ACL − hacia el canal + condensación/evaporación | `res_view[:, 4]` |
| c_H2, AGDL | desde la ACL − hacia el canal | `res_view[:, 5]` |
| s, AGDL | desde la ACL − hacia el canal (p_liq = 103 kPa) − condensación/evaporación | `res_view[:, 6]` |

El potencial electrónico del ánodo es la referencia (Φ_el,an = 0). `update_BC_membrane_anode`
(línea 734) fija Φ_io del lado ánodo de la membrana en −η − E_an − η_protón, donde η_protón es
la solución analítica de un electrodo poroso con cinética lineal (ver la guía de P03), y la
actividad de agua en a = RH + 2 s_r.

## 4.7 Interfaz cátodo/membrana: `coupling_residual`, línea 596

`update_BC_cathode_membrane` (línea 689) trata los dos medios volúmenes como resistencias en serie:

```math
i_{io} = -\frac{g_{ca}\,g_m}{g_{ca} + g_m}\left(\Phi_{io,m,0} - \Phi_{io,ca,N}\right),
\qquad g = \frac{\kappa}{\Delta x/2}
```

La incógnita algebraica extra c_H* (H⁺ del lado cátodo de la cara) se fija exigiendo que el
flujo de protones por esa cara sea i_io/F. El vanadio y el HSO₄⁻ ven una pared de flujo cero
(Neumann) en la membrana: **esta versión no modela el cruce de vanadio**.

## 4.8 Comprobaciones que todo modelo modificado debe pasar

| Prueba | Cómo | Valor de la base |
| --- | --- | --- |
| Arranque consistente | `max(abs(cell.dae_residual_full_cell(0, y0, yp0)))` | 8,4 × 10⁻⁴ |
| Conservación de carga | ∫ j dx en el cátodo vs. i_io vs. corriente aplicada (`perfiles`) | −500,00 / −500,00 / −500 A/m² |
| Faraday en el inventario de V(V) | Δn_V5 (tanque + poros) vs. I t/F | diferencia < 10⁻⁹ mol en 600 s |
| Malla | NC = 100 vs. 250 | ΔV < 0,1 mV |
| Interruptor en cero | su modificación con parámetro nulo debe reproducir la base | idéntico |
