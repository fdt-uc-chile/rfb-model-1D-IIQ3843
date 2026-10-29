# Proyectos individuales con rfbmodel1D (HVRFB)

Versión 2 de las guías (29-10-2026). Reemplaza la sección "Guía para realizarlo" de la versión 1
(28-10-2026) para los seis proyectos que usan este modelo. Lo demás de su guía (contexto,
retroalimentación de la CR1, hitos) sigue vigente.

## Qué cambió respecto de la versión 1

**Se resolvieron los pendientes del equipo docente:**

| Pendiente en la v1 | Estado |
| --- | --- |
| Tutorial de la corrida base de rfb-model | [docs/02_corrida_base.md](../02_corrida_base.md) y `tutoriales/01_corrida_base_perfiles.py` |
| Descripción de la estructura del código | [docs/03_arquitectura.md](../03_arquitectura.md), con diagramas UML |
| Confirmar que el registro Zenodo 22274872 sea público | Sí: acceso abierto, CC BY 4.0, publicado el 03-09-2026; código con licencia MIT |
| Nombres de archivos y funciones (antes, términos de búsqueda) | En la sección "Dónde está su proyecto en el código" de cada guía |
| Tesis de J. Barros en `refs/electrochem` | **Sigue pendiente.** Las guías no dependen de ella |

**Los alcances se redujeron.** Quedan unas dos semanas hasta la EPI (13-11) y el código es
complejo. Cada guía trae ahora:

- una pregunta y tres objetivos más acotados;
- un script de inicio que **ya corre** y deja marcado con `TODO` lo que usted programa;
- un paso a paso con archivo, función y el resultado numérico que debe obtener en cada paso;
- ningún cambio dentro de `rfbmodel1D/` en el alcance mínimo: todo se hace con una subclase de
  `Cell` o con un script externo (ver [§3.5 de la arquitectura](../03_arquitectura.md#35-el-patrón-de-extensión-del-curso-subclase-de-cell)).

## Resumen

| Proyecto | Pregunta reducida | Enganche con el código | Script |
| --- | --- | --- | --- |
| [P01](P01_uptake_vanadio.md) De Col | ¿Cuánto baja κ de la membrana entre SoC 0,2 y 0,9 por el vanadio absorbido (Donnan) y cuánto cambia la caída óhmica a 500 A/m²? | Subclase: `update_kappa_io_m` | `P01_inicio.py` |
| [P02](P02_pseudo2D.md) Rifo | ¿Cuánto se desvía la corriente total de una celda con Ny segmentos en serie respecto del modelo 1D, según caudal? | Script externo que llama a la celda Ny veces | `P02_inicio.py` |
| [P03](P03_anodo.md) Rebolledo | ¿Cómo cambian η_act y η_protón del ánodo con la carga de Pt y el espesor de ACL, y cuánto agrega la AGDL? | Barrido de parámetros + modelo 1D del ánodo aparte | `P03_inicio.py` |
| [P28](P28_agua.md) Carcelén | ¿Cuánta agua gana el catolito por ciclo y cuánta habría que reponer al día? | Post-proceso: `flujos_agua_membrana` | `P28_inicio.py` |
| [P29](P29_termico.md) Navarro | ¿Supera el electrolito 40 °C en un día de Calama con enfriamiento pasivo? | Post-proceso: Q_gen = i·A·(V − OCV) + balance de energía aparte | `P29_inicio.py` |
| [P31](P31_corriente_terciaria.md) Mariné | ¿Cuánto cambian j(x) y η al pasar del límite secundario a la capa límite δ? | Subclase: `butler_volmer` | `P31_inicio.py` |

## Flujo de trabajo común

1. Instale y corra el Tutorial 1 ([docs/01](../01_instalacion.md), [docs/02](../02_corrida_base.md)).
2. Cree su rama: `git checkout -b P0X-apellido`.
3. Corra su script de inicio sin cambios y compare con los valores de su guía.
4. Copie el script a `src/` (o trabaje sobre él) y complete los `TODO` en el orden de la guía.
5. Antes de cada figura, pase las pruebas de la sección "Verificación" de su guía.

## Plan de tareas reprogramado (32 h, mismas horas por fase)

| Fase | Horas | Semana | Contenido común |
| --- | --- | --- | --- |
| F1 | 4 | S2 (26-10) | Instalar, Tutorial 1, script de inicio, ubicar su función en el código |
| F2 | 6 | S2 (26-10) | Ficha matemática y tabla de parámetros |
| F3 | 10 | S3 (02-11) | Programar y verificar |
| F4 | 8 | S3–S4 | Corridas y figuras |
| F5 | 4 | S4 (09-11) | Informe EPI (13-11, 17:30) y CR2 (17-11) |

Copie las tareas detalladas de su guía en la hoja Gantt de su planilla `plan_<Apellido>_P<ID>.xlsx`.
