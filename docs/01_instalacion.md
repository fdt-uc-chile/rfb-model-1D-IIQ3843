# 1. Instalación

Necesita Python 3.10 o más reciente. Los tutoriales se probaron con Python 3.13 (Miniconda,
Windows 11) el 29-10-2026. No necesita WSL ni compilar nada.

## Paso 1 — Obtener el código

Con git (recomendado, así podrá trabajar en una rama propia):

```bash
git clone <URL del repositorio del curso>
cd rfb-model-1D-IIQ3843
```

Sin git: descargue el ZIP desde GitHub (botón **Code → Download ZIP**) y descomprímalo.

## Paso 2 — Crear un entorno e instalar

Con conda (Anaconda o Miniconda):

```bash
conda create -n rfb python=3.12 -y
conda activate rfb
pip install -e .
pip install -r requirements.txt
```

Con venv:

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate      macOS/Linux:  source .venv/bin/activate
pip install -e .
pip install -r requirements.txt
```

`pip install -e .` instala el paquete `rfbmodel1D` en modo editable: los cambios que haga en
`rfbmodel1D/` se usan de inmediato, sin reinstalar. `solve_dae` es el integrador BDF para
sistemas algebraico-diferenciales que usa el modelo.

## Paso 3 — Verificar

```bash
python -c "import rfbmodel1D, solve_dae; print(rfbmodel1D.__file__)"
python -m pytest -q tests
python tutoriales/01_corrida_base_perfiles.py
```

1. La ruta impresa debe terminar en `rfb-model-1D-IIQ3843/rfbmodel1D/__init__.py`. Si apunta a
   otra carpeta (por ejemplo, una copia antigua de `rfbmodel1D-master.zip`), desinstálela con
   `pip uninstall rfbmodel1D` y repita el paso 2.
2. `pytest` debe informar `1 passed`.
3. El tutorial debe imprimir `éxito=True` y crear `resultados/tutorial1_perfiles.png`.

## Problemas frecuentes

| Síntoma | Causa | Solución |
| --- | --- | --- |
| `Python was not found; run without arguments to install from the Microsoft Store` | Windows usa el alias de la tienda | Active su entorno conda o desactive el alias en *Configuración → Aplicaciones → Alias de ejecución* |
| `ModuleNotFoundError: No module named 'rfbmodel1D'` | Instaló en otro entorno, o no corrió `pip install -e .` | Active el entorno correcto y repita el paso 2 |
| `ModuleNotFoundError: No module named 'iiq3843_utils'` | Ejecutó un tutorial copiado fuera de `tutoriales/` | Ejecute los scripts desde `tutoriales/` o copie también `iiq3843_utils.py` |
| La ventana de la figura no aparece | Backend sin interfaz (servidor, WSL) | La figura igual se guarda en `resultados/`; o use `MPLBACKEND=Agg` |
| Advertencias `RuntimeWarning: invalid value encountered in scalar power` al inicio | El resolvedor de condiciones iniciales prueba estados no físicos | Son inofensivas si al final aparece `éxito=True` |

## Editor

Cualquier editor sirve. En VS Code: abra la carpeta del repositorio, elija el intérprete del
entorno `rfb` (*Python: Select Interpreter*) y use *Run Python File*. Para ubicar una función,
use *Go to Definition* (F12) sobre su nombre o la búsqueda global (Ctrl+Shift+F).
