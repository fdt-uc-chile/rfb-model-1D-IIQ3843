"""Generate the membrane lambda-profile result files used by the plotting notebook.

This script mirrors the publication notebook's parameter-refresh pattern:
1) load the latest HFR fit JSON for each test,
2) load the latest matching global E_cell fit JSON,
3) update the model parameters in-memory,
4) re-run the DAE at the selected operating points,
5) save the lambda profiles to the expected npz file names.

The saved arrays match what the notebook expects:
    x_mem, lam_ch, lam_dis, j_charge_lam, j_discharge_lam

python scripts/utils/generate_lambda_profiles.py --save
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from solve_dae.integrate import solve_dae

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(ROOT / "scripts" / "utils"))
from simulate_HFR import build_cell, build_atol, soc_cutoff_event  # noqa: E402
from rfbmodel1D.parameters.params import (  # noqa: E402
    anodeparams,
    cathodeparams,
    cellparams,
    membraneparams,
)

TF = 1800
RTOL = 1e-5

DEFAULT_TEST_NOS = [5, 16, 25]


def _representative_indices(n_pts: int, n_curves: int = 4) -> list[int]:
    """Return evenly spaced representative indices spanning the full range.

    This keeps a minimum, a maximum, and intermediate values while avoiding stale
    hard-coded indices that may no longer be valid for the current CSV layout.
    """
    if n_pts <= 0:
        return []
    if n_pts <= n_curves:
        return list(range(n_pts))
    idx = np.linspace(0, n_pts - 1, n_curves, dtype=int)
    return sorted({int(i) for i in idx})


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tests", type=int, nargs="+", default=None,
                   help="Restrict generation to these test numbers (default: all of [5, 16, 25]).")
    p.add_argument("--save", action="store_true",
                   help="Optionally persist npz outputs after the forward re-simulation. Default is to run the solver only.")
    return p.parse_args()


def _valid_indices(indices: list[int], n_pts: int) -> list[int]:
    valid = [idx for idx in indices if 0 <= idx < n_pts]
    if len(valid) != len(indices):
        dropped = sorted(set(indices) - set(valid))
        print(f"  note: dropping invalid indices {dropped} for {n_pts} available points")
    return valid


def _latest_global_fit(test_nos: list[int]) -> tuple[Path, dict]:
    global_dir = ROOT / "results" / "test fitting" / "global"
    candidates = sorted(global_dir.glob("global_k_ca_fit_*.json"))
    for path in reversed(candidates):
        with open(path, "r") as f:
            fit = json.load(f)
        if fit.get("algorithm") != "grid_refine_v2":
            continue
        fit_tests = {int(t) for t in fit.get("tests", [])}
        if set(test_nos).issubset(fit_tests):
            return path, fit
    raise FileNotFoundError(f"No compatible global fit found for tests {test_nos} in {global_dir}")


def _latest_hfr_fit(test_no: int) -> Path:
    folder = ROOT / "results" / "test fitting" / f"test {test_no}"
    candidates = sorted(folder.glob(f"HFR_fit_test{test_no}_*.json"))
    if not candidates:
        raise FileNotFoundError(f"No HFR fit JSON found for test {test_no} in {folder}")
    return candidates[-1]


def _apply_fit_params(test_no: int, k_ca_global: float, j0_an_by_test: dict[int, float], hfr_fit_path: Path) -> None:
    with open(hfr_fit_path, "r") as f:
        hfr_fit = json.load(f)
    p, s = hfr_fit["params"], hfr_fit["settings"]

    membraneparams.k_sorp_vw = p["k_vw_sorp"]
    membraneparams.k_desorp_vw = p["S_desorp"] * p["k_vw_sorp"]
    membraneparams.k_sorp_lw = p["k_lw_sorp"]
    membraneparams.k_desorp_lw = p["S_desorp"] * p["k_lw_sorp"]
    membraneparams.k_ca_lw = p["k_ca"]
    membraneparams.sf_D = s["sf_D"]
    cellparams.vanadium_uptake = bool(s.get("vanadium_uptake", False))

    cathodeparams.reaction_rate_constant = float(k_ca_global)
    anodeparams.j_0 = float(j0_an_by_test[int(test_no)])


def _solve_profile(soc: float, T: float, Q_v: float, j_appl: float):
    cell, y0, yp0 = build_cell(soc, j_appl, T, Q_v)
    atol = build_atol(cell)
    jac_sparsity = cell.compute_jac_sparsity(y0, yp0)
    sol = solve_dae(
        cell.dae_residual_full_cell,
        (0.0, TF), y0, yp0,
        method="BDF",
        atol=atol,
        rtol=RTOL,
        events=soc_cutoff_event,
        jac_sparsity=jac_sparsity,
    )
    if sol.t[-1] < 0.99 * TF:
        print(f"  warning: early stop at t={sol.t[-1]:.1f}s for j={j_appl:+.0f} A/m²")
    state = cell.split_state(sol.y[:, -1])
    lam = np.asarray(state["membrane"][:, 0], dtype=float)
    s_acl = float(state["anode"][2])
    s_agdl = float(state["anode"][6])
    x_mem = np.asarray(cell.mesh.membrane_mesh.x_grid, dtype=float)
    return x_mem, lam, s_acl, s_agdl


def main():
    args = parse_args()
    TEST_NOS = args.tests if args.tests else DEFAULT_TEST_NOS

    ecell_csv = ROOT / "data" / "experimental" / "E_cell.csv"
    ecell_data = pd.read_csv(ecell_csv)
    global_fit_path, global_fit = _latest_global_fit(TEST_NOS)
    k_ca_global = float(global_fit["params"]["k_ca"])
    j0_an_by_test = {int(k): float(v) for k, v in global_fit["params"]["j0_an"].items()}
    print(f"Using global fit: {global_fit_path.name}")
    print(f"  shared k_ca = {k_ca_global:.4e} m/s")

    for test_no in TEST_NOS:
        test_rows = ecell_data[ecell_data["test_no"] == test_no]
        soc = float(test_rows["SOC"].iloc[0])
        T = float(test_rows["T(°C)"].iloc[0])
        Q_v = float(test_rows["Q_v (ml/min)"].iloc[0])

        charge_j = test_rows[test_rows["mode"] == "charge"]["j_appl (A/m2)"].to_numpy(dtype=float)
        discharge_j = test_rows[test_rows["mode"] == "discharge"]["j_appl (A/m2)"].to_numpy(dtype=float)

        hfr_fit_path = _latest_hfr_fit(test_no)
        _apply_fit_params(test_no, k_ca_global, j0_an_by_test, hfr_fit_path)

        x_mem = None
        lam_ch = []
        lam_dis = []
        s_acl_ch = []
        s_agdl_ch = []
        s_acl_dis = []
        s_agdl_dis = []
        j_charge_lam = []
        j_discharge_lam = []

        print(f"\nTest {test_no}: SOC={soc}  T={T}°C  Q_v={Q_v} ml/min")
        charge_indices = _valid_indices(_representative_indices(len(charge_j), 4), len(charge_j))
        for idx in charge_indices:
            j = float(charge_j[idx])
            x_mem, lam, s_acl, s_agdl = _solve_profile(soc, T, Q_v, float(j))
            lam_ch.append(lam)
            s_acl_ch.append(s_acl)
            s_agdl_ch.append(s_agdl)
            j_charge_lam.append(j)
            print(f"  charge   j={j:+.0f} A/m²  lambda0={lam[0]:.3f}  lambdaL={lam[-1]:.3f}  s_acl={s_acl:.3f}  s_agdl={s_agdl:.3f}")

        discharge_indices = _valid_indices(_representative_indices(len(discharge_j), 4), len(discharge_j))
        for idx in discharge_indices:
            j = float(discharge_j[idx])
            x_mem, lam, s_acl, s_agdl = _solve_profile(soc, T, Q_v, -float(j))
            lam_dis.append(lam)
            s_acl_dis.append(s_acl)
            s_agdl_dis.append(s_agdl)
            j_discharge_lam.append(float(j))
            print(f"  discharge j={-j:+.0f} A/m²  lambda0={lam[0]:.3f}  lambdaL={lam[-1]:.3f}  s_acl={s_acl:.3f}  s_agdl={s_agdl:.3f}")

        if args.save:
            out_dir = ROOT / "results" / "test fitting" / f"test {test_no}"
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"test{test_no}_results.npz"
            payload = {
                "x_mem": np.asarray(x_mem, dtype=float),
                "lam_ch": np.asarray(lam_ch, dtype=float),
                "lam_dis": np.asarray(lam_dis, dtype=float),
                "s_acl_ch": np.asarray(s_acl_ch, dtype=float),
                "s_agdl_ch": np.asarray(s_agdl_ch, dtype=float),
                "s_acl_dis": np.asarray(s_acl_dis, dtype=float),
                "s_agdl_dis": np.asarray(s_agdl_dis, dtype=float),
                "j_charge_lam": np.asarray(j_charge_lam, dtype=float),
                "j_discharge_lam": np.asarray(j_discharge_lam, dtype=float),
            }
            np.savez(out_path, **payload)
            print(f"Saved -> {out_path}")
        else:
            print(f"Simulation rerun complete for test {test_no} without writing .npz cache files.")


if __name__ == "__main__":
    main()
