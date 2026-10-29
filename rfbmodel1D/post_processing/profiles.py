import numpy as np


def extract_cathode_profiles(cell, y):
    """Extract the cathode state block at a single time."""
    Nc = cell.mesh.cathode_mesh.n_x
    n_ca = cell.layout.n_vars_cathode
    n_tank = cell.layout.n_vars_tank
    off_ca = n_tank
    off_mem = off_ca + Nc * n_ca
    y_ca = y[off_ca:off_mem].reshape(Nc, n_ca)
    return {
        "x": np.asarray(cell.mesh.cathode_mesh.x_grid),
        "c_V4": y_ca[:, 0],
        "c_V5": y_ca[:, 1],
        "c_H": y_ca[:, 2],
        "c_HSO4": y_ca[:, 3],
        "phi_el": y_ca[:, 4],
        "phi_io": y_ca[:, 5],
    }


def extract_membrane_profiles(cell, y):
    """Extract the membrane state block at a single time."""
    Nc = cell.mesh.cathode_mesh.n_x
    Nm = cell.mesh.membrane_mesh.n_x
    n_ca = cell.layout.n_vars_cathode
    n_mem = cell.layout.n_vars_membrane
    n_tank = cell.layout.n_vars_tank
    off_ca = n_tank
    off_mem = off_ca + Nc * n_ca
    y_mem = y[off_mem:off_mem + Nm * n_mem].reshape(Nm, n_mem)
    return {
        "x": np.asarray(cell.mesh.membrane_mesh.x_grid),
        "lambda": y_mem[:, 0],
        "phi_io": y_mem[:, 1],
    }


def extract_tank_time_series(sol, cell):
    """Extract tank variables over time from a solution object."""
    tank = sol.y[: cell.layout.n_vars_tank, :]
    return {
        "t": sol.t,
        "c_V4": tank[0, :],
        "c_V5": tank[1, :],
        "c_H": tank[2, :],
        "c_HSO4": tank[3, :],
        "V_tank": tank[4, :],
    }


def extract_full_cell_time_series(sol, cell):
    """Extract common time-series vectors used in the full-cell notebook."""
    Nc = cell.mesh.cathode_mesh.n_x
    Nm = cell.mesh.membrane_mesh.n_x
    n_tank = cell.layout.n_vars_tank
    n_ca = cell.layout.n_vars_cathode
    n_mem = cell.layout.n_vars_membrane
    n_an = cell.layout.n_vars_anode
    off_ca = n_tank
    off_mem = off_ca + Nc * n_ca
    off_an = off_mem + Nm * n_mem
    y_ca_all = sol.y[off_ca:off_mem, :].reshape(Nc, n_ca, -1)
    y_mem_all = sol.y[off_mem:off_an, :].reshape(Nm, n_mem, -1)
    y_an_all = sol.y[off_an:off_an + n_an, :]
    tank = sol.y[:n_tank, :]

    c_H_ca = y_ca_all[:, 2, :]
    soc = np.divide(tank[1, :], tank[0, :] + tank[1, :], out=np.zeros_like(tank[1, :]), where=(tank[0, :] + tank[1, :]) > 0)

    return {
        "t": sol.t,
        "tank": {
            "c_V4": tank[0, :],
            "c_V5": tank[1, :],
            "c_H": tank[2, :],
            "c_HSO4": tank[3, :],
            "V_tank": tank[4, :],
            "soc": soc,
        },
        "cathode": {
            "c_H_first": c_H_ca[0, :],
            "c_H_last": c_H_ca[-1, :],
        },
        "membrane": {
            "lambda_mean": np.mean(y_mem_all[:, 0, :], axis=0),
            "lambda_min": np.min(y_mem_all[:, 0, :], axis=0),
            "lambda_max": np.max(y_mem_all[:, 0, :], axis=0),
        },
        "anode": {
            "c_vw_acl": y_an_all[0, :],
            "c_H2_acl": y_an_all[1, :],
            "s_acl": y_an_all[2, :],
            "eta_acl": y_an_all[3, :],
            "c_vw_agdl": y_an_all[4, :],
            "c_H2_agdl": y_an_all[5, :],
            "s_agdl": y_an_all[6, :],
        },
    }
