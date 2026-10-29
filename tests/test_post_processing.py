import numpy as np
from types import SimpleNamespace

from rfbmodel1D.parameters.params import Vars
from rfbmodel1D.post_processing import calculate_current_density


class DummyBC:
    def __init__(self, value=0.0, bc_type="neumann"):
        self.value = value
        self.type = bc_type


class DummyCell:
    def __init__(self):
        self.cathode = SimpleNamespace(thickness=1.0, eff_conductivity=1.0, kappa_io=np.array([1.0, 1.0]))
        self.membrane = SimpleNamespace(thickness=1.0, kappa_io=np.array([1.0, 1.0]))
        self.mesh = SimpleNamespace(
            cathode_mesh=SimpleNamespace(n_x=2, dx=np.array([0.5, 0.5]), x_grid=np.array([0.0, 0.5, 1.0])),
            membrane_mesh=SimpleNamespace(n_x=2, dx=np.array([0.5, 0.5]), x_grid=np.array([0.0, 0.5, 1.0])),
        )
        self.bc_mode = "galvanostatic"
        self.bc = {
            Vars.phi_el: {"cathode": {"left": DummyBC(0.0), "right": DummyBC(0.0)}, "membrane": {"left": DummyBC(0.0), "right": DummyBC(0.0)}},
            Vars.phi_io: {"cathode": {"left": DummyBC(0.0), "right": DummyBC(0.0)}, "membrane": {"left": DummyBC(0.0), "right": DummyBC(0.0)}},
        }

    def split_state(self, y):
        cathode = np.array([
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ])
        membrane = np.array([
            [0.0, 0.0],
            [0.0, 0.0],
        ])
        return {"cathode": cathode, "membrane": membrane}


def test_calculate_current_density_returns_spatial_arrays():
    cell = DummyCell()
    result = calculate_current_density(cell, np.zeros(0))

    assert set(result.keys()) >= {"x_ca", "x_mem", "i_el_ca", "i_io_ca", "i_io_mem"}
    assert result["i_el_ca"].shape == (3,)
    assert result["i_io_ca"].shape == (3,)
    assert result["i_io_mem"].shape == (3,)
