"""Reusable post-processing tools for extracted model states and solution outputs."""

from .profiles import (
    extract_cathode_profiles,
    extract_membrane_profiles,
    extract_tank_time_series,
    extract_full_cell_time_series,
)
from .currents import calculate_current_density
from .plotting import (
    plot_tank_time_series,
    plot_cathode_profiles,
    plot_current_density_profile,
)

__all__ = [
    "extract_cathode_profiles",
    "extract_membrane_profiles",
    "extract_tank_time_series",
    "extract_full_cell_time_series",
    "calculate_current_density",
    "plot_tank_time_series",
    "plot_cathode_profiles",
    "plot_current_density_profile",
]
