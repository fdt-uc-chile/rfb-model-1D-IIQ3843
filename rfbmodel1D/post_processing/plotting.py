import matplotlib.pyplot as plt


def plot_tank_time_series(t, tank, ax=None):
    """Simple time-series plot for tank species and SOC."""
    if ax is None:
        _, ax = plt.subplots(2, 2, figsize=(10, 7))

    ax[0, 0].plot(t, tank["c_V4"], lw=2, label="V(IV)")
    ax[0, 0].plot(t, tank["c_V5"], lw=2, label="V(V)")
    ax[0, 0].set_ylabel("Concentration [mol/m^3]")
    ax[0, 0].set_title("Tank vanadium species")
    ax[0, 0].legend()
    ax[0, 0].grid(True, alpha=0.3)

    ax[0, 1].plot(t, tank["c_H"], lw=2, label="H+")
    ax[0, 1].plot(t, tank["c_HSO4"], lw=2, label="HSO4-")
    ax[0, 1].set_ylabel("Concentration [mol/m^3]")
    ax[0, 1].set_title("Tank acid species")
    ax[0, 1].legend()
    ax[0, 1].grid(True, alpha=0.3)

    if "soc" in tank:
        ax[1, 0].plot(t, tank["soc"] * 100, lw=2)
        ax[1, 0].set_ylabel("SOC [%]")
        ax[1, 0].set_title("State of charge")
        ax[1, 0].grid(True, alpha=0.3)
    else:
        ax[1, 0].axis("off")

    ax[1, 1].plot(t, tank["V_tank"], lw=2, color="k")
    ax[1, 1].set_ylabel("Volume [m^3]")
    ax[1, 1].set_title("Tank volume")
    ax[1, 1].grid(True, alpha=0.3)

    for axis in ax.flat:
        axis.set_xlabel("Time [s]")

    return ax


def plot_cathode_profiles(x, cathode, ax=None):
    """Simple cathode profile plot for major species and potentials."""
    if ax is None:
        _, ax = plt.subplots(2, 2, figsize=(10, 7))

    ax[0, 0].plot(x, cathode["c_V4"], "o-", lw=2, label="V(IV)")
    ax[0, 0].plot(x, cathode["c_V5"], "s-", lw=2, label="V(V)")
    ax[0, 0].set_title("Cathode vanadium concentrations")
    ax[0, 0].set_ylabel("Concentration [mol/m^3]")
    ax[0, 0].legend()
    ax[0, 0].grid(True, alpha=0.3)

    ax[0, 1].plot(x, cathode["c_H"], "o-", lw=2, label="H+")
    ax[0, 1].plot(x, cathode["c_HSO4"], "s-", lw=2, label="HSO4-")
    ax[0, 1].set_title("Cathode acid concentrations")
    ax[0, 1].set_ylabel("Concentration [mol/m^3]")
    ax[0, 1].legend()
    ax[0, 1].grid(True, alpha=0.3)

    ax[1, 0].plot(x, cathode["phi_el"], "o-", lw=2, label="phi_el")
    ax[1, 0].plot(x, cathode["phi_io"], "s-", lw=2, label="phi_io")
    ax[1, 0].set_title("Cathode potentials")
    ax[1, 0].set_ylabel("Potential [V]")
    ax[1, 0].legend()
    ax[1, 0].grid(True, alpha=0.3)

    ax[1, 1].axis("off")
    return ax


def plot_current_density_profile(x_ca, x_mem, i_io_ca, i_io_mem, ax=None):
    """Simple current-density profile plot across cathode and membrane."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))

    ax.plot(x_ca, i_io_ca, "o-", lw=2, label="Cathode ionic current")
    ax.plot(x_mem, i_io_mem, "s-", lw=2, label="Membrane ionic current")
    ax.set_xlabel("Position")
    ax.set_ylabel("Current density [A/m^2]")
    ax.set_title("Ionic current density profile")
    ax.grid(True, alpha=0.3)
    ax.legend()
    return ax
