"""Functions to visualize grid, parameter and output."""
# Third-party
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# Local
from .fluxes import hydraulics_mvg
from .fluxes import hydraulics_rjitema
from .grid import calc_jacobian
from .grid import setup_grid
from .parameter import grid_params
from .parameter import hydraulic_params
from .parameter import hydraulic_params_mvg
from .parameter import options


def plot_grid():

    z, z_m, _, _ = setup_grid(grid_params)

    plt.hlines(z, xmin=0.0, xmax=1.0, linestyle="--")
    plt.hlines(z_m, xmin=0.0, xmax=1.0)

    for layer in z:
        plt.text(0.1, layer, str(layer) + " m", va="top")

    plt.title(
        f"Grid with $\\Delta z_1 = {grid_params['z1']}$ m, $b = {grid_params['b']}$"
        f"and $nz = {grid_params['nz']}$"
    )

    plt.ylabel(r"$z$ [m]")

    # remove xticks
    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.show()


def plot_jacobian():

    z, z_m, _, _ = setup_grid(grid_params)

    metric_zmid = calc_jacobian(z_m, grid_params)
    metric_zhalf = calc_jacobian(z, grid_params)
    plt.plot(
        z_m,
        metric_zmid,
        linestyle="",
        marker="x",
        label=r"$\frac{\partial \xi}{\partial z}$ at Full Layers",
    )
    plt.plot(
        z,
        metric_zhalf,
        linestyle="",
        marker="^",
        label=r"$\frac{\partial \xi}{\partial z}$ at Half Layers",
    )
    plt.xlabel(r"$z$ [m]")
    plt.ylabel(r"Metric Term $[\mathrm{m}^{-1}]$")
    plt.grid(True)
    plt.title("Metric Term for the Prescribed Grid")
    plt.legend()
    plt.close()


# pylint: disable=too-many-locals
def plot_conductivities():

    theta_max = hydraulic_params["pore_volume"]
    theta = np.linspace(0.0, theta_max, 10000)
    j = slice(1, len(theta))
    jm1 = slice(0, len(theta) - 1)
    theta_interpol = np.zeros(len(theta) + 1)
    theta_interpol[j] = 0.5 * (theta[j] + theta[jm1])
    theta_interpol[len(theta)] = theta[-1]

    theta_max_m = hydraulic_params_mvg["pore_volume"]
    theta_m = np.linspace(0.0, theta_max_m, 10000)
    j = slice(1, len(theta_m))
    jm1 = slice(0, len(theta_m) - 1)
    theta_interpol_m = np.zeros(len(theta_m) + 1)
    theta_interpol_m[j] = 0.5 * (theta_m[j] + theta_m[jm1])
    theta_interpol_m[len(theta_m)] = theta_m[-1]

    k0_rij = np.zeros(len(theta) + 1)
    k_rij = np.zeros(len(theta) + 1)
    d_rij = np.zeros(len(theta) + 1)
    k0_mvg = np.zeros(len(theta_m) + 1)
    k_mvg = np.zeros(len(theta_m) + 1)
    d_mvg = np.zeros(len(theta_m) + 1)

    k0_rij[:] = hydraulic_params["k0"]
    k0_mvg[:] = hydraulic_params_mvg["k0"]

    k_rij, d_rij = hydraulics_rjitema(theta, k0_rij, k_rij, d_rij, hydraulic_params)
    k_mvg, d_mvg = hydraulics_mvg(theta_m, k0_mvg, k_mvg, d_mvg, hydraulic_params_mvg)

    # because the functions implement bdries, we have to reset them for a proper
    # plot, the resulting inaccuracy is not visible on the plot.
    k_rij[0] = k_rij[1]
    d_rij[0] = d_rij[1]
    k_rij[-1] = k_rij[-2]
    d_rij[-1] = d_rij[-2]
    k_mvg[0] = k_mvg[1]
    d_mvg[0] = d_mvg[1]
    k_mvg[-1] = k_mvg[-2]
    d_mvg[-1] = d_mvg[-2]

    plt.plot(
        theta_interpol,
        k_rij,
        label=r"$K_\mathrm{RIJTEMA}$",
        color="red",
        alpha=0.5,
        linestyle="--",
    )
    plt.plot(
        theta_interpol,
        d_rij,
        label=r"$D_\mathrm{RIJTEMA}$",
        color="blue",
        alpha=0.5,
        linestyle="--",
    )
    plt.plot(theta_interpol_m, k_mvg, label=r"$K_\mathrm{MVG}$", color="red")
    plt.plot(theta_interpol_m, d_mvg, label=r"$D_\mathrm{MVG}$", color="blue")
    plt.xlabel(r"$\theta$ $[\mathrm{m^3}/\mathrm{m^3}]$")
    plt.ylabel(r"$K$ $[\mathrm{m}/\mathrm{s}]$, $D$ $[\mathrm{m}^2/\mathrm{s}]$")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.show()


def prepare_netcdf(output):

    dataset = xr.Dataset(
        coords=dict(zip(("t", "z_mid"), (output["t_ax"], output["z_mid"]))),
    )

    return dataset


def plot_var_1d(output, out_ds, varname, dimension):

    var = output[varname][:, 0]

    out_da = xr.DataArray(
        var,
        dims=[dimension.split("_")[0]],
        coords={dimension.split("_")[0]: output[dimension]},
    )

    out_da.plot()
    ax = plt.gca()
    if dimension in ["z", "z_mid"]:
        ax.set_yticks(output["z_mid"])
        ax.set_ylim(output["z"][-1], output["z"][0])
    plt.grid(True)
    plt.title(varname)
    plt.close()
    out_ds[varname] = out_da

    return out_ds


def plot_var_massgrid(output, out_ds, varname):

    var = output[varname]

    out_da = xr.DataArray(
        var,
        coords=(output["t_ax"], output["z_mid"]),
        dims=("t", "z"),
    )

    out_da.T.plot()
    ax = plt.gca()
    ax.set_yticks(output["z_mid"])
    ax.set_ylim(output["z"][-1], output["z"][0])
    plt.title(varname)
    plt.close()
    out_ds[varname] = out_da

    return out_ds


def plot_saturation(output, out_ds):
    if options["hydparam"] == "rjitema":
        output.update({"sat": (output["w_vol"] / hydraulic_params["pore_volume"])})
    elif options["hydparam"] == "mvg":
        output.update({"sat": (output["w_vol"] / hydraulic_params_mvg["pore_volume"])})
    plot_var_massgrid(output, out_ds, "sat")

    return out_ds


def plot_var_fluxgrid(output, out_ds, varname):

    var = output[varname]

    out_da = xr.DataArray(
        var,
        coords=(output["t_ax"], output["z"]),
        dims=("t", "z_h"),
    )

    out_da.T.plot()
    ax = plt.gca()
    ax.set_yticks(output["z"])
    ax.set_ylim(output["z"][-1], output["z"][0])
    plt.title(varname)
    plt.close()
    out_ds[varname] = out_da

    return out_ds


if __name__ == "__main__":

    plot_grid()
    plot_jacobian()
