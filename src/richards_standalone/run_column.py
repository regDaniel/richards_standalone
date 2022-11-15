"""Run the single column model.

Currently only 1D available

"""
# Third-party
import numpy as np

# First-party
from richards_standalone.initialize import initial_condition

# Local
from .fluxes import calc_k_decharme
from .fluxes import est_flux_fg
from .fluxes import est_qground_fg
from .fluxes import flux_corr
from .fluxes import flux_corr_transformed
from .fluxes import flux_fg_transformed
from .fluxes import hydraulics_mvg
from .fluxes import hydraulics_rjitema
from .grid import calc_dzeta
from .grid import calc_jacobian
from .grid import setup_grid
from .grid import transform_back
from .grid import transform_scalar
from .output import init_output
from .output import save_timelevel
from .parameter import hydraulic_params_mvg
from .parameter import options
from .parameter import plant_params
from .parameter import precision_params
from .parameter import time_params
from .parameter import update_parameter
from .plots import plot_saturation
from .plots import plot_var_1d
from .plots import plot_var_fluxgrid
from .plots import plot_var_massgrid
from .plots import prepare_netcdf
from .solver import update_implicit
from .solver import update_implicit_tr

# set precision globally
FloatWP = precision_params["working_precision"]


def run_nolateral(infil, update_dict={}):

    # check user parameter
    time_params, grid_params, hydraulic_params, runoff_params = update_parameter(
        update_dict
    )

    # preparation
    z, z_mid, dz, dz_h = setup_grid(grid_params)
    metric_zmid = calc_jacobian(z_mid, grid_params)
    metric_zhalf = calc_jacobian(z, grid_params)
    dzeta = calc_dzeta(grid_params)

    # initial profile, arrays
    w_vol, w_vol_new, k, d, d_flux, qground, qsurf, wt_depth = initial_condition(
        nsat_layers=0, sat_init=FloatWP(0.8)
    )
    k_flux = np.zeros(len(z), dtype=FloatWP)
    diff_tend = np.zeros(len(w_vol), dtype=FloatWP)
    k_tend = np.zeros(len(w_vol), dtype=FloatWP)

    # retrieve some parameter
    if options["hydparam"] == "rjitema":
        k_sat = hydraulic_params["k0"]
    elif options["hydparam"] == "mvg":
        k_sat = hydraulic_params_mvg["k0"]
    z_root = plant_params["z_root"]
    f = plant_params["f"]

    # saturated hydraulic conductivity
    k0_half = calc_k_decharme(z, k_sat, f, z_root)
    k0_mid = calc_k_decharme(z_mid, k_sat, f, z_root)

    # output
    output_list = [
        "w_vol",
        "qground",
        "wt_depth",
        "k_flux",
        "k_tend",
        "d_flux",
        "diff_tend",
        "qsurf",
    ]
    output_vars = [w_vol, qground, wt_depth, k_flux, k_tend, d_flux, diff_tend, qsurf]
    outvars_dict = dict(zip(output_list, output_vars))
    nout = int(time_params["t_max"] / time_params["dt_out"]) + 1
    output = init_output(outvars_dict, nout, grid_params["nz"])
    print("Output Variables {}".format(output.keys()))
    output = save_timelevel(output, 0, outvars_dict)  # save initial conds
    tlevel = 1  # 0 is reserved for initial conds

    t_max = time_params["t_max"]
    dt = time_params["dt"]
    nt = int(t_max / dt)
    out_interval = int(time_params["dt_out"] / dt)

    for timestep in range(nt):

        # First Guess Fluxes and Ground Runoff
        if options["hydparam"] == "rjitema":
            k, d = hydraulics_rjitema(w_vol, dz, k0_half, k, d, hydraulic_params)
        elif options["hydparam"] == "mvg":
            k, d = hydraulics_mvg(w_vol, dz, k0_half, k, d, hydraulic_params_mvg)

        if options["discretization"] == "transformed":
            d_flux[:], k_flux[:] = flux_fg_transformed(k, d, d_flux, w_vol, dzeta)
        else:
            d_flux[:], k_flux[:] = est_flux_fg(k, d, d_flux, w_vol, dz_h)

        if options["hydparam"] == "rjitema":
            qground, wt_depth[0] = est_qground_fg(
                qground, w_vol, dz, z, k0_mid, runoff_params, hydraulic_params
            )
        elif options["hydparam"] == "mvg":
            qground, wt_depth[0] = est_qground_fg(
                qground, w_vol, dz, z, k0_mid, runoff_params, hydraulic_params_mvg
            )

        # Flux Corrections
        if (
            options["hydparam"] == "rjitema"
            and not options["discretization"] == "transformed"
        ):
            k_flux[:], d, qground, qsurf[0] = flux_corr(
                k_flux,
                d_flux,
                d,
                infil[timestep],
                qground,
                w_vol,
                dz,
                dz_h,
                time_params,
                hydraulic_params,
                runoff_params,
            )
        elif (
            options["hydparam"] == "mvg"
            and not options["discretization"] == "transformed"
        ):
            k_flux[:], d, qground, qsurf[0] = flux_corr(
                k_flux,
                d_flux,
                d,
                infil[timestep],
                qground,
                w_vol,
                dz,
                dz_h,
                time_params,
                hydraulic_params_mvg,
                runoff_params,
            )
        elif (
            options["hydparam"] == "rjitema"
            and options["discretization"] == "transformed"
        ):
            qground_tr = transform_scalar(qground, metric_zmid)
            w_vol_tr = transform_scalar(w_vol, metric_zmid)
            k_flux[:], d, qground_tr, qsurf[0] = flux_corr_transformed(
                k_flux,
                d_flux,
                d,
                infil[timestep],
                qground_tr,
                w_vol_tr,
                dzeta,
                time_params,
                hydraulic_params,
                runoff_params,
                metric_zmid,
            )
        elif (
            options["hydparam"] == "mvg" and options["discretization"] == "transformed"
        ):
            qground_tr = transform_scalar(qground, metric_zmid)
            w_vol_tr = transform_scalar(w_vol, metric_zmid)
            k_flux[:], d, qground_tr, qsurf[0] = flux_corr_transformed(
                k_flux,
                d_flux,
                d,
                infil[timestep],
                qground_tr,
                w_vol_tr,
                dzeta,
                time_params,
                hydraulic_params_mvg,
                runoff_params,
                metric_zmid,
            )

        # Solve LSOE
        if options["discretization"] == "transformed":
            w_vol_new, k_tend[:], diff_tend[:] = update_implicit_tr(
                dt,
                dzeta,
                metric_zhalf,
                metric_zmid,
                d,
                k_flux,
                w_vol,
                qground_tr,
                w_vol_new,
            )
        else:
            w_vol_new, k_tend[:], diff_tend[:] = update_implicit(
                dt, dz, dz_h, d, k_flux, qground, w_vol, w_vol_new
            )

        # time step
        w_vol[:] = np.copy(w_vol_new)
        # for output:
        if options["discretization"] == "transformed":
            qground[:] = transform_back(qground_tr, metric_zmid)

        # Save output fields
        if timestep % out_interval == 0:
            output = save_timelevel(output, tlevel, outvars_dict)
            tlevel += 1

    # Save Grid Information
    output.update({"z": z})
    output.update({"z_mid": z_mid})
    t_ax = np.linspace(FloatWP(0.0), t_max, nout)
    output.update({"t_ax": t_ax})

    return output


if __name__ == "__main__":

    nt = int(time_params["t_max"] / time_params["dt"])
    infil = np.zeros(nt, dtype=FloatWP)
    infil[0 : int(nt / 4)] = FloatWP(10.0e-6)

    output = run_nolateral(infil)
    out_ds = prepare_netcdf(output)
    out_ds = plot_var_massgrid(output, out_ds, "qground")
    out_ds = plot_var_1d(output, out_ds, "wt_depth", "t_ax")
    out_ds = plot_var_1d(output, out_ds, "qsurf", "t_ax")
    out_ds = plot_var_fluxgrid(output, out_ds, "k_flux")
    out_ds = plot_var_massgrid(output, out_ds, "k_tend")
    out_ds = plot_var_fluxgrid(output, out_ds, "d_flux")
    out_ds = plot_var_massgrid(output, out_ds, "diff_tend")
    out_ds = plot_saturation(output, out_ds)
    out_ds.to_netcdf("test.nc")
