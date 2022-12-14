"""Functions to calculate the fluxes and flux corrections.

Inside the soil matrix and for the interfaces atmosphere - soil and
soil water - discharge.

"""
# Third-party
import numpy as np

# Local
from .grid import transform_scalar
from .parameter import precision_params

# set precision globally
FloatWP = precision_params["working_precision"]


def calc_k_decharme(z_array, k_sat, f, z_root):

    k0 = k_sat * np.exp(-f * (z_array - z_root))

    return k0


def hydraulics_rjitema(w_vol, k0, k, d, hydraulic_params):

    d0 = hydraulic_params["d0"]
    k1 = hydraulic_params["k1"]
    d1 = hydraulic_params["d1"]
    povo = hydraulic_params["pore_volume"]
    adp = hydraulic_params["air_dry_pt"]

    nz = len(w_vol)

    j = slice(1, nz)
    jm1 = slice(0, nz - 1)

    if not (len(k0) == nz + 1 or len(k0) == 1):
        raise RuntimeError(
            "k0 must either be an array of length nz+1 \
                            or a scalar"
        )
    w_vol_interpol = np.zeros(nz + 1, dtype=FloatWP)
    w_vol_interpol[j] = FloatWP(0.5) * (w_vol[j] + w_vol[jm1])
    w_vol_interpol[0] = w_vol[0]

    j = slice(0, nz)
    k[j] = k0[j] * np.exp(k1 * (povo - w_vol_interpol[j]) / (povo - adp))
    d[j] = d0 * np.exp(d1 * (povo - w_vol_interpol[j]) / (povo - adp))
    d[0] = FloatWP(0.0)  # boundary condition
    k[0] = k0[0]  # boundary condition
    d[-1] = FloatWP(0.0)
    k[-1] = FloatWP(0.0)

    return k, d


def hydraulics_mvg(w_vol, k0, k, d, hydraulic_params_mvg):

    povo = hydraulic_params_mvg["pore_volume"]
    adp = hydraulic_params_mvg["air_dry_pt"]
    alp_mvg = hydraulic_params_mvg["alp_mvg"]
    m_mvg = hydraulic_params_mvg["m_mvg"]

    nz = len(w_vol)
    j = slice(1, nz)
    jm1 = slice(0, nz - 1)

    if not (len(k0) == nz + 1 or len(k0) == 1):
        raise RuntimeError(
            "k0 must either be an array of length nz+1 \
                            or a scalar"
        )
    w_vol_interpol = np.zeros(nz + 1, dtype=FloatWP)
    w_vol_interpol[j] = FloatWP(0.5) * (w_vol[j] + w_vol[jm1])
    w_vol_interpol[0] = w_vol[0]

    rh = (w_vol_interpol - adp) / (povo - adp)  # relative humidity
    k[j] = (
        k0[j]
        * np.sqrt(rh[j])
        * (
            FloatWP(1.0)
            - np.abs((1.0) - (np.abs(rh[j])) ** (FloatWP(1.0) / m_mvg)) ** m_mvg
        )
        ** FloatWP(2.0)
    )
    d[j] = (
        (FloatWP(1.0) - m_mvg)
        * k[j]
        / (m_mvg * alp_mvg)
        * FloatWP(1.0)
        / (w_vol_interpol[j] - adp)
        * (np.abs(rh[j])) ** (-FloatWP(1.0) / m_mvg)
        * (np.abs((np.abs(rh[j])) ** (-FloatWP(1.0) / m_mvg) - FloatWP(1.0)))
        ** (-m_mvg)
    )

    # limit diffusion, in case of NaN it will drop to max value.
    diffmax = FloatWP(10.0 * 5.31e-6)
    d = np.nan_to_num(d, diffmax)
    d = np.minimum(d, diffmax)

    # boundary condition
    d[0] = FloatWP(0.0)
    k[0] = k0[0]

    return k, d


def est_flux_fg(k, d, dflux_fg, w_vol, dz_h):

    nz = len(w_vol)
    # the last layer interface must be set to the boundary condition.
    j = slice(0, nz - 1)
    jp1 = slice(1, nz)

    dflux_fg[jp1] = -d[jp1] * FloatWP(1.0) / dz_h[jp1] * (w_vol[jp1] - w_vol[j])
    kflux_fg = k

    dflux_fg[-1] = FloatWP(0.0)
    kflux_fg[-1] = FloatWP(0.0)  # boundary condition

    return dflux_fg, kflux_fg


def flux_fg_transformed(k, d, dflux_fg, w_vol, dzeta):

    nz = len(w_vol)
    # the last layer interface must be set to the boundary condition.
    j = slice(0, nz - 1)
    jp1 = slice(1, nz)

    dflux_fg[jp1] = -d[jp1] * FloatWP(1.0) / dzeta * (w_vol[jp1] - w_vol[j])
    kflux_fg = k

    dflux_fg[-1] = FloatWP(0.0)
    kflux_fg[-1] = FloatWP(0.0)  # boundary condition

    return dflux_fg, kflux_fg


# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
def est_qground_fg(qground, w_vol, dz, z, k0_m, runoff_params, hydraulic_params):

    nz = len(w_vol)
    qground[:] = FloatWP(0.0)  # reset

    pore_volume = hydraulic_params["pore_volume"]
    s_oro = runoff_params["s_oro"]
    lg1 = runoff_params["lg1"]
    eps = runoff_params["epsilon"]

    j_unsat = nz - 1
    while (w_vol[j_unsat] >= pore_volume - eps) and j_unsat >= 0:
        j_unsat = j_unsat - 1

    j_sat = j_unsat + 1
    if j_sat == 0:
        j_unsat = j_sat

    if j_sat > 0:
        d1 = min(
            FloatWP(1.0),
            max(FloatWP(0.0), (w_vol[j_unsat] - w_vol[j_unsat - 1]))
            / (pore_volume - w_vol[j_unsat - 1] + eps),
        )
        wt_depth = z[j_sat] - d1 * dz[j_unsat]
    else:
        wt_depth = FloatWP(0.0)

    # discharge from saturated layers
    qground[j_sat:nz] = s_oro * lg1 * k0_m[j_sat:nz]
    if j_unsat > 0:
        qground[j_unsat] = s_oro * lg1 * k0_m[j_unsat] * d1

    return qground, wt_depth


# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
def flux_corr(
    k_flux,
    d,
    infil,
    qground,
    w_vol,
    dz,
    time_params,
    hydraulic_params,
):

    pore_volume = hydraulic_params["pore_volume"]
    air_dry_pt = hydraulic_params["air_dry_pt"]
    dt = time_params["dt"]
    nz = len(qground)

    # overdepletion
    for j in range(nz - 1, -1, -1):
        alpha_out = max(
            FloatWP(0.0),
            min(
                FloatWP(1.0),
                (w_vol[j] - air_dry_pt) / (dt * (k_flux[j + 1] / dz[j] + qground[j])),
            ),
        )
        k_flux[j + 1] = k_flux[j + 1] * alpha_out
        qground[j] = qground[j] * alpha_out

    # overfill
    for j in range(nz - 1, -1, -1):
        k_limit = k_flux[j + 1] + dz[j] * ((pore_volume - w_vol[j]) / dt + qground[j])
        k_flux[j] = max(FloatWP(0.0), min(k_flux[j], k_limit))

    # surface runoff
    qsurf = max(FloatWP(0.0), infil - k_flux[0])
    k_flux[0] = min(k_flux[0], infil)

    return k_flux, d, qground, qsurf


# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
def flx_corr_full(
    k_flux,
    d_flux,
    d,
    infil,
    qground,
    w_vol,
    dz,
    time_params,
    hydraulic_params,
    runoff_params,
):

    pore_volume = hydraulic_params["pore_volume"]
    air_dry_pt = hydraulic_params["air_dry_pt"]
    eps = runoff_params["epsilon"]
    dt = time_params["dt"]
    nz = len(qground)

    # overdepletion
    for j in range(nz - 1, -1, -1):
        flux_out = (
            k_flux[j + 1]
            + max(FloatWP(0.0), d_flux[j + 1])
            - min(FloatWP(0.0), d_flux[j])
        )
        alpha_out = max(
            FloatWP(0.0),
            min(
                FloatWP(1.0),
                (w_vol[j] - air_dry_pt) / (dt * (flux_out / dz[j] + qground[j]) + eps),
            ),
        )
        k_flux[j + 1] = k_flux[j + 1] * alpha_out
        qground[j] = qground[j] * alpha_out

    # overfill
    for j in range(nz - 1, -1, -1):
        flux_out2 = k_flux[j + 1] + d_flux[j + 1] - d_flux[j]
        k_limit = flux_out2 + dz[j] * ((pore_volume - w_vol[j]) / dt + qground[j])
        k_flux[j] = max(FloatWP(0.0), min(k_flux[j], k_limit + eps))

    # surface runoff
    qsurf = max(FloatWP(0.0), infil - k_flux[0])
    k_flux[0] = min(k_flux[0], infil)

    return k_flux, d, qground, qsurf


# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
def flux_corr_transformed(
    k_flux,
    d,
    infil,
    qground_tr,
    w_vol_tr,
    dzeta,
    time_params,
    hydraulic_params,
    runoff_params,
    metric_zmid,
):

    pore_volume_tr = transform_scalar(hydraulic_params["pore_volume"], metric_zmid)
    air_dry_pt_tr = transform_scalar(hydraulic_params["air_dry_pt"], metric_zmid)
    dt = time_params["dt"]
    nz = len(qground_tr)

    # overdepletion
    for j in range(nz - 1, -1, -1):
        alpha_out = max(
            FloatWP(0.0),
            min(
                FloatWP(1.0),
                (w_vol_tr[j] - air_dry_pt_tr[j])
                / (dt * (k_flux[j + 1] / dzeta + qground_tr[j] + runoff_params["eps"])),
            ),
        )
        k_flux[j + 1] = k_flux[j + 1] * alpha_out
        qground_tr[j] = qground_tr[j] * alpha_out

    # overfill
    for j in range(nz - 1, -1, -1):
        k_limit = k_flux[j + 1] + dzeta * (
            (pore_volume_tr[j] - w_vol_tr[j]) / dt + qground_tr[j]
        )
        k_flux[j] = max(FloatWP(0.0), min(k_flux[j], k_limit))

    # surface runoff
    qsurf = max(FloatWP(0.0), infil - k_flux[0])
    k_flux[0] = min(k_flux[0], infil)

    return k_flux, d, qground_tr, qsurf


def overflow_clip(w_vol, qground, hydraulic_params, time_params):
    """Clip water content to a max of pore volume."""
    pore_volume = hydraulic_params["pore_volume"]

    qground = (
        qground + np.maximum(FloatWP(0.0), w_vol - pore_volume) / time_params["dt"]
    )
    w_vol = np.minimum(w_vol, pore_volume)

    return w_vol, qground


def underflow_clip(w_vol, hydraulic_params):
    """Clip water content to a min of air dry pt."""
    air_dry_pt = hydraulic_params["air_dry_pt"]

    w_vol = np.maximum(w_vol, air_dry_pt)

    return w_vol
