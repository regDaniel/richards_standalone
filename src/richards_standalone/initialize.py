# Third-party
import numpy as np
import parameter as p
# set precision globally
from parameter import precision_params

float_wp = precision_params["working_precision"]


def initial_condition(nsat_layers=0, sat_init=float_wp(0.8), smooth_init_profile=False):

    nz = p.grid_params["nz"]
    nz1 = nz + 1
    pore_volume = p.hydraulic_params["pore_volume"]

    # initial moisture profile
    w_vol = np.zeros(nz, dtype=float_wp)
    w_vol[:] = sat_init * pore_volume
    if nsat_layers > 0:
        w_vol[(nz - nsat_layers) : nz] = pore_volume
    if smooth_init_profile:
        pass

    wt_depth = np.array([2.0], dtype=float_wp)  # change accordingly
    qsurf = np.zeros(1, dtype=float_wp)

    # arrays
    w_vol_new = np.zeros(nz, dtype=float_wp)
    k = np.zeros(nz1, dtype=float_wp)
    d = np.zeros(nz1, dtype=float_wp)
    d_flux = np.zeros(nz1, dtype=float_wp)
    qground = np.zeros(nz, dtype=float_wp)

    return w_vol, w_vol_new, k, d, d_flux, qground, qsurf, wt_depth
