"""Provide initial conditions from user settings."""
# Third-party
import numpy as np

# Local
from .parameter import grid_params
from .parameter import hydraulic_params
from .parameter import hydraulic_params_mvg
from .parameter import options
from .parameter import precision_params

FloatWP = precision_params["working_precision"]


def initial_condition(nsat_layers=0, sat_init=FloatWP(0.8), smooth_init_profile=False):

    nz = grid_params["nz"]
    nz1 = nz + 1
    if options["hydparam"] == "rjitema":
        pore_volume = hydraulic_params["pore_volume"]
    elif options["hydparam"] == "mvg":
        pore_volume = hydraulic_params_mvg["pore_volume"]

    # initial moisture profile
    w_vol = np.zeros(nz, dtype=FloatWP)
    w_vol[:] = sat_init * pore_volume
    if nsat_layers > 0:
        w_vol[(nz - nsat_layers) : nz] = pore_volume
    if smooth_init_profile:
        pass

    wt_depth = np.array([2.0], dtype=FloatWP)  # change accordingly
    qsurf = np.zeros(1, dtype=FloatWP)

    # arrays
    w_vol_new = np.zeros(nz, dtype=FloatWP)
    k = np.zeros(nz1, dtype=FloatWP)
    d = np.zeros(nz1, dtype=FloatWP)
    d_flux = np.zeros(nz1, dtype=FloatWP)
    qground = np.zeros(nz, dtype=FloatWP)

    return w_vol, w_vol_new, k, d, d_flux, qground, qsurf, wt_depth
