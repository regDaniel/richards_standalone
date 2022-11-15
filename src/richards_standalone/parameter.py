"""Control model parameter space."""
# Third-party
import numpy as np

precision_params = {
    "working_precision": np.double,
}

# set precision globally
FloatWP = precision_params["working_precision"]

time_params = {
    "dt": FloatWP(10.0),
    "t_max": FloatWP(3600.0 * 67.2),
    "dt_out": FloatWP(10.0),
}

grid_params = {
    "z1": FloatWP(0.01),
    "b": FloatWP(1.0e-1),
    "nz": int(32),
}

hydraulic_params = {
    "k0": FloatWP(531.0 * 10 ** (-8)),
    "k1": FloatWP(-19.88),
    "d0": FloatWP(357.0 * 10 ** (-8)),
    "d1": FloatWP(-7.44),
    "pore_volume": FloatWP(0.455),
    "air_dry_pt": FloatWP(0.035),
}

# SAND
# hydraulic_params = {
# 'k0'                : FloatWP(4.79E-5),
# 'k1'                : FloatWP(-19.27),
# 'd0'                : FloatWP(1.84E-5),
# 'd1'                : FloatWP(-8.45),
# 'pore_volume'       : FloatWP(0.364),
# 'air_dry_pt'        : FloatWP(0.012),
#   }

hydraulic_params_mvg = {
    "k0": FloatWP(2.89e-6),
    "pore_volume": FloatWP(0.43),
    "air_dry_pt": FloatWP(0.078),
    "alp_mvg": FloatWP(3.6),
    "m_mvg": FloatWP(0.359),
}

# SAND
# hydraulic_params_mvg = {
# 'k0'                : FloatWP(8.25E-5),
# 'pore_volume'       : FloatWP(0.43),
# 'air_dry_pt'        : FloatWP(0.045),
# 'alp_mvg'           : FloatWP(14.5),
# 'm_mvg'             : FloatWP(0.627)
#  }

plant_params = {"f": FloatWP(0.0), "z_root": FloatWP(1.0)}

runoff_params = {
    "s_oro": FloatWP(0.2),
    "lg1": FloatWP(0.4),
    "epsilon": np.finfo(FloatWP).eps,
}

options = {"hydparam": "rjitema", "discretization": "transformed"}


def update_parameter(update_dict):

    for u_key, _ in update_dict:
        if u_key in time_params:
            time_params[u_key] = update_dict[u_key]
        elif u_key in grid_params:
            grid_params[u_key] = update_dict[u_key]
        elif u_key in hydraulic_params:
            hydraulic_params[u_key] = update_dict[u_key]
        elif u_key in plant_params:
            plant_params[u_key] = update_dict[u_key]
        elif u_key in runoff_params:
            runoff_params[u_key] = update_dict[u_key]
        elif u_key in options:
            options[u_key] = update_dict[u_key]
        else:
            raise RuntimeError(
                f"Parameter key '{u_key}' could not be found." "Please check spelling"
            )

    return time_params, grid_params, hydraulic_params, runoff_params
