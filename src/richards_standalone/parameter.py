# Third-party
import numpy as np

precision_params = {
    "working_precision": np.double,
}

# set precision globally
float_wp = precision_params["working_precision"]

time_params = {
    "dt": float_wp(10.0),
    "t_max": float_wp(3600.0 * 67.2),
    "dt_out": float_wp(10.0),
}

grid_params = {
    "z1": float_wp(0.01),
    "b": float_wp(1.0e-1),
    "nz": int(32),
}

hydraulic_params = {
    "k0": float_wp(531.0 * 10 ** (-8)),
    "k1": float_wp(-19.88),
    "d0": float_wp(357.0 * 10 ** (-8)),
    "d1": float_wp(-7.44),
    "pore_volume": float_wp(0.455),
    "air_dry_pt": float_wp(0.035),
}

# SAND
# hydraulic_params = {
# 'k0'                : float_wp(4.79E-5),
# 'k1'                : float_wp(-19.27),
# 'd0'                : float_wp(1.84E-5),
# 'd1'                : float_wp(-8.45),
# 'pore_volume'       : float_wp(0.364),
# 'air_dry_pt'        : float_wp(0.012),
#   }

hydraulic_params_mvg = {
    "k0": float_wp(2.89e-6),
    "pore_volume": float_wp(0.43),
    "air_dry_pt": float_wp(0.078),
    "alp_mvg": float_wp(3.6),
    "m_mvg": float_wp(0.359),
}

# SAND
# hydraulic_params_mvg = {
# 'k0'                : float_wp(8.25E-5),
# 'pore_volume'       : float_wp(0.43),
# 'air_dry_pt'        : float_wp(0.045),
# 'alp_mvg'           : float_wp(14.5),
# 'm_mvg'             : float_wp(0.627)
#  }

plant_params = {"f": float_wp(0.0), "z_root": float_wp(1.0)}

runoff_params = {
    "s_oro": float_wp(0.2),
    "lg1": float_wp(0.4),
    "epsilon": np.finfo(float_wp).eps,
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
