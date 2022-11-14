# Third-party
import numpy as np
from parameter import precision_params

# set precision globally
float_wp = precision_params["working_precision"]


def init_output(outvars_dict, nout, nz):

    # create another dict with arrays for output.
    output = {
        entry: np.zeros((nout, len(outvars_dict[entry])), dtype=float_wp)
        for entry in outvars_dict.keys()
    }

    return output


def save_timelevel(output, tlevel, outvars_dict):

    for variable in output.keys():
        output[variable][tlevel, :] = outvars_dict[variable][:]

    return output
