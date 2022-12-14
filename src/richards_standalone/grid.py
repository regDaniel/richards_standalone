"""Functions to set up the grid for the numerical approximations."""
# Standard library
import warnings as wa

# Third-party
import numpy as np

# Local
# set precision globally
from .parameter import precision_params

# from parameter import z1, b, nz


FloatWP = precision_params["working_precision"]


def setup_grid(grid_params):

    z1 = grid_params["z1"]
    b = grid_params["b"]
    nz = grid_params["nz"]

    # number of layer interfaces
    nz1 = nz + 1

    # vertical distance of new coord
    dzeta = np.log(b + 1)

    # layer interfaces
    zeta = np.arange(0.0, nz1 * dzeta, dzeta, dtype=FloatWP)
    dz = np.zeros(nz, dtype=FloatWP)
    z = np.zeros(nz1, dtype=FloatWP)
    dz[0] = z1
    for k in range(1, nz1):
        z[k] = FloatWP(1.0) / b * z1 * (np.exp(zeta[k]) - 1)
    for k in range(1, nz):
        dz[k] = z[k + 1] - z[k]

    # mass grid
    z_m = np.zeros(nz, dtype=FloatWP)
    for k in range(0, nz):
        z_m[k] = (
            FloatWP(1.0) / b * z1 * (np.exp((k + FloatWP(0.5)) * dzeta) - FloatWP(1.0))
        )

    dz_h = np.zeros(nz1, dtype=FloatWP)
    dz_h[0] = z_m[0]  # this value will never be used if bdries are correct
    for k in range(1, nz):
        dz_h[k] = z_m[k] - z_m[k - 1]
    dz_h[nz] = z[nz] - z_m[nz - 1]  # same, should not appear!

    if b < FloatWP(1.0e-3):
        wa.warn(
            "b is too small! For stability reasons a constant grid spacing\
        is used and the metric is set to 1."
        )
        dz[:] = z1
        dz_h[:] = z1
        for k in range(1, nz1):
            z[k] = z[k - 1] + dz[k - 1]
        for k in range(0, nz):
            z_m[k] = (k + 0.5) * z1

    return z, z_m, dz, dz_h


def calc_dzeta(grid_params):

    dzeta = np.log(grid_params["b"] + 1.0)

    return dzeta


def calc_jacobian(z_array, grid_params):

    z1 = grid_params["z1"]
    b = grid_params["b"]

    jacobian = FloatWP(1.0) / (z_array + z1 / b)
    if b < FloatWP(1.0e-3):
        jacobian[:] = FloatWP(1.0)

    return jacobian


def transform_scalar(x, jacobian):

    y = (1.0 / jacobian) * x

    return y


def transform_back(x, jacobian):

    y = jacobian * x

    return y
