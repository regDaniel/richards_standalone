"""Numerical solvers."""
# Third-party
import numpy as np
from scipy.linalg import solve_banded

# Local
# set precision globally
from .parameter import precision_params

FloatWP = precision_params["working_precision"]


# pylint: disable=too-many-arguments
def update_euler(k, d, qground, w_vol, dz, dz_h, dt):

    nz = len(w_vol)
    j = slice(1, nz - 1)
    jp1 = slice(2, nz)
    jm1 = slice(0, nz - 2)

    flux_sum = np.zeros(nz, dtype=FloatWP)

    flux_sum[j] = (
        (
            d[jp1] * (w_vol[jp1] - w_vol[j]) / dz_h[jp1]
            - d[j] * (w_vol[j] - w_vol[jm1]) / dz_h[j]
        )
        / dz[j]
        - (k[jp1] - k[j]) / dz[j]
        - qground[j]
    )

    flux_sum[0] = (
        d[1] * (w_vol[1] - w_vol[0]) / (dz_h[1] * dz[0])
        - (k[1] - k[0]) / dz[0]
        - qground[0]
    )
    flux_sum[nz - 1] = (
        -d[nz - 1] * (w_vol[nz - 1] - w_vol[nz - 2]) / (dz_h[nz - 1] * dz[nz - 1])
        - (k[nz - 1] - k[nz - 2]) / dz[nz - 1]
        - qground[nz - 1]
    )

    w_vol_upd = w_vol + dt * flux_sum

    return w_vol_upd


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def update_implicit(dt, dz, dz_h, d, k, qground, w_vol, w_vol_new):

    # matrix dimensions
    nz = len(w_vol)
    len_diag = nz

    # Indices
    j = slice(0, nz)
    jp1 = slice(1, nz + 1)

    # matrix coefficients
    a = np.zeros(nz, dtype=FloatWP)
    b = np.zeros(nz, dtype=FloatWP)
    c = np.zeros(nz, dtype=FloatWP)

    b[j] = -dt / dz[j] * d[jp1] / dz_h[jp1]
    c[j] = -dt / dz[j] * d[j] / dz_h[j]
    a[j] = (dt / dz[j]) * (d[jp1] / dz_h[jp1] + d[j] / dz_h[j])
    # RHS
    y = np.zeros(nz, FloatWP)
    y[j] = w_vol[j] - dt / dz[j] * (k[jp1] - k[j]) - dt * qground[j]

    # diagonal and off-diag vectors
    diag_main = FloatWP(1.0) + a[j]
    diag_lower = c[(1):(nz)]
    diag_upper = b[: (nz - 1)]

    # set up array with matrix diagonals
    lhs_banded = np.zeros((3, len_diag), dtype=FloatWP)
    lhs_banded[0, 1:] = diag_upper
    lhs_banded[1, :] = diag_main
    lhs_banded[2, :-1] = diag_lower

    # solve lsoe
    w_vol_new[0:nz] = solve_banded((1, 1), lhs_banded, y[j], check_finite=True)
    diff_tend = w_vol_new - y
    k_tend = y - w_vol - dt * qground

    # update solution
    return w_vol_new, k_tend, diff_tend


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def update_implicit_tr(dt, dzeta, j_h, j_mid, d, k, w_vol, qground_tr, w_vol_new):

    # matrix dimensions
    nz = len(w_vol)
    len_diag = nz

    # Indices
    j = slice(0, nz)
    jp1 = slice(1, nz + 1)

    # matrix coefficients
    a = np.zeros(nz, dtype=FloatWP)
    b = np.zeros(nz, dtype=FloatWP)
    c = np.zeros(nz, dtype=FloatWP)

    b[j] = -dt / dzeta * d[jp1] / dzeta * j_h[jp1]
    c[j] = -dt / dzeta * d[j] / dzeta * j_h[j]
    a[j] = (dt / dzeta) * (d[jp1] / dzeta * j_h[jp1] + d[j] / dzeta * j_h[j])
    # RHS
    y = np.zeros(nz, FloatWP)
    y[j] = w_vol[j] / j_mid[j] - dt / dzeta * (k[jp1] - k[j]) - dt * qground_tr[j]

    # diagonal and off-diag vectors
    diag_main = FloatWP(1.0) / j_mid[j] + a[j]
    diag_lower = c[(1):(nz)]
    diag_upper = b[: (nz - 1)]

    # set up array with matrix diagonals
    lhs_banded = np.zeros((3, len_diag), dtype=FloatWP)
    lhs_banded[0, 1:] = diag_upper
    lhs_banded[1, :] = diag_main
    lhs_banded[2, :-1] = diag_lower

    # solve lsoe
    w_vol_new[0:nz] = solve_banded((1, 1), lhs_banded, y[j], check_finite=True)
    diff_tend = w_vol_new - y * j_mid
    k_tend = y * j_mid - w_vol - dt * qground_tr * j_mid

    # update solution
    return w_vol_new, k_tend, diff_tend
