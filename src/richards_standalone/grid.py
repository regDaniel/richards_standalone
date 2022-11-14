import numpy as np
import warnings as wa
#from parameter import z1, b, nz

#set precision globally
from parameter import precision_params
float_wp = precision_params['working_precision']



def setup_grid(grid_params):

    z1 = grid_params['z1']
    b  = grid_params['b']
    nz = grid_params['nz']

    #number of layer interfaces
    nz1 = nz + 1

    #vertical distance of new coord
    dzeta = np.log(b+1)

    #layer interfaces
    zeta = np.arange(0., nz1*dzeta, dzeta, dtype=float_wp)
    dz = np.zeros(nz, dtype=float_wp)
    z = np.zeros(nz1, dtype=float_wp)
    dz[0] = z1
    for k in range(1, nz1):
        z[k] = float_wp(1.)/b * z1 * (np.exp(zeta[k]) - 1)
    for k in range(1, nz):
        dz[k] = z[k+1] - z[k]

    #mass grid
    z_m = np.zeros(nz, dtype=float_wp)
    for k in range(0,nz):
        z_m[k] = float_wp(1.)/b * z1 * (np.exp((k+float_wp(0.5))*dzeta) - float_wp(1.))

    dz_h = np.zeros(nz1, dtype=float_wp)
    dz_h[0] = z_m[0] #this value will never be used if bdries are correct
    for k in range(1, nz):
        dz_h[k] = z_m[k] - z_m[k-1]
    dz_h[nz] = z[nz] - z_m[nz-1] #same, should not appear!

    if b < float_wp(1.E-3) :
        wa.warn('b is too small! For stability reasons a constant grid spacing\
        is used and the metric is set to 1.')
        dz[:] = z1
        dz_h[:] = z1
        for k in range(1, nz1):
            z[k] = z[k-1] + dz[k-1]
        for k in range(0,nz):
            z_m[k] = (k + 0.5) * z1



    return z, z_m, dz, dz_h

def calc_dzeta(grid_params):

    dzeta = np.log(grid_params['b'] + 1.)

    return dzeta


def calc_jacobian(z_array, grid_params):

    z1 = grid_params['z1']
    b  = grid_params['b']

    jacobian = float_wp(1.) / (z_array + z1 / b)
    if b < float_wp(1.E-3):
        jacobian[:] = float_wp(1.)

    return jacobian



def transform_scalar(x, jacobian):

    y = (1. / jacobian) * x

    return y


def transform_back(x, jacobian):

    y = jacobian * x

    return y
