# -*- coding: utf-8 -*-
"""
This module contains the core logic derived in my Master's thesis for systems of
vibrations interacting with a single electron level by directly modulating the
hopping amplitudes. In this file I take a Hungarian notation with dim x dim
matrices denoted by the prefix `mat_`, dim^2 x dim^2 matrices denoted by the
prefix `smat_`, and 1xdim^2 vectors denoted by the prefix `vec_`. 
"""

import invibro.make_phi
import invibro.utils
import numpy
from numpy import log, pi, exp

"""The universal wideband phi function defined in my writeup."""
phi = invibro.make_phi.phi

"""A convenience class for lumping together several parameters which describe
a lead."""
class Lead:
    def __init__(self, mu, T, gamma, vib, band_size):
        """Create a lead. parameters are chemical potential, temperature,
        coupling rate $\Gamma_r$, vibration matrix $f^r_{mn}$, and band width W_r."""
        self.phi = lambda omega: \
            (phi((omega - mu) / T, band_size / T) - log(1 - mu / (band_size + omega))) / (2 * pi)
        self.mat_f = mat_from_fn(vib)
        self.mat_fbar = self.mat_f.T.conjugate()
        self.gamma = gamma
        self.im0 = 0.5j * self.gamma * self.mat_f.dot(self.mat_fbar)

"""The number of phonon levels currently being tracked. If this is changed, you
should cancel all calculations and recreate the Lead objects."""
dim = 20

"""An example params dict for calculations. This is also authoritative: if you
don't explicitly provide your own params dict, this one will be used instead."""
params = {
    'postprocess': None,
    'leads': [],
    'e_level': 0.0,
    'ph_energy': lambda n: n,
    'ph_state': lambda m, n: 1.0 if m == n == 0 else 0.0
}

def mat_from_fn(fn):
    """Produce a (dim x dim) matrix from a function which gives its indices.

    If this is fed a numpy array it will simply copy that array, so it is always
    safe to run on things which 'might be' functions."""
    if type(fn) == numpy.ndarray:
        return fn.copy()
    else:
        out = numpy.zeros((dim, dim), dtype=complex)
        for u in range(0, dim):
            for v in range(0, dim):
                out[u][v] = fn(u, v)
        return out

def smat_Y(omega, mat_eps, leads):
    """Calculate the self-energy tensor."""
    cached = {}
    out = numpy.zeros((dim**2, dim**2), dtype=complex)
    for m in range(0, dim):
        for n in range(0, dim):
            e = mat_eps[m][n]
            if e not in cached:
                Y0, Y1 = 0.0, 0.0
                for lead in leads:
                    # we multiply f by phi and fbar by gamma.
                    aug_f = lead.gamma * lead.mat_f * lead.phi(omega + mat_eps - e)
                    Y0 -= aug_f.dot(lead.mat_fbar) + lead.im0
                    Y1 += lead.mat_fbar.dot(aug_f)
                cached[e] = (Y0, Y1)
            Y = cached[e]
            # Y[m][n][a][b] = delta[b][n] Y0[m][a] + delta[a][m] Y1[b][n]
            for a in range(0, dim):
                out[dim*m + n][dim*a + n] += Y[0][m][a]
            for b in range(0, dim):
                out[dim*m + n][dim*m + b] += Y[1][b][n]
    return out

def tensor_to_matrix(fn, omega):
    return numpy.array(tuple(tuple(
            # upper indices (a, b) increment first if Y^{ab}_{mn} G_{ab} is the product.
            fn(v / dim, v % dim, u / dim, u % dim, omega)
            for v in range(0, dim ** 2))
        for u in range(0, dim ** 2)), dtype=complex)

def vec_from_mat(fn):
    return mat_from_fn(fn).reshape(dim**2)

def mat_from_vec(v):
    return v.reshape((dim, dim))

def mat_eps(ph_energy):
    return mat_from_fn(lambda m, n: float(ph_energy(m) - ph_energy(n)))

def density_of_states(omegas, params=params):
    """Calculate the density of states for the dot level."""
    mat_E = mat_eps(params['ph_energy'])
    smat_E = numpy.diag(mat_E.reshape(dim ** 2))
    e0 = params['e_level']
    if params.get('postprocess') != None:
        mat_pp = mat_from_fn(params['postprocess'])
    else:
        mat_pp = numpy.diag(numpy.ones(dim))
    rho_b = mat_from_fn(params['ph_state'])
    rho_b /= rho_b.trace()
    @numpy.vectorize
    def calc(w):
        smat_Y_w = smat_Y(w, mat_E, params['leads'])
        mat_A = mat_from_vec(numpy.linalg.solve(
            numpy.diag((w - e0) * numpy.ones(dim**2)) - smat_E - smat_Y_w,
            vec_from_mat(rho_b)
        ))
        return -2 * (mat_pp.dot(mat_A).trace().imag)
    return calc(omegas)

