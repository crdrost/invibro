# -*- coding: utf-8 -*-
"""
This submodule implements certain auxiliary functions which were very helpful in
my calculated examples but which are not needed for the core logic implemented
in `invibro/__init__.py`.
"""

from decimal import Decimal
from numpy import exp

def delta(i, j):
    """Calculate matrix elements for an identity operator."""
    return 1 if i == j else 0

def annihilator(u, v):
    """
    Calculate matrix elements for the case of the annihilation operator:
        b |n> = sqrt(n) |n - 1>
    """
    return v ** 0.5 if u == v - 1 else 0.0

def x_quadrature(u, v):
    """Calculate matrix elements for the x-quadrature $b^\dagger + b$."""
    return annihilator(u, v) + annihilator(v, u)

def linear_coupling(k):
    """Calculate matrix elements for the case of linear coupling $1 + kx$."""
    return lambda u, v: delta(u, v) + k * x_quadrature(u, v)
    
def lorentzian(xs, params):
    """Calculate the Lorentzian curve for the current leads in the limit of no
    electron-phonon coupling."""
    gamma = sum(lead.gamma for lead in params['leads'])
    return gamma / ((xs - params['e_level'])**2 + (gamma/2)**2)

def harmonic_energies(hf=1.0):
    """Compute energy values for a normal harmonic oscillator; E_n = n * hf."""
    return lambda n: hf * (n + 0.5)
    
def thermal_dist(temp, energies):
    """
    The matrix for a thermal distribution with the given energies function. The
    matrix will be normalized later, so here it is presented as unnormalized.
    """
    return lambda u, v: exp(-energies(u) / temp) if u == v else 0.0

def fact(n, k=0):
    """Compute n factorial, and partial factorials n!/k! = n * n - 1 * ... * k + 1."""
    prod = 1
    for i in xrange(k + 1, n + 1):
        prod *= i
    return prod

def sqrt_fact(n, k=0):
    """Compute the square root of n! if 0 <= n < 300 to double precision."""
    return float(Decimal(fact(n, k)).sqrt())

def coherent_dist(z):
    """
    A matrix for a coherent state distribution with the given value b |z> = z |z>.
    """
    # convert to magnitude, phasor representation so the result is always real.
    mag, phasor = abs(z), z / abs(z)
    norm = exp(mag**2)
    return lambda m, n: mag ** (m + n) * phasor ** (m - n) / (
        norm * sqrt_fact(m) * sqrt_fact(n)
    )