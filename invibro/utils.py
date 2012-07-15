# -*- coding: utf-8 -*-
"""
This submodule implements utilities used within invibro itself which are not
necessarily intended to be exposed to other programs (see `common.py` for the
ones which actively are). These include calculating the dimensionless
Fermi-Dirac function and its derivatives; generating integration lattices which
can scale logarithmically with the integration bounds rather than linearly; and
serializing numpy floats into base64 strings for inclusion in .py files.
"""

from numpy import linspace, array, log, exp, sinh, cosh, vectorize
from struct import pack, unpack
from base64 import b64encode, b64decode
from zlib import compress, decompress

@vectorize
def n(x):
    """Calculate the Fermi-Dirac occupation function."""
    return 1.0 if x.real < -40 else 0.0 if x.real > 40 else \
            0.5 * exp(-0.5 * x.real) / cosh(0.5 * x.real)

n.__doc__ = """Calculate the Fermi-Dirac occupation function."""

def dndx(x):
    """First derivative of the Fermi-Dirac function"""
    return -0.25 / cosh(x / 2.0) ** 2 if -40 < x.real < 40 else 0.0

def d2ndx2(x):
    """Second derivative of the Fermi-Dirac function"""
    return 0.25 * sinh(x / 2.0) / cosh(x / 2.0) ** 3 if -40 < x.real < 40 else 0.0

def d3ndx3(x):
    """Third derivative of the Fermi-Dirac function"""
    return 0.125 * (2.0 - cosh(x)) / cosh(x / 2.0) ** 4 if -40 < x.real < 40 else 0.0


def logspace(bound, spacing, shape):
    """
    Generate a lattice of floats `x` with density function:

        f(x) ~ 1 / (spacing * (1 + x/shape))

    This means that for values less than `shape` it stays mostly evenly spaced
    (there is some linear dependence but it's not too bad) but at 10,000 you
    see 1/10th as many points as you did at 1,000, and so forth; the density
    goes like 1/x.

    The number of points `N` in this lattice grows asymptotically like:

        N ~= shape * ln(Z) / spacing

    This is why I have called it a logspace; going to a larger `Z` 
    only takes time logarithmic in `Z`. It is more efficient than linspace for
    the same error if your second derivative falls like 1/x^2.

    The exact formula was computed by converting f(x) on [0, Z] into the
    cumulative distribution function `F(x) = log(1 + x/shape)/log(1 + Z/shape)`
    and inverting this into `x = shape * ((1 + Z/shape)^F - 1)`. We sample this
    evenly on the interval `[0, 1]`. The distance between `0` and `1/N` gives a
    formula for `N` in terms of `spacing`.
    """
    shape = float(shape) # force float conversion.
    N = int(log(1 + bound / shape) / log(1 + spacing / shape)) + 1
    points = linspace(0.0, 1.0, N)
    return shape * (((1.0 + bound/shape) ** points) - 1)

def str_from_vec(fv):
    """Convert a vector of numpy floats into a string."""
    s = b64encode(compress(pack("<%sd" % len(fv), *fv), 9))
    out = "\n"
    for i in range(0, len(s), 80):
        out += s[i : i + 80] + "\n"
    return out

def vec_from_str(s):
    """Convert a string made by str_from_vec into a numpy vector."""
    base = decompress(b64decode(s.replace("\n", "")))
    return array(unpack("<%sd" % (len(base) / 8), base))