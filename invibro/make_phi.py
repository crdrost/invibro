# -*- coding: utf-8 -*-
"""
This submodule is devoted to implementing the universal function `phi(x, Z)` 
described in my writeup. The function `phi` implements the spectral properties
of a self-energy in the wide-band limit for these sorts of modulated couplings.
These get stored in a base64-binary format in the file `cached.py` so that they
may be rapidly used for new calculations without recalculating the whole
function; we use scipy's interpolation functions to interpolate the values.
"""

from scipy.interpolate import interp1d
from scipy import integrate
from numpy import exp, cosh, sinh, pi, log, piecewise, vectorize
from invibro.utils import n, dndx, d2ndx2, d3ndx3, logspace, str_from_vec

@vectorize
def raw_phi0(x, Z0):
    """Calculate \phi_0(x) for some array of x."""
    n_x = n(x)
    # taylor expand n(z) about x to see that these are the correct
    # parameters to fill in the singularity near x:
    y0, m, k = -dndx(x), -d2ndx2(x)/2.0, -d3ndx3(x)/6.0
    def integrand(z):
        return (n(z) - n_x) / (x - z) if abs(x - z) > 1e-4 else \
            y0 + m * (z - x) + k * (z - x) ** 2
    return integrate.quad(integrand, -Z0, Z0)[0]

def make_phi0_cache(Z0, x_bound, spacing, logshape):
    """
    Cache values of \phi_0(x) for 0 <= x <= x_bound. 
    
    See `invibro.utils` for the docs on how the cache works and what
    `logshape` means. This cache carries an extra attribute `Z0` reflecting
    the value of Z0 handed in above.
    """
    xs = logspace(x_bound, spacing, logshape)
    ys = raw_phi0(xs, Z0)
    return {
        "Z0": Z0,
        "xs": (x_bound, spacing, logshape),
        "ys": ys,
        "interp": interp1d(xs, ys, 'linear')
    }

# Use cached values if possible.
#from invibro import cached
#phi0_cache = cached.phi0

try:
    from invibro import cached
    phi0_cache = cached.phi0
except ImportError:
    print("No cachefile, so creating our own cache. This may take a while.")
    phi0_cache = make_phi0_cache(200.0, 150.0, 0.001, 2.0)

def phi(xs, Z):
    """Evaluate phi(x) on a numpy list.
    
    This uses invibro.phi.phi0_cache to interpolate values."""
    Z0, x_bound = phi0_cache['Z0'], phi0_cache['xs'][0]
    c1 = pi ** 2 /6; c2 = 7 * pi ** 4 / 60.0
    neg = lambda x: log((Z - x) / (Z + x))
    interp = lambda x: phi0_cache['interp'](x) + log((Z + x) / (Z0 + x))
    large = lambda x: log(1 + Z / x) + c1 / x ** 2 + c2 / x ** 4
    return (
        piecewise(xs, [xs < 0], [neg, 0.0]) +
        piecewise(abs(xs), [abs(xs) < x_bound], [interp, large]) +
        complex(0, -0.5) * n(xs)
    )

def write_phi0_cache(out_name, cache):
    """Write a cache file similar to `invibro/cached.py`."""
    template = '''# -*- coding: utf-8 -*-
        from invibro import utils
        from scipy.interpolate import interp1d
        _Z0 = %s
        _xs = %s
        _ys = utils.vec_from_str("""%s""")
        phi0 = {
            "Z0": _Z0, "xs": _xs, "ys": _ys,
            "interp": interp1d(utils.logspace(*_xs), _ys, 'linear')
        }
        '''.replace("\n        ", "\n")
    with open(out_name, 'w') as f:
        f.write(template  % 
            (cache['Z0'], repr(cache['xs']), str_from_vec(cache['ys']))
        )
