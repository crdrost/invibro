# -*- coding: utf-8 -*-

"""
This file generates graphs which calculate the parameters in Flensberg's paper
"Tunneling broadening of vibrational sidebands in molecular transistors",
published in Phys. Rev. B 68 205323 (2003). However, since the calculation
scheme we've derived ultimately originates from this paper, the parameters
below instead calculate the case for Galperin, Nitzan, and Ratner's paper
"Resonant inelastic tunneling in molecular junctions", PRB 73: 045314 (2006).
"""

import invibro
from invibro.common import thermal_dist, harmonic_energies, fact
from numpy import exp, linspace, pi

def choose(n, k):
    """Choose function, n choose k."""
    return fact(n, n - k) / fact(k)

def laguerre(n, a, x):
    """Laguerre polynomial L_n^(a)(x) as defined on A&S p. 775."""
    return sum(choose(n + a, n - m) * (-x) ** m / fact(m) for m in range(0, n + 1))

def displace(L):
    g = L ** 2
    def sgn(x):
        return -1 if x < 0 else 1
    def matrix(m, n, L):
        if m > n: 
            return matrix(n, m, -L)
        else:
            return (exp(-g) * g ** (n - m) * fact(m)/fact(n)) ** 0.5 * \
                laguerre(m, n - m, g) * sgn(L) ** (n - m)
    return lambda m, n: matrix(m, n, L)

w0 = 0.02 # eV
T = 0.0258520269 # 300K in eV
invibro.dim=8

energies = harmonic_energies(0.02)

q = invibro.mat_from_fn(displace(1.0))
qbar = q.T.conjugate()

rho_0 = invibro.mat_from_fn(thermal_dist(T, energies))
rho_0 /= rho_0.trace()
rho_1 = qbar.dot(rho_0).dot(q)

from datetime import datetime
import Gnuplot
def plot(xs, ys, name, first=False):
    g = Gnuplot.Gnuplot()
    g('set xrange [1.88:2.08]')
    g('set yrange [0:200]')
    g('set xtics 1.92,0.04,2.04')
    g('set size 0.5,0.5')
    if first:
        g('set ytics 0,40,320')
    else:
        g('set ytics -80,400,320')
        g('set mytics 10')
    g.plot(Gnuplot.Data(xs, ys / (2 * pi), **{'with': 'lines lc rgbcolor "blue" lw 3'}))
    g.hardcopy(name, color=1, terminal='postscript')
    print("File printed to %s" % name)

matrices = (rho_0, 0.5 * rho_0 + 0.5 * rho_1, rho_1)
fillings = (0.0, 0.5, 1.0)
for n in range(len(matrices)):
    print "\ncalc %i" % (n + 1);
    t0 = datetime.utcnow()
    xs = linspace(1.88, 2.08, 1000)
    ys = invibro.density_of_states(xs, {
        'postprocess': qbar,
        'leads': [
            invibro.Lead(-0.02, T, 0.002, q, 50.0),
            invibro.Lead(-0.02, T, 0.002, q, 50.0)
        ],
        'e_level': 1.98,
        'ph_energy': energies,
        'ph_state': q.dot(matrices[n])
    })
    plot(xs, ys, '/tmp/gnr_filling=%s.ps' % fillings[n], first=(n == 2))
    print "finished in: %s" % (datetime.utcnow() - t0)
