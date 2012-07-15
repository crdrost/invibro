# -*- coding: utf-8 -*-

"""
This file generates graphs for linearly-dependent coupling 1 + 0.5 x, where the
ground state has an x-variance of 1.0 (i.e. x is dimensionless). 
"""

import invibro
from invibro.common import lorentzian, linear_coupling, thermal_dist, harmonic_energies
from numpy import exp, linspace
from datetime import datetime
import Gnuplot

invibro.dim = 30

def plot(xs, ys, ref, name):
    g = Gnuplot.Gnuplot()
    #g('set xtics 1,0.5,3')
    g('set size 0.6,0.6')
    if max(ref) > 10:
        g('set yrange [0:20]')
        g('set ytics 0,5,20')
    g.plot(
        Gnuplot.Data(xs, ref, **{'with': 'lines lt 2 lc rgbcolor "blue" lw 2'}),
        Gnuplot.Data(xs, ys, **{'with': 'lines lt 1 lc rgbcolor "red" lw 2'})
    )
    g.hardcopy(name, color=1, terminal='postscript')
    print("File printed to %s" % name)

kelvin = 8.6173423e-5 * 1000 # meV
n = 1
from matplotlib import pyplot
for k in (1.0, 0.1, 0.05):
    xs = linspace(-4.0, 4.0, 2000) if k == 1.0 else linspace(-2.0, 2.0, 2000)
    T = 5.2 * kelvin
    energies = harmonic_energies(hf=0.5)
    print("calc %s" % n); n += 1
    t0 = datetime.utcnow()
    calc = {
        'leads': [
            invibro.Lead(0.0, T, k, linear_coupling(0.5), 50000.0),
            invibro.Lead(0.0, T, k, linear_coupling(0.5), 50000.0)
        ],
        'e_level': 0.0,
        'ph_energy': energies,
        'ph_state': thermal_dist(T, energies)
    }
    ref = lorentzian(xs, calc)
    ys = invibro.density_of_states(xs, calc)
    plot(xs, ys, ref, '/tmp/linear_gamma=%s.ps' % k)
    print("finished in: %s" % (datetime.utcnow() - t0))
    print("")
