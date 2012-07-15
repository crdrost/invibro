# -*- coding: utf-8 -*-

"""
This file generates graphs for linearly-dependent coupling 1 + 0.5 x, but with a
very high temperature and a very clear investigation of the apparent falloff at
E ~= +/- h f.
"""

import invibro
from invibro.common import lorentzian, linear_coupling, thermal_dist, harmonic_energies
from numpy import exp, linspace
from datetime import datetime
import Gnuplot

invibro.dim = 8

def render(xs, ys, ref, name):
    g = Gnuplot.Gnuplot()
    #g('set xtics 1,0.5,3')
    g('set xrange [%s:%s]' % (xs[0], xs[-1]))
    g('set yrange [0:%s]' % max(ys))
    g('set size 0.6,0.6')
    g.xlabel('E (meV)')
    g.ylabel('A (1/meV)')
    g.plot(
        Gnuplot.Data(xs, ref, **{'with': 'lines lt 2 lc rgbcolor "blue" lw 2'}),
        Gnuplot.Data(xs, ys, **{'with': 'lines lt 1 lc rgbcolor "red" lw 2'})
    )
    g.hardcopy(name, color=1, terminal='postscript')
    print("File printed to %s" % name)

kelvin = 8.6173423e-5 * 1000 # 1K in meV
n = 1
for k in (0.6, -0.3):
    xs = linspace(-0.6, k, 2000)
    energies = harmonic_energies(hf=0.5)
    T = 5 * 0.5
    print("calc %s" % n); n += 1
    t0 = datetime.utcnow()
    calc = {
        'leads': [
            invibro.Lead(0.0, T, 0.025, linear_coupling(0.5), 50000.0),
            invibro.Lead(0.0, T, 0.025, linear_coupling(0.5), 50000.0)
        ],
        'e_level': 0.0,
        'ph_energy': energies,
        'ph_state': thermal_dist(T, energies)
    }
    ref = lorentzian(xs, calc)
    ys = invibro.density_of_states(xs, calc)
    render(xs, ys, ref, '/tmp/falloff-highT_%s.ps' % k)
    print("finished in: %s" % (datetime.utcnow() - t0))
    print("")
