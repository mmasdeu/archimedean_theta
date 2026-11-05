# Copyright (C) 2026 Marc Masdeu
#
# This file is part of archimedean_theta.
#
# archimedean_theta is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# archimedean_theta is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with archimedean_theta.  If not, see <https://www.gnu.org/licenses/>.

# Truncation error of the period matrix of the ring family (Section 7.4), direct solver at every p.
# usage: sage benchmarks/rates_ring.sage g
import os, sys, time
sys.path.insert(0, os.getcwd())
from schottky_theta import *
g = int(sys.argv[-1])
K = ComplexField(260)
def ring(g, K):
    eta = K.pi()/(6*g); theta = K.pi()/(2*g); I = K(0,1)
    z1 = (I*((theta+eta)/2)).exp()/((theta-eta)/2).cos(); r = ((theta-eta)/2).tan()
    D = {}
    for j in range(g):
        D[j+1] = Disc((2*K.pi()*I*j/g).exp()*z1.conjugate(), r); D[-(j+1)] = Disc((2*K.pi()*I*j/g).exp()*z1, r)
    return D
D = ring(g, K)
P = {}
for p in [10, 20, 30, 40, 60, 90]:
    t0 = time.time(); S = SchottkyTheta(D, p); F = S.solve(); P[p] = S.period_matrix(F)
    print('g=%d p=%d done (%.0f s)' % (g, p, time.time()-t0), flush=True)
ref = P[90]
for p in [10, 20, 30, 40, 60]:
    print('g=%d  p=%d  max |B(p) - B(90)| = %.1e' % (g, p, float(max(abs(P[p][a,b]-ref[a,b]) for a in range(g) for b in range(g)))), flush=True)
