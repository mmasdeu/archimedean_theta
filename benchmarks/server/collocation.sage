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

# Collocation (Crowdy--Marshall 2007) versus the algorithm of this repository, on the square X = omega^2 of the
# Schottky--Klein prime function of the unit disc with N discs of radius 1/10 removed, centred at the N-th roots of
# unity times 1/2.  Usage (from the repository root):  sage benchmarks/server/collocation.sage N
import sys, os, time
sys.path.insert(0, os.getcwd())
from theta_arch import *
from collocation import Domain, prime_square

N = int(sys.argv[1])
PREC_C = 200
CC = ComplexField(400); I = CC.gen()
centers = [CC(1)/2 * (2*CC(pi)*I*k/N).exp() for k in range(N)]
balls = {k + 1: Ball(centers[k], CC(1)/10) for k in range(N)}
z = (1 + I)/10; w = CC(2)/10 - I/10
PREF = int(sys.argv[2]) if len(sys.argv) > 2 else 45
t0 = time.time(); ref = SKprime(balls, z, w, PREF)
print('N=%d reference SKprime p=%d time=%.1f X=%s' % (N, PREF, time.time() - t0, ref.n(digits=30)), flush=True)
for p in [10, 20, 30, 40]:
    t0 = time.time(); X = SKprime(balls, z, w, p); dt = time.time() - t0
    print('N=%d new p=%d time=%.2f err=%.1e' % (N, p, dt, float(abs(X - ref))), flush=True)
K = ComplexBallField(PREC_C)
Kc = lambda x: K(x.real(), x.imag())
D = Domain([Kc(c) for c in centers], [K(1)/10] * N)
for N1 in [10, 20, 30, 40]:
    t0 = time.time()
    X = prime_square(D, N1, PREC_C, [Kc(z)], Kc(w))[0]
    dt = time.time() - t0
    Xc = ComplexField(400)(X.mid())
    print('N=%d colloc N1=%d time=%.2f err=%.1e' % (N, N1, dt, float(abs(Xc - ref))), flush=True)
