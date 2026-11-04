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

from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
def invs(p):
    phis = load('examples/fro/phis_%d.sobj' % p); R = phis[0].parent()
    x = PolynomialRing(R,'x').gen(); f = prod(x - ph for ph in phis)
    I2, I4, I6, I10 = igusa_clebsch_invariants(f)
    return [I2**5/I10, I2**3*I4/I10, I2**2*I6/I10]
a = invs(200); b = invs(240)
for k in range(3):
    print('i%d: agreement p200/p240 = %.1f digits' % (k+1, float(-log(abs((a[k]-b[k])/b[k]), 10))))
D = 100
R = RealField(int(D*3.33))
for k, v in enumerate(b):
    v = R(v); m = floor(log(abs(v),10)); v = v/10**m
    print('i%d / 10^%d = %s' % (k+1, m, v.n(digits=40)))
    for d in range(1, 9):
        P = algdep(v, d)
        H = max(abs(c) for c in P.coefficients())
        print('   degree %d: best polynomial height 10^%.1f   (generic expectation 10^%.1f), residual %.1e' % (d, float(log(H,10)), D/(d+1), float(abs(P(v)))))
