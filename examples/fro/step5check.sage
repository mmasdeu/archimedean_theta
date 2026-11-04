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
phis = load('examples/fro/phis_240.sobj'); R0 = phis[0].parent()
x = PolynomialRing(R0, 'x').gen(); f = prod(x - ph for ph in phis)
I2, I4, I6, I10 = igusa_clebsch_invariants(f)
D = 100; R = RealField(int(D*3.33)); PI = R.pi()
v = R(I2**5/I10); w = v/10**floor(log(abs(v),10))
xs = [w] + [PI**j for j in range(3)]
print('pari realprecision default:', pari.default('realprecision'))
rel = pari(xs).lindep()
print('relation', rel, ' residual', sum(R(ZZ(c))*t for c,t in zip(rel,xs)))
pari.set_real_precision(D)
rel = pari(xs).lindep()
print('with realprecision', D, ':', rel, ' residual', sum(R(ZZ(c))*t for c,t in zip(rel,xs)), ' log10 height', float(log(max(abs(ZZ(c)) for c in rel),10)))
