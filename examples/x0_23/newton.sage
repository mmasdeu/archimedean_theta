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

load('examples/x0_23/model.sage')
# unknowns y = (r1, c2, r2), c1 = 1; solve J(y) = TARGET by chord-Newton with increasing truncation order
def F(y, p, bits):
    R = RealField(bits)
    J, _ = invariants([1, y[1]], [y[0], y[2]], p, bits)
    return vector(R, [R(v.real()) for v in J]) - TARGET.change_ring(R)
def bits_for(p): return int(3.33*0.26*p) + 64
y = vector(RealField(bits_for(200)), [0.16675309408461045, 1.2883496775840126, 0.07913698863104919])
# Jacobian at p = 60 by central differences
p0 = 60; b0 = bits_for(p0); h = RealField(b0)(10)**-8
Jac = matrix(RealField(b0), 3, 3)
for k in range(3):
    e = vector(RealField(b0), [0,0,0]); e[k] = h
    Jac.set_column(k, (F(y+e, p0, b0) - F(y-e, p0, b0))/(2*h))
Jinv = Jac.inverse()
print('Jacobian computed', flush=True)
for p in [60, 120, 200]:
    b = bits_for(p); R = RealField(b); y = y.change_ring(R)
    for it in range(12):
        t0 = time.time(); Fy = F(y, p, b); y = y - Jinv.change_ring(R)*Fy
        res = float(max(abs(v) for v in Fy))
        print('p=%d it=%d residual %.2e time %.0f s  y=%s' % (p, it, res, time.time()-t0, [v.n(digits=25) for v in y]), flush=True)
        if res < 10**(-0.24*p + 1.5): break
    save(list(y), 'examples/x0_23/y_p%d.sobj' % p)
print('DONE', flush=True)
