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

# Independent check (run with the system sage, which has mpmath): the Riemann matrix of the Schottky group
# (Algorithm 2, period differences) must be symplectically equivalent to the Riemann matrix of the curve
# y^2 = f(x) computed by numerical integration (Sage's RiemannSurface).
load('examples/x0_23/model.sage')
from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface, integer_matrix_relations
import sys
p = int(sys.argv[-1]); bits = int(3.33*0.26*p) + 64
y = load('examples/x0_23/y_p%d.sobj' % p)
K = ComplexField(bits)
discs, gens = datum([1, y[1]], [y[0], y[2]], K)
S = SchottkyTheta(discs, p, gens); F = S.solve()
B = S.period_matrix(F)
tau_s = B/(2*K.pi()*K.gen())
# entries of tau are defined modulo 1 (logarithms); make the matrix symmetric by an integer shift
tau_s[0,1] = tau_s[0,1] + (tau_s[1,0] - tau_s[0,1]).real().round()
print('Schottky tau (p=%d):' % p); print(tau_s.apply_map(lambda z: z.n(digits=20)))
print('symmetry defect %.1e' % float(max(abs(B[a,b]-B[b,a]) for a in range(2) for b in range(2))))
A2.<xx, yy> = QQ[]
f = (xx^3 - xx + 1)*(xx^3 - 8*xx^2 + 3*xx - 7)
RS = RiemannSurface(yy^2 - f, prec=200)
tau_c = RS.riemann_matrix()
print('curve tau (numerical integration):'); print(tau_c.apply_map(lambda z: z.n(digits=20)))
CF = ComplexField(90)   # our tau is accurate to about 0.24*p digits; do not ask for more
rels = integer_matrix_relations(tau_c.change_ring(CF), tau_s.change_ring(CF))
print('number of integer relations found:', len(rels))
Jstd = matrix(ZZ, 4, 4, [0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0])
for M in rels:
    print(M, ' det =', M.det(), ' symplectic:', M.transpose()*Jstd*M == Jstd or M.transpose()*Jstd*M == -Jstd)
