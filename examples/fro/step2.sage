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

load('examples/fro/fro.sage')
K = ComplexField(250)
eta, theta = K.pi()/12, K.pi()/4
discs, gens = fro_datum(eta, theta, 3, K)
res = {}
for p in [30, 50]:
    t0 = time.time(); S = SchottkyTheta(discs, p, gens); F = S.solve(); P = S.period_matrix(F); tt = time.time()-t0
    res[p] = P/(2*K.pi()*K.gen())
    print('p', p, 'time %.1f'%tt, 'symm defect %.1e'%max(float(abs(P[a,b]-P[b,a])) for a in range(3) for b in range(3)), flush=True)
print('max |tau30 - tau50| = %.1e' % max(float(abs(res[30][a,b]-res[50][a,b])) for a in range(3) for b in range(3)))
print(res[50].apply_map(lambda x: x.n(digits=25)))
fro = matrix(CC, 3, 3, [0.39980895499614727, -0.004497125148067202, -0.004497125148229181, -0.004497125148067202, 0.3998089549956627, -0.004497125149408407, -0.004497125148229181, -0.004497125149408407, 0.3998089549963685])*CC.gen()
print('max |tau - FRO| = %.1e' % max(abs(CC(res[50][a,b]) - fro[a,b]) for a in range(3) for b in range(3)))
