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

import os
# Numbers of the "Correctness checks" subsection (genus 2 datum, discs of radius 1 at +-2+-2i).
import sys
sys.path.insert(0, os.getcwd())
from schottky_theta import *
def config(K):
    i = K(0,1); r = K(1)
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
for bits, p in [(150, 30), (150, 60), (200, 70)]:
    K = ComplexField(bits); S = SchottkyTheta(config(K), p); F = S.solve()
    B = S.period_matrix(F); B2 = S.period_matrix(F, angle=2.1)
    print('bits=%d p=%d  symmetry defect %.1e  base point dependence %.1e' % (bits, p, float(abs(B[0,1]-B[1,0])), float(max(abs(B[a,b]-B2[a,b]) for a in range(2) for b in range(2)))))
    print('   B11 =', B[0,0].real().n(digits=30), ' B12 =', B[0,1].real().n(digits=30))
z = K(3)/10 + K(0,1)/10
for L in [6, 8]:
    v, c = naive_Omega(S, 1, z, L)
    print('words of length <= %d: %d terms, |difference| = %.1e' % (L, c, float(abs(v - S.Omega(1, F.column(0), z)))))
