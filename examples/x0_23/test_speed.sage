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
for p in [16, 30]:
    t0 = time.time(); J, vals = invariants([1, 3], [0.3, 1.0], p, 100); print('p', p, 'time %.2f' % (time.time()-t0), [v.n(digits=10) for v in J])
print('weierstrass values', [v.n(digits=8) for v in vals])
print('target', TARGET.n(digits=10))
