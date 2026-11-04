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
c, r = [1, 1.2883496775840126], [0.16675309408461045, 0.07913698863104919]
ref = None; out = []
for p in [20, 40, 60, 80, 100]:
    t0 = time.time(); J, _ = invariants(c, r, p, 300); tt = time.time() - t0
    out.append((p, tt, J))
ref = out[-1][2]
for p, tt, J in out:
    print('p=%d time=%.1f  |J - J(p=100)| = %.1e' % (p, tt, float(max(abs(a-b) for a, b in zip(J, ref)))), flush=True)
print('J(p=100) =', [v.n(digits=20) for v in ref]); print('target   =', TARGET.n(digits=20))
