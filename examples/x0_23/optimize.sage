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
from scipy.optimize import minimize
import numpy as np
T = [float(t) for t in TARGET]
def params_from(x):
    a = 1/(1+np.exp(-x[0]))          # r1 in (0,1)
    s, r2 = np.exp(x[1]), np.exp(x[2])
    return [1, 1 + a + s + r2], [a, r2]
def misfit(x, p=14, bits=80):
    c, r = params_from(x)
    try:
        J, _ = invariants(c, r, p, bits)
    except Exception:
        return 1e6
    return float(sum((np.arcsinh(float(v.real())/abs(t)) - np.arcsinh(t/abs(t)))**2 for v, t in zip(J, T)))
starts = [(0.1,0.1,0.1), (0.1,0.3,0.1), (0.1,0.1,0.3), (0.2,0.3,0.3), (0.65,0.1,1.5), (0.35,0.3,0.7)]
results = []
for a, s, r2 in starts:
    x0 = np.array([np.log(a/(1-a)), np.log(s), np.log(r2)])
    res = minimize(misfit, x0, method='Nelder-Mead', options={'xatol':1e-6, 'fatol':1e-12, 'maxiter':400})
    c, r = params_from(res.x)
    print('start', (a,s,r2), '-> misfit %.3e' % res.fun, ' c =', [float(v) for v in c], ' r =', [float(v) for v in r], flush=True)
    results.append((res.fun, res.x))
results.sort(key=lambda t: t[0])
save(results[0][1].tolist(), 'examples/x0_23/best_x.sobj')
