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

# Algebraicity up to powers of pi of the absolute Igusa--Clebsch invariants of S_{pi/12, pi/4},
# compared with random controls of the same precision (identical PARI calls).
from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
phis = load('examples/fro/phis_240.sobj'); R0 = phis[0].parent()
x = PolynomialRing(R0, 'x').gen(); f = prod(x - ph for ph in phis)
I2, I4, I6, I10 = igusa_clebsch_invariants(f)
D = 100
R = RealField(int(D*3.33)); PI = R.pi()
pari.set_real_precision(D)
set_random_seed(1)
def norm(v):
    return v / 10**floor(log(abs(v), 10))
def lh(rel): return float(log(max(abs(ZZ(c)) for c in rel), 10))
def tests(w):
    out = {}
    for d in range(1, 5):   # (a) w/pi^n algebraic of degree d, |n| <= 30; min over n
        out['a%d' % d] = min(lh(pari(norm(w / PI**n)).algdep(d)) for n in range(-30, 31))
    for M in [2, 4, 6, 8]:  # (b) w in Q[pi]
        out['b%d' % M] = lh(pari([norm(w)] + [PI**j for j in range(M+1)]).lindep())
    for M in [1, 2, 3, 4]:  # (c) w in Q(pi)
        out['c%d' % M] = lh(pari([norm(w)*PI**j for j in range(M+1)] + [PI**j for j in range(M+1)]).lindep())
    return out
names = ['i1', 'i2', 'i3']
res = {nm: tests(R(v)) for nm, v in zip(names, [I2**5/I10, I2**3*I4/I10, I2**2*I6/I10])}
ctrl = [tests(R.random_element(1, 10)) for _ in range(40)]
keys = list(res['i1'].keys())
print('test   ' + '  '.join('%6s' % n for n in names) + '   controls: min / median / max  (log10 of height)')
for k in keys:
    cs = sorted(c[k] for c in ctrl)
    print('%-5s  ' % k + '  '.join('%6.1f' % res[n][k] for n in names) + '   %6.1f / %6.1f / %6.1f' % (cs[0], cs[len(cs)//2], cs[-1]))
