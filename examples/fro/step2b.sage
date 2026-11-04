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
S = SchottkyTheta(discs, 50, gens); F = S.solve(); tau = S.period_matrix(F)/(2*K.pi()*K.gen())
def cr(a,b,c,d): return (a-c)*(b-d)/((a-d)*(b-c))
for (m,n) in [(1,2),(1,3),(2,3)]:
    Am,Bm = S.fixed[m]; An,Bn = S.fixed[n]
    t = cr(Bm,Am,Bn,An).log()/(2*K.pi()*K.gen())
    print(m,n,'tau =', tau[m-1,n-1].n(digits=16), ' identity term =', t.n(digits=16), ' tau - identity term =', (tau[m-1,n-1]-t).n(digits=16))
# also the naive word-length truncation of the double coset sum without the identity, L<=8
def word_sum(m, n, L, include_identity):
    tot = K(0)
    Am,Bm = S.fixed[m]; An,Bn = S.fixed[n]
    def rec(g, wd):
        nonlocal tot
        if (wd and abs(wd[0]) != m and abs(wd[-1]) != n) or (not wd and include_identity):
            tot += cr(Bm, Am, mobius(g,Bn), mobius(g,An)).log()
        if len(wd) == L: return
        for j in S.I:
            if wd and wd[-1] == -j: continue
            rec(g*S.sigma[j], wd+[j])
    rec(S.sigma[1]**0, [])
    return tot/(2*K.pi()*K.gen())
for L in [6, 8]:
    print('L', L, 'without identity', word_sum(1,2,L,False).n(digits=16), ' with identity', word_sum(1,2,L,True).n(digits=16))
