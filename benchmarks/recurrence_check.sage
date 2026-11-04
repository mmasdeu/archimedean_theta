import os
# Stability check of the three-term recurrence of Lemma 4.1: compare the transfer matrices with those obtained by
# repeated multiplication of power series, at the working precision and at p = 60 and 120.
import sys
sys.path.insert(0, os.getcwd())
from schottky_theta import *
load('examples/fro/fro.sage')
def by_series(p, sigma, ti, tj, K):
    a, b, c, d = (tj*sigma.adjugate()*ti.adjugate()).list()
    R = PowerSeriesRing(K, 'u', default_prec=p); u = R.gen()
    phi = ((a*u + b)/(c*u + d)).add_bigoh(p); M = Matrix(K, p, p); pw = R(1).add_bigoh(p)
    for k in range(p):
        co = pw.padded_list(p)
        for j in range(p): M[k, j] = co[j]
        pw = (pw*phi).add_bigoh(p)
    return M
def check(name, discs, gens, p, bits):
    K = discs[1].center.parent()
    S = SchottkyTheta(discs, 2, gens)
    worst = 0; ac_max = 0
    for i in S.I:
        for j in S.I:
            if j == -i: continue
            a, b, c, d = (S.t[j]*S.sigma[i].adjugate()*S.t[i].adjugate()).list()
            ac_max = max(ac_max, float(abs(a*c)/abs(b*d)))
            M1 = substitution_matrix(p, S.sigma[i], S.t[i], S.t[j]); M2 = by_series(p, S.sigma[i], S.t[i], S.t[j], K)
            scale = max(1, float(max(abs(x) for x in M2.list())))
            worst = max(worst, float(max(abs(x - y) for x, y in zip(M1.list(), M2.list()))) / scale)
    print('%-28s p=%3d bits=%d  max |ac/bd| = %.2f   max relative difference recurrence vs series = %.1e  (2^-bits = %.1e)' % (name, p, bits, ac_max, worst, 2.0**-bits), flush=True)
for p, bits in [(60, 120), (120, 200)]:
    K = ComplexField(bits); i_ = K(0, 1); r = K(1)
    D = {1:Disc(-2-2*i_,r),-1:Disc(-2+2*i_,r),2:Disc(2-2*i_,r),-2:Disc(2+2*i_,r)}
    check('genus 2 datum (Section 7)', D, None, p, bits)
    discs, gens = fro_datum(K.pi()/12, K.pi()/4, 2, K)
    check('Fairchild--Rios Ortiz, g=2', discs, gens, p, bits)
