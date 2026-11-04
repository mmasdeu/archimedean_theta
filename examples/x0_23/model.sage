import os
# Bogatyrev's real Schottky model of genus 2 with two circles of type Sigma^- (one real oval):
# circles C_j = C(c_j, r_j) on the positive half-line, C_{-j} = -C_j, generators
# S_j(u) = G_j(-u) = c_j + r_j^2/(u + c_j), G_j(u) = c_j - r_j^2/(u - c_j) fixing c_j +- i r_j.
# The hyperelliptic involution is u -> -u; the Weierstrass points are 0, oo, c_j +- i r_j.
import sys, time
sys.path.insert(0, os.getcwd())
from schottky_theta import *
from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants

def datum(c, r, K):
    discs, gens = {}, {}
    for j in (1, 2):
        cj, rj = K(c[j-1]), K(r[j-1])
        discs[j] = Disc(cj, rj); discs[-j] = Disc(-cj, rj)
        gens[j] = Matrix(K, 2, 2, [cj, cj**2 + rj**2, 1, cj])
    return discs, gens

def omega_at_infinity(S, n, Fn):
    # coefficient of 1/z^2 in omega_n(z) as z -> oo
    A, B = S.fixed[n]
    return (B - A) - sum(S.discs[i].radius * Fn[S.idx[i]*S.q] for i in S.I)

def weierstrass_values(c, r, p, bits):
    K = ComplexField(bits); I = K.gen()
    discs, gens = datum(c, r, K)
    S = SchottkyTheta(discs, p, gens)
    F = S.solve()
    F1, F2 = F.column(0), F.column(1)
    vals = [omega_at_infinity(S, 2, F2)/omega_at_infinity(S, 1, F1)]
    for z in [K(0)] + [K(c[j]) + s*I*K(r[j]) for j in (0, 1) for s in (1, -1)]:
        vals.append(S.omega(2, F2, z)/S.omega(1, F1, z))
    return vals, S, F

def invariants(c, r, p, bits):
    vals, S, F = weierstrass_values(c, r, p, bits)
    K = vals[0].parent()
    x = PolynomialRing(K, 'x').gen()
    I2, I4, I6, I10 = igusa_clebsch_invariants(prod(x - v for v in vals))
    return vector(K, [I4/I2**2, I2*I4/I6, I10/I2**5]), vals

TARGET = vector(QQ, [669124/(-1144)**2, (-1144)*669124/(-17115576), 606355001344/(-1144)**5])
