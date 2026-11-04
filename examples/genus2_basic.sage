# Basic usage of schottky_theta.py (run from the repository root:  sage examples/genus2_basic.sage)
#
# Genus 2 Schottky datum: discs of radius 1 centred at -2-+2i (D_{+-1}) and 2-+2i (D_{+-2}),
# generators sigma_i mapping the complement of D_{-i} onto the closure of D_i.
import os, sys
sys.path.insert(0, os.getcwd())
from schottky_theta import *

K = ComplexField(200)                 # working precision (bits)
i = K.gen()
discs = {1: Disc(-2-2*i, K(1)), -1: Disc(-2+2*i, K(1)),
         2: Disc( 2-2*i, K(1)), -2: Disc( 2+2*i, K(1))}

p = 60                                 # truncation order (number of Taylor coefficients)
S = SchottkyTheta(discs, p)            # precomputes the transfer matrices
F = S.solve()                          # coefficients of F_n^{(i)} for n = 1, 2 (one column per n)

z = K(3)/10 + K(1)/10*i
print('Omega_1(z) =', S.Omega(1, F.column(0), z).n(digits=30))
print('omega_1(z) =', S.omega(1, F.column(0), z).n(digits=30), '(coefficient of dz)')

B = S.period_matrix(F)                 # B_{mn} = Omega_n(sigma_m z) - Omega_n(z), defined modulo 2 pi i
print('period matrix B =')
print(B.apply_map(lambda x: x.n(digits=25)))
print('symmetry defect:', abs(B[0, 1] - B[1, 0]))
print('Riemann matrix tau = B/(2 pi i):')
print((B/(2*K.pi()*i)).apply_map(lambda x: x.n(digits=25)))

# Comparison with the truncated Poincare series (words of length <= 6).
value, terms = naive_Omega(S, 1, z, 6)
print('truncated series with %d terms differs by %.1e' % (terms, abs(value - S.Omega(1, F.column(0), z))))

# Your own generators: pass generators={1: matrix, 2: matrix}, each mapping the complement of
# discs[-i] onto the closure of discs[i]  (see examples/fro/fro.sage).
