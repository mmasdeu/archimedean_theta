# Check that y^2 = (x^3-x+1)(x^3-8x^2+3x-7) is a model of X_0(23): compare #C(F_p) with the newforms of level 23.
from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
R.<x> = QQ[]
f = (x^3 - x + 1)*(x^3 - 8*x^2 + 3*x - 7)
print('f =', f, '  disc =', factor(f.discriminant()))
nf = CuspForms(Gamma0(23), 2).newforms('a')
print('newforms:', nf)
ok = True
for p in prime_range(3, 60):
    if p == 23: continue
    C = HyperellipticCurve(f.change_ring(GF(p)))
    Np = C.count_points(1)[0]
    tr = sum(F[p] for F in nf) if len(nf) == 2 else nf[0][p].trace()
    if Np != p + 1 - tr: ok = False; print('mismatch at', p, Np, p+1-tr)
print('point counts match the level-23 newforms for 3 <= p < 60:', ok)
IC = igusa_clebsch_invariants(f)
print('Igusa-Clebsch invariants:', IC)
I2, I4, I6, I10 = IC
print('absolute (I2^5/I10, I2^3 I4/I10, I2^2 I6/I10):', [I2^5/I10, I2^3*I4/I10, I2^2*I6/I10])
print('real roots of f:', [r for r in f.roots(RR)], ' leading coeff', f.leading_coefficient())
