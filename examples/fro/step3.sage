load('examples/fro/fro.sage')
from sage.schemes.hyperelliptic_curves.invariants import igusa_clebsch_invariants
import sys
eta_q, theta_q = QQ(sys.argv[-3]), QQ(sys.argv[-2]); p = int(sys.argv[-1])   # angles as multiples of pi
bits = int(p*0.6*3.33) + 64
K = ComplexField(bits)
eta, theta = eta_q*K.pi(), theta_q*K.pi()
discs, gens = fro_datum(eta, theta, 2, K)
t0 = time.time()
S = SchottkyTheta(discs, p, gens)
F = Matrix(K, [S.iterate(n, p) for n in (1, 2)]).transpose()
t1 = time.time() - t0
I = K.gen()
xs = [K(1), K(-1), (I*theta).exp(), (I*eta).exp(), (I*(K.pi()-theta)).exp(), (I*(K.pi()-eta)).exp()]
phis = [S.omega(2, F.column(1), x)/S.omega(1, F.column(0), x) for x in xs]
R = RealField(bits)
print('p=%d bits=%d time=%.1f' % (p, bits, t1))
print('max imaginary part of phi_j: %.1e' % max(float(abs(f.imag())) for f in phis))
print('phi_j:', [ph.real().n(digits=15) for ph in phis])
print('pair products:', [(a+1,b+1,(phis[a]*phis[b]).real().n(digits=12)) for a in range(6) for b in range(a+1,6)])
save([ph.real() for ph in phis], 'examples/fro/phis_%d.sobj' % p)
Rx = PolynomialRing(R, 'x'); x = Rx.gen()
f = prod(x - R(ph.real()) for ph in phis)
I2, I4, I6, I10 = igusa_clebsch_invariants(f)
i1, i2, i3 = I2**5/I10, I2**3*I4/I10, I2**2*I6/I10
for name, v in [('i1',i1),('i2',i2),('i3',i3)]:
    print(name, v.str(digits=int(bits/3.33)-15, no_sci=False))
