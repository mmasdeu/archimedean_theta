load('examples/x0_23/model.sage')
y120 = load('examples/x0_23/y_p120.sobj'); y200 = load('examples/x0_23/y_p200.sobj')
R = y200[0].parent()
print('precision bits', R.prec())
for name, a, b in zip(['r1', 'c2', 'r2'], y120, y200):
    print('%s = %s   (|p120 - p200| = %.1e)' % (name, b.str(digits=48), float(abs(R(a)-b))))
# sensitivity: Jacobian at p=60
def F(y, p, bits):
    RR_ = RealField(bits); J, _ = invariants([1, y[1]], [y[0], y[2]], p, bits)
    return vector(RR_, [RR_(v.real()) for v in J])
b0 = 120; h = RealField(b0)(10)**-10; yy = vector(RealField(b0), y200)
Jac = matrix(RealField(b0), 3, 3)
for k in range(3):
    e = vector(RealField(b0), [0,0,0]); e[k] = h
    Jac.set_column(k, (F(yy+e, 60, b0) - F(yy-e, 60, b0))/(2*h))
print('||J^{-1}|| (max row sum) = %.2e' % float(Jac.inverse().apply_map(abs).norm(Infinity)))
# Schottky data
K = ComplexField(R.prec()); discs, gens = datum([1, y200[1]], [y200[0], y200[2]], K)
for j in (1, 2):
    A, B = fixed_points(gens[j]); a, b, c, d = gens[j].list()
    mu = (a*d - b*c)/(c*B + d)**2
    print('S_%d = %s' % (j, [v.real().n(digits=30) for v in gens[j].list()]))
    print('   A_%d = %s, B_%d = %s, mu_%d = %s' % (j, A.n(digits=30), j, B.n(digits=30), j, mu.n(digits=30)))
