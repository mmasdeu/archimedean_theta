load('examples/x0_23/model.sage')
for p in [16, 30]:
    t0 = time.time(); J, vals = invariants([1, 3], [0.3, 1.0], p, 100); print('p', p, 'time %.2f' % (time.time()-t0), [v.n(digits=10) for v in J])
print('weierstrass values', [v.n(digits=8) for v in vals])
print('target', TARGET.n(digits=10))
