load('examples/x0_23/model.sage')
# parametrization: c1 = 1, r1 = a, c2 = (1 + a) + s + r2, r2 = b*(...) ; variables a in (0,1), s > 0 gap, r2 > 0
def params(a, s, r2):
    return [1, 1 + a + s + r2], [a, r2]
T = TARGET.change_ring(RR)
def score(J):
    J = [RR(v.real()) for v in J]
    sc = 0
    for v, t in zip(J, T):
        if v == 0 or sign(v) != sign(t): sc += 10
        else: sc += abs(log(abs(v)) - log(abs(t)))
    return sc
best = []
t0 = time.time()
for a in [0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9]:
    for s in [0.1, 0.3, 0.7, 1.5, 3, 6]:
        for r2 in [0.1, 0.3, 0.7, 1.5, 3, 6]:
            c, r = params(a, s, r2)
            try:
                J, _ = invariants(c, r, 14, 80)
                best.append((score(J), a, s, r2))
            except Exception as e:
                pass
best.sort()
print('time %.0f s' % (time.time()-t0))
for b in best[:10]: print(b)
