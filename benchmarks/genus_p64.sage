import os
# Accuracy of p=40 for the ring family: comparison with p=64 (genus 2 and 3), 200 bits.
import sys
sys.path.insert(0, os.getcwd())
from schottky_theta import *
K = ComplexField(200)
def ring(g, K):
    eta = K.pi()/(6*g); theta = K.pi()/(2*g); I = K(0,1)
    z1 = (I*((theta+eta)/2)).exp()/((theta-eta)/2).cos(); r = ((theta-eta)/2).tan()
    D = {}
    for j in range(g):
        D[j+1] = Disc((2*K.pi()*I*j/g).exp()*z1.conjugate(), r); D[-(j+1)] = Disc((2*K.pi()*I*j/g).exp()*z1, r)
    return D
for g in [2, 3]:
    P = {}
    for p in [40, 64]:
        S = SchottkyTheta(ring(g, K), p); F = Matrix(K, [S.iterate(n, p) for n in range(1, g+1)]).transpose(); P[p] = S.period_matrix(F)
    print('g=%d  max |B(p=40) - B(p=64)| = %.1e' % (g, float(max(abs(P[40][a,b]-P[64][a,b]) for a in range(g) for b in range(g)))), flush=True)
