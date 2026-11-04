# The first implementation (theta_arch.py, needs the darmonpoints package), genus 2.
# Run from the repository root:  sage examples/theta_arch_genus2.sage
import os, sys, time
sys.path.insert(0, os.getcwd())
from theta_arch import *

prec = 200
CC = ComplexField(prec)
i = CC.gen()
balls = {1: Ball(-2-2*i, 1), -1: Ball(-2+2*i, 1), 2: Ball(2-2*i, 1), -2: Ball(2+2*i, 1)}
mat = [matrix_from_balls(balls[j], balls[-j]) for j in [1, 2]]
samplepts = [CC(0), CC(.5), .5*i]

t0 = time.time()
A = ptheta(samplepts, mat, 2, 8)       # truncated Poincare series (words of length <= 8)
print('truncated series: %.2f s' % (time.time() - t0))
for k in [10, 20, 40]:
    t0 = time.time()
    D1, Phi1, param1 = ptheta_oc(mat, 1, balls, k, show_progress=False)
    D2, Phi2, param2 = ptheta_oc(mat, 2, balls, k, show_progress=False)
    B = Matrix(CC, len(samplepts), 2, 0)
    for j, s in enumerate(samplepts):
        B[j, 0] = eval(D1, Phi1, param1, s)
        B[j, 1] = eval(D2, Phi2, param2, s)
    print('p=%d: %.2f s, difference with the truncated series %.1e' % (k, time.time() - t0, (B - A).norm()))
