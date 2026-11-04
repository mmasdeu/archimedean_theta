# Genus 2 member of the Fairchild--Rios Ortiz family with the first implementation (needs darmonpoints).
# Run from the repository root:  sage examples/fro_family_theta_arch.sage
import os, sys, time
sys.path.insert(0, os.getcwd())
from theta_arch import *

prec = 100
CC = ComplexField(prec)
I = CC.gen()
eta, theta = CC(pi)/12, CC(pi)/4
g = 2
balls = point_construct(eta, theta, g, prec)
mat = [matrix_construct(balls[j], balls[-j]) for j in [1, 2]]
samplepts = [CC(1), CC(-1), (theta*I).exp(), (CC(pi)*I - theta*I).exp(), (eta*I).exp(), (CC(pi)*I - eta*I).exp()]

t0 = time.time(); A = ptheta(samplepts, mat, g, 8); print('truncated series: %.1f s' % (time.time() - t0))
B0 = Matrix(CC, len(samplepts), 2, 0)
for k in [10, 20, 30]:
    D1, Phi1, param1 = ptheta_oc(mat, 1, balls, k, show_progress=False)
    D2, Phi2, param2 = ptheta_oc(mat, 2, balls, k, show_progress=False)
    B = Matrix(CC, len(samplepts), 2, 0)
    for j, s in enumerate(samplepts):
        B[j, 0] = eval(D1, Phi1, param1, s)
        B[j, 1] = eval(D2, Phi2, param2, s)
    print('p=%d: change %.1e, difference with truncated series %.1e' % (k, (B - B0).norm(), (B - A).norm()))
    B0 = B
