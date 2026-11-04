# Square X = omega(z, w)^2 of the Schottky-Klein prime function of the triply connected domain
# |z| < 1, |z - 1/2| > 1/10, |z - i/2| > 1/10, at z = (-1-i)/2, w = 1  (Crowdy--Marshall 2007, Section 5).
# Needs darmonpoints.  Run from the repository root:  sage examples/crowdy_marshall_prime.sage
import os, sys, time
sys.path.insert(0, os.getcwd())
from theta_arch import *

CC = ComplexField(400)
i = CC.gen()
balls = {1: Ball(CC(1)/2, CC(1)/10), 2: Ball(CC.gen()/2, CC(1)/10)}
z = (-1 - i)/2
w = CC(1)
published = CC('2.39754812221763') + CC('1.76164377547306')*i
for p in [10, 20, 40]:
    t0 = time.time()
    X = SKprime(balls, z, w, p)
    print('p=%d (%.1f s): X = %s   |X - published| = %.1e' % (p, time.time() - t0, X.n(digits=30), abs(X - published)))
