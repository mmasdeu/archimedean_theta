import os
import sys, time
sys.path.insert(0, os.getcwd())
from schottky_theta import *

def fro_datum(eta, theta, g, K):
    # point_construct and matrix_construct of theta_arch.py (translated from the Julia code of Fairchild--Rios Ortiz)
    I = K.gen(); pi_ = K.pi()
    z1 = (I*((theta+eta)/2)).exp() / ((theta-eta)/2).cos()
    r = ((theta-eta)/2).tan()
    zs = [((2*pi_*I*j)/g).exp()*z1 for j in range(g)]
    ws = [((2*pi_*I*j)/g).exp()*z1.conjugate() for j in range(g)]
    discs = {}
    for j in range(g):
        discs[j+1] = Disc(ws[j], r); discs[-(j+1)] = Disc(zs[j], r)
    gens = {}
    for j in range(1, g+1):
        w = discs[j].center; z = discs[-j].center
        arg = (w - z).arg(); e = (-I*arg).exp()
        A = Matrix(K, 2, 2, [(e*w)/r, -(w*z*e**2 + r**2)/(r*e), e/r, (-z*e)/r])
        gens[j] = A
    return discs, gens

def check_orientation(discs, gens):
    for j, A in gens.items():
        a, b, c, d = A.list()
        img_inf = a/c
        inDj = abs(img_inf - discs[j].center) < discs[j].radius
        inDmj = abs(img_inf - discs[-j].center) < discs[-j].radius
        print('gen', j, 'sigma(inf) in D_j:', inDj, ' in D_-j:', inDmj, ' trace:', (a+d)/sqrt(a*d-b*c))
