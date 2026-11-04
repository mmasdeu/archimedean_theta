# Copyright (C) 2026 Marc Masdeu
#
# This file is part of archimedean_theta.
#
# archimedean_theta is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# archimedean_theta is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with archimedean_theta.  If not, see <https://www.gnu.org/licenses/>.

from itertools import chain
from sage.all import *
from darmonpoints.divisors import Divisors
from darmonpoints.meromorphic import divisor_to_pseries
from tqdm import tqdm
# Functions translated from Julia code in "Crossing the Transcendental Divide"
# Samantha Fairchild, Ángel David Ríos Ortiz

def matrix_from_balls(Bp, Bm):
    r'''
    Returns a matrix that takes Bm^c to the closure of Bp.
    '''
    a = Bp.center
    CC = a.parent()
    r = Bp.radius
    b = Bm.center
    s = Bm.radius
    ans = Matrix(CC,2,2,[1, -a, 0, r])**-1 * Matrix(CC,2,2,[0, 1, 1, 0]) * Matrix(CC, 2, 2, [1, -b, 0, s])
    ans.set_immutable()
    return ans

##Function for constructing centers and radii of the circles
def point_construct(eta, theta, g, prec):
    CC = ComplexField(prec)
    I = CC.gen()
    pi = 4 * CC(1).real().arctan()
    z1 = (I * ((theta+eta)/2)).exp() * (CC(theta - eta) / 2).real().sec() # Center of generating circle
    zs = [((j*2*pi*I) / g).exp() * z1 for j in srange(g)] # input first tuple of centers z1,....,zg
    ws = [((j*2*pi*I) / g).exp() * z1.conjugate() for j in srange(g)] #input second tuple of center w1,...., wg
    rs = [CC((theta-eta)/2).real().tan() for _ in srange(g)] #input radii of circles r1,...,rg
    balls = {i : Ball(w, r) for i, (w, r) in enumerate(zip(ws, rs), start=1)}
    balls.update({-i : Ball(z, r) for i, (z, r) in enumerate(zip(zs, rs), start=1)})
    return balls

##Function for Constructing Matrices of the Schottky Group
def matrix_construct(B1, B2):
    w = B1.center
    z = B2.center
    if B1.radius != B2.radius:
        raise NotImplementedError
    r = B1.radius
    CC = z.parent()
    I = CC.gen()
    argpt =  (w - z).arg() #Angle between isometric circle centers
    eiargpt = (-I*argpt).exp()
    A = matrix(CC,2,2,[(eiargpt*w)/r, -(w*z*eiargpt**2 + r**2)/(r* eiargpt), eiargpt/r, (-z* eiargpt)/r])
    return A
    
def enumerate_group_elements(gens, invgens, length):
    if length == 0:
        yield [], gens[0].parent().one()
    elif length == 1:
        for i, (v, vinv) in enumerate(zip(gens, invgens)):
            yield [i + 1], v
            yield [-i - 1], vinv
    else:
        for wd, g in enumerate_group_elements(gens, invgens, length - 1):
            for i, (v, vinv) in enumerate(zip(gens, invgens)):
                if wd[-1] != -(i + 1):
                    yield wd + [i + 1], g * v
                if wd[-1] != i + 1:
                    yield wd + [-i - 1], g * vinv


def all_elements_up_to_length(gens, invgens, N):
    return chain(*[enumerate_group_elements(gens, invgens, i) for i in srange(N + 1)])

def act(g, x):
    a, b, c, d = g.list()
    K = x.parent()
    if x == Infinity:
        return a / c
    if K(c) * x + K(d) == 0:
        return Infinity
    else:
        return (K(a) * x + K(b)) / (K(c) * x + K(d))

class Ball():
    def __init__(self, center, radius):
        self.center = center
        self.radius = radius
    def __contains__(self, z):
        return (z - self.center).abs() < self.radius

##Function for evaluating truncated series
def ptheta(samplepts, mat,g, L): #samples Poincare Theta series for all words of length < L, given as a function of z
    g = len(mat)
    CC = samplepts[0].parent()
    sampleproj = matrix(CC, len(samplepts), g)
    invmat = [~m for m in mat]
    list_of_words = list(all_elements_up_to_length(mat, invmat, L))
    for j, m in enumerate(mat, start=1): #create the jth coset G/ [mat[j]]
        jCoset_mat = []
        for wd, sig in list_of_words:
            if len(wd) == 0 or ZZ(wd[-1]).abs() != j:
                jCoset_mat.append(sig)
        ##Aj and Bj are the fixed points
        Aj, Bj = fixed_points(m)
        
        #Now for each point in sapmleproj, evaluate truncated series and place into sampleproj[k][j]
        for k, z in enumerate(samplepts):
            sampleproj[k, j-1] = sum(1 / (z - act(sig, Bj)) - 1 / (z - act(sig, Aj))
                for sig in jCoset_mat)
    return sampleproj

@cached_function
def get_action_data(Ps, g, oldparam, param, prec):
    verb_level = get_verbose()
    set_verbose(0)
    g.set_immutable()
    K = Ps.base_ring()
    tg = param.change_ring(K)
    a, b, c, d = (oldparam * (tg * g).adjugate()).list()
    zz = (Ps([b, a]) / Ps([d, c])).truncate(prec)
    ans = [zz.parent()(1), zz]
    while len(ans) < prec:
        zz_ps = (zz * ans[-1]).truncate(prec)
        ans.append(zz_ps)
    set_verbose(verb_level)
    m = Matrix(K, prec, prec, 0)
    for i, zz_ps in enumerate(ans):
        # zz_ps_log = zz_ps.log()
        for j, aij in enumerate(zz_ps):
            m[i, j] = aij  # i, j entry contains j-th term of zz^i
    return m

@cached_function
def find_parameter(Bp): # g takes Bp^c to unit disk
    a = Bp.center
    CC = a.parent()
    r = Bp.radius
    ans = Matrix(CC,2,2, [0,1,1,0]) * Matrix(CC,2,2,[1, -a, 0, r])
    ans.set_immutable()
    return ans

def fixed_points(g):
    a, b, c, d = g.list()
    jDisc = sqrt((a - d)**2 + 4*b*c) #discriminant is sqrt((a-d)^2 + 4bc)
    return ((a - d) + jDisc)/2/c, ((a - d) - jDisc)/2/c

def divisor_to_dict(D, balls):
    Div = D.parent()
    ans = {i : Div([]) for i in balls}
    for P, n in D:
        i = next((i for i, B in balls.items() if P in B))
        ans[i] += Div([(ZZ(n), P)])
    return ans

##Function for evaluating truncated series
def ptheta_oc(mat, idx, balls, prec, show_progress=True):
    idx = ZZ(idx)
    g = len(mat)
    CC = mat[0][0,0].parent()
    Ps = PowerSeriesRing(CC,names='t',default_prec=prec+1)
    Div = Divisors(CC)
    invmat = [~m for m in mat]
    G = [None] + mat + list(reversed(invmat))
    gens_ext = [(i, G[i]) for i in srange(1, g+1)] + [(-i, G[-i]) for i in srange(1,g+1)]
    param = {i : find_parameter(balls[i]) for i, _ in gens_ext}

    A, B = fixed_points(G[idx])
    D = Div([(ZZ(1),B),(-ZZ(1),A)]) # for (z-B)/(z-A)
    D1 = sum((gi * D for i, gi in gens_ext if ZZ(i).abs() != idx), Div([]))
    action_data = {
        (i, j) : get_action_data(Ps, gi, param[j], param[i], prec)
            for i, gi in gens_ext for j, _ in gens_ext if i != -j
        }
    F0 = {}
    for i, g in gens_ext:
        gD = sum((g * val for j, val in divisor_to_dict(D1, balls).items() if j != -i))
        F0[i] = divisor_to_pseries(param[i], Ps, gD, prec).log()
    Fnlist = [F0]
    progress = tqdm if show_progress else lambda x : x
    for _ in progress(range(prec)):
        tmp = {}
        for j, Fj in Fnlist[-1].items():        
            Fjexpvec = vector(CC, Fj.exp().list())
            for i, gi in gens_ext:
                if i != -j:
                    vl0 = Ps((Fjexpvec * action_data[i,j]).list())
                    vl = (vl0 / vl0[0]).log().add_bigoh(prec) # Fj must be a vector, action_data a matrix
                    try:
                        tmp[i] += vl
                    except KeyError:
                        tmp[i] = vl
        Fnlist.append(tmp)
    Fnlist = {ky : sum((F[ky] for F in Fnlist[1:]), Fnlist[0][ky]) for ky in Fnlist[0]}
    return D+D1, Fnlist, param

## Function for evaluating truncated series for Theta(z,A,B) = prod (z-gA) / (z-gB)
def theta_oc(mat, balls, A, B, prec, show_progress=True):
    g = len(mat)
    CC = mat[0][0,0].parent()
    Ps = PowerSeriesRing(CC,names='t',default_prec=prec+1)
    Div = Divisors(CC)
    invmat = [~m for m in mat]
    G = [None] + mat + list(reversed(invmat))
    gens_ext = [(i, G[i]) for i in srange(1, g+1)] + [(-i, G[-i]) for i in srange(1,g+1)]
    param = {i : find_parameter(balls[i]) for i, _ in gens_ext}

    D = Div([(ZZ(1),B),(-ZZ(1),A)]) # for (z-B)/(z-A)
    D1 = sum((gi * D for i, gi in gens_ext), Div([]))
    action_data = {
        (i, j) : get_action_data(Ps, gi, param[j], param[i], prec)
            for i, gi in gens_ext for j, _ in gens_ext if i != -j
        }
    F0 = {}
    for i, g in gens_ext:
        gD = sum((g * val for j, val in divisor_to_dict(D1, balls).items() if j != -i))
        F0[i] = divisor_to_pseries(param[i], Ps, gD, prec).log()
    Fnlist = [F0]
    progress = tqdm if show_progress else lambda x : x
    for _ in progress(range(prec)):
        tmp = {}
        for j, Fj in Fnlist[-1].items():        
            Fjexpvec = vector(CC, Fj.exp().list())
            for i, gi in gens_ext:
                if i != -j:
                    vl0 = Ps((Fjexpvec * action_data[i,j]).list())
                    vl = (vl0 / vl0[0]).log().add_bigoh(prec) # Fj must be a vector, action_data a matrix
                    try:
                        tmp[i] += vl
                    except KeyError:
                        tmp[i] = vl
        Fnlist.append(tmp)
    Fnlist = {ky : sum((F[ky] for F in Fnlist[1:]), Fnlist[0][ky]) for ky in Fnlist[0]}
    return D, D1, Fnlist, param

def eval(D, Phi, param, z):
    ans = sum(n/(z-P) for P, n in D)
    for ky, F in Phi.items():
        tt = act(param[ky], z)
        assert tt.abs() < 1.001, tt.abs()
        a, b, c, d = param[ky].list()
        ans += (a*d-b*c) / (c*z+d)**2 * F.derivative().polynomial()(tt)
    return ans

def eval_mul(D, Phi, param, z):
    ans = prod((z-P)**n for P, n in D)
    for ky, F in Phi.items():
        tt = act(param[ky], z)
        assert tt.abs() < 1.001, tt.abs()
        ans *= F.polynomial()(tt).exp()
    return ans


def SKprime(balls, z, w, prec):
    CC = z.parent()
    all_balls = {}
    mat = []
    for i, Bi in balls.items():
        all_balls[i] = Bi
        delta, q = Bi.center, Bi.radius
        den = delta.abs()**2 - q**2
        delta1 = delta / den
        q1 = q / den.abs()
        all_balls[-i] = Ball(delta1, q1)
        g = matrix(CC,2,2,[den, -delta, delta.conjugate(), -1])
        mat.append(g)
    _, D1, Fnlist, param = theta_oc(mat, all_balls, w, z, prec, show_progress=True)
    return (z-w)**2  * eval_mul(D1, Fnlist, param, w) / eval_mul(D1, Fnlist, param, z) 
