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

import sys,time; sys.path.insert(0, os.getcwd())
from schottky_theta import *
import math
K=ComplexField(160)
def ring(g,K):
    # 2g discs: centres on the unit circle? use the Fairchild--Rios Ortiz construction
    eta=K(pi)/(6*g); theta=K(pi)/(2*g)
    I=K(0,1)
    z1=(I*((theta+eta)/2)).exp()*(((theta-eta)/2).cos())**-1
    r=((theta-eta)/2).tan()
    zs=[(2*K(pi)*I*j/g).exp()*z1 for j in range(g)]
    ws=[(2*K(pi)*I*j/g).exp()*z1.conjugate() for j in range(g)]
    D={}
    for j in range(g):
        D[j+1]=Disc(ws[j],r); D[-(j+1)]=Disc(zs[j],r)
    return D
for g in [2,3,4,5,6]:
    D=ring(g,K)
    ks=list(D); gap=min(float(abs(D[a].center-D[b].center)-D[a].radius-D[b].radius) for a in ks for b in ks if a<b)
    res=[]
    Fs={}
    for p in [20,40]:
        t0=time.time(); S=SchottkyTheta(D,p); tb=time.time()-t0
        t0=time.time(); F=Matrix(K,[S.iterate(n,p) for n in range(1,g+1)]).transpose(); ts=time.time()-t0
        P=S.period_matrix(F); Fs[p]=P
        asym=max(float(abs(P[a,b]-P[b,a])) for a in range(g) for b in range(g))
        res.append((p,round(tb,1),round(ts,1),'%.1e'%asym))
    diff=max(float(abs(Fs[20][a,b]-Fs[40][a,b])) for a in range(g) for b in range(g))
    print('g',g,'gap %.3f'%gap,res,'|P20-P40| %.1e'%diff,flush=True)
