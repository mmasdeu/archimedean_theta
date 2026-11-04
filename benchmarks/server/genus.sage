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

load('benchmarks/server/common.sage')
K=ComplexField(160)
def ring(g,K):
    eta=K(pi)/(6*g); theta=K(pi)/(2*g); I=K(0,1)
    z1=(I*((theta+eta)/2)).exp()*(((theta-eta)/2).cos())**-1
    r=((theta-eta)/2).tan()
    D={}
    for j in range(g):
        D[j+1]=Disc((2*K(pi)*I*j/g).exp()*z1.conjugate(),r); D[-(j+1)]=Disc((2*K(pi)*I*j/g).exp()*z1,r)
    return D
def period(S,F,m,n):
    pass
for g in [2,3,4,5,6]:
    D=ring(g,K); row=[]
    for p in [20,40]:
        t0=time.time(); S=SchottkyTheta(D,p); F=Matrix(K,[S.iterate(n,p) for n in range(1,g+1)]).transpose(); tt=time.time()-t0
        P=S.period_matrix(F); asym=max(float(abs(P[a,b]-P[b,a])) for a in range(g) for b in range(g))
        row.append((p,round(tt,1),'%.1e'%asym))
    print('g',g,row,flush=True)
