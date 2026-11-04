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
def config(K,kind):
    i=K(0,1)
    r=K(1) if kind=='g2' else K(16)/10
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
for kind in ['g2','g2close']:
    print('=====',kind,flush=True)
    Kr=ComplexBallField(500)
    z0=Kr(3)/10+Kr(1)/10*Kr(0,1)
    Sref=SchottkyTheta(config(Kr,kind),150); Fr=Sref.solve()
    ref=Sref.Omega(1,Fr.column(0),z0); print('ref done',flush=True)
    # decay rate: how many digits per coefficient
    for d in [5,10,15,20,40,80]:
        K=ComplexBallField(int(d*3.33)+40); z=K(3)/10+K(0,1)/10
        for p in range(4,152,2):
            t0=time.time(); S=SchottkyTheta(config(K,kind),p); F=S.solve(); tt=time.time()-t0
            err=abs(Kr(S.Omega(1,F.column(0),z)) - ref)
            if err<10^(-d): print('ours  d=%d p=%d time(all n, incl. precomputation)=%.2fs err=%.1e'%(d,p,tt,float(err)),flush=True); break
        if d<=10:
            K=ComplexField(int(d*3.33)+40); z=K(3)/10+K(0,1)/10
            S=SchottkyTheta(config(K,kind),3)
            for mu in [10^(-d+1)*10^(-k) for k in range(0,9)]:
                t0=time.time(); a,c=adaptive_Omega(S,1,z,K(mu)); tt=time.time()-t0
                err=abs(Kr(a)-ref)
                if err<10^(-d): print('  adaptive d=%d mu=%.0e terms=%d time(one n)=%.2fs err=%.1e'%(d,mu,c,tt,float(err)),flush=True); break
            else: print('  adaptive failed d=%d'%d)
