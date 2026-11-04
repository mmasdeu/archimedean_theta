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
def config(K,r):
    i=K(0,1); r=K(r)
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
K=ComplexField(200); z=K(3)/10+K(0,1)/10
for r in [K(16)/10, K(18)/10]:
    S=SchottkyTheta(config(K,r),100); F=Matrix(K,[S.iterate(n,100) for n in (1,2)]).transpose(); ref=S.Omega(1,F.column(0),z)
    out=[]
    for p in [10,20,40,60]:
        Sp=SchottkyTheta(config(K,r),p); Fp=Matrix(K,[Sp.iterate(n,p) for n in (1,2)]).transpose(); out.append((p,'%.1e'%float(abs(Sp.Omega(1,Fp.column(0),z)-ref))))
    print(float(r),out,flush=True)
