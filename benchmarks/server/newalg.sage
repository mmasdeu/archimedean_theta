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

# Table 2, new algorithm: smallest p with error < 10^-d; time = min over 3 runs of build+solve (all n)
load('benchmarks/server/common.sage')
Kr=ComplexBallField(500)
load('benchmarks/server/refload.sage')
ref=load_ref(Kr)
for d in [5,10,11,12,15,20,40,80]:
    K=ComplexBallField(int(d*3.33)+40); z=K(3)/10+K(0,1)/10
    for p in range(4,200):
        S=SchottkyTheta(config(K),p); F=S.solve(); err=abs(Kr(S.Omega(1,F.column(0),z))-ref)
        if err<10^(-d): break
    ts=[]
    for _ in range(3):
        t0=time.time(); S=SchottkyTheta(config(K),p); F=S.solve(); ts.append(time.time()-t0)
    print('new d=%d p=%d time=%.3f err=%.1e'%(d,p,min(ts),float(err)),flush=True)
