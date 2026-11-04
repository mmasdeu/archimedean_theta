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
Kr=ComplexBallField(500); z0=Kr(3)/10+Kr(1)/10*Kr(0,1)
t0=time.time(); S=SchottkyTheta(config(Kr),150); F=S.solve(); v=S.Omega(1,F.column(0),z0)
re=v.real().mid().exact_rational(); im=v.imag().mid().exact_rational()
open('benchmarks/server/ref.tmp','w').write('%s\n%s\n'%(re,im)); os.rename('benchmarks/server/ref.tmp','benchmarks/server/ref.txt')
print('ref written, time',time.time()-t0,'radius',float(v.real().rad()),flush=True)
