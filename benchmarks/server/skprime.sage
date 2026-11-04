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

import sys,time,os
sys.path.insert(0, os.getcwd())   # run from the repository root
from theta_arch import *
CC=ComplexField(400); i=CC.gen()
balls={1:Ball(CC(1)/2,CC(1)/10),2:Ball(CC.gen()/2,CC(1)/10)}
z=(-1-i)/2; w=CC(1)
for prec in [10,20,40,60]:
    t0=time.time(); X=SKprime(balls,z,w,prec); print('SKprime p=%d time=%.2f X=%s'%(prec,time.time()-t0,X.n(digits=30)),flush=True)
