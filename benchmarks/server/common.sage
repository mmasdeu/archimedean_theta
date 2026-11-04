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
from schottky_theta import *
def config(K,r=1):
    i=K(0,1); r=K(r)
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
