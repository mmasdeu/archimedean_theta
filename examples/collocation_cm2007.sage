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

# The collocation algorithm of Crowdy--Marshall (2007), reimplemented in collocation.py, on their example:
# X = omega^2 at z = (-1-i)/2, w = 1 for |z| < 1, |z - 1/2| > 1/10, |z - i/2| > 1/10.
# Run from the repository root:  sage examples/collocation_cm2007.sage     (about 30 s)
import os, sys, time
sys.path.insert(0, os.getcwd())
from collocation import Domain, prime_square
K = ComplexBallField(200); i = K(0, 1)
D = Domain([K(1)/2, K(0, 1)/2], [K(1)/10, K(1)/10])
z = (-1 - i)/2; w = K(1)
published = ComplexField(200)('2.39754812221763') + ComplexField(200)('1.76164377547306')*ComplexField(200).gen()
for N1 in [10, 20, 30]:
    t0 = time.time()
    X = ComplexField(200)(prime_square(D, N1, 200, [z], w)[0].mid())
    print('N1=%d (%.1f s): X = %s   |X - published| = %.1e' % (N1, time.time() - t0, X.n(digits=20), abs(X - published)))
