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

r'''
Multiprecision implementation (Sage) of the collocation algorithm of Crowdy and Marshall,
"Computing the Schottky-Klein prime function on the Schottky double of planar domains",
Comput. Methods Funct. Theory 7 (2007), 293-308, Section 4, written to compare its cost with the
algorithm of this repository.  It computes X(z, w) = omega(z, w)^2 for the domain
|z| < 1, |z - delta_k| > q_k (k = 1..M).

Differences with the reference (which is a Matlab script working in double precision):
  * everything is done at a prescribed precision, the least-squares problems being solved through the
    normal equations at about twice the working precision (the loss of digits being
    controlled by comparing different truncation orders);
  * the linear algebra is done with ball matrices (Flint), the midpoints being used.
'''
from sage.all import (QQ, ComplexBallField, RealBallField, ComplexField, Matrix, vector, pi, exp, log)


class Domain:
    def __init__(self, centers, radii):
        self.delta = list(centers)           # centres delta_k of the circles C_k
        self.q = list(radii)                 # radii q_k
        self.M = len(self.delta)
        # circles C_k' : reflections of C_k in the unit circle
        self.deltap = []
        self.qp = []
        for d, q in zip(self.delta, self.q):
            den = abs(d) ** 2 - q ** 2
            self.deltap.append(d / den)
            self.qp.append(q / abs(den))
        self.qtilde = [qp if float(abs(qp).mid()) < 1 else qp.parent()(1) for qp in self.qp]

    def theta(self, k, z):
        d, q = self.delta[k], self.q[k]
        return d + q ** 2 * z / (1 - d.conjugate() * z)

    def dtheta(self, k, z):
        d, q = self.delta[k], self.q[k]
        return q ** 2 / (1 - d.conjugate() * z) ** 2

    def circles(self):
        r'''List of the 2M circles (centre, radius): C_1..C_M, C_1'..C_M\'.'''
        return list(zip(self.delta + self.deltap, self.q + self.qp))


def basis(D, N1, z):
    r'''The 2*M*N1 functions (q_k/(z-delta_k))^m and (qtilde_k/(z-delta_k'))^m, in this order
    (k outer loop, m inner loop), evaluated at z.'''
    out = []
    for dd, qq in list(zip(D.delta, D.q)) + list(zip(D.deltap, D.qtilde)):
        u = qq / (z - dd)
        pw = u
        for m in range(N1):
            out.append(pw)
            pw = pw * u
    return out


def circle_points(c, r, n, K, phase):
    I = K(0, 1)
    return [c + r * (2 * K(pi) * I * (l + phase) / n).exp() for l in range(n)]


def lsq_real(A, b):
    r'''Least squares solution of A x = b for a real ball matrix A, through the normal equations.'''
    At = A.transpose()
    G = At * A
    rhs = At * b
    return G.solve_right(rhs)


def lsq_complex(A, b):
    r'''Idem for a complex matrix, A^H A x = A^H b.'''
    Ah = A.conjugate_transpose() if hasattr(A, 'conjugate_transpose') else A.transpose().conjugate()
    G = Ah * A
    rhs = Ah * b
    return G.solve_right(rhs)


def step1(D, N1, prec, N2=None, phase=QQ(37)/100):
    r'''Coefficients (a, b) of the integrals v_j, j = 1..M, see (12) of the paper of Crowdy--Marshall.
    The matrix of the least-squares problem does not depend on j, so it is factored once.  Returns
    a list (indexed by j) of lists of the 2*M*N1 complex coefficients, in the order of `basis`.'''
    K = ComplexBallField(prec)
    R = RealBallField(prec)
    M = D.M
    N2 = N2 or (2 * N1 + 1)
    circs = D.circles()
    ncoef = 2 * M * N1
    rows = []
    rhs = [[] for _ in range(M)]
    for m, (c, r) in enumerate(circs):
        for z in circle_points(c, r, N2, K, phase):
            Bz = basis(D, N1, z)
            # Im(a B) = Re(a) Im(B) + Im(a) Re(B)
            row = [x.imag() for x in Bz] + [x.real() for x in Bz] + [R(0)] * (2 * M)
            row[2 * ncoef + m] = R(-1)
            rows.append(row)
            for j in range(M):
                # Im( (1/2 pi i) log((z-dj)/(z-dpj)) ) = -(1/2 pi) log|...|
                ratio = (z - D.delta[j]) / (z - D.deltap[j])
                rhs[j].append(ratio.abs().log() / (2 * R(pi)))
    A = Matrix(R, rows)
    B = Matrix(R, rhs).transpose()
    At = A.transpose()
    X = (At * A).solve_right(At * B)
    return [[K(X[i, j], X[ncoef + i, j]) for i in range(ncoef)] for j in range(M)]


def eval_v(D, N1, j, coef, z, K):
    r'''v_j(z) from its coefficients.'''
    ans = ((z - D.delta[j]) / (z - D.deltap[j])).log() / (2 * K(pi) * K(0, 1))
    Bz = basis(D, N1, z)
    return ans + sum(c * b for c, b in zip(coef, Bz))


def prime_square(D, N1, prec, targets, w, N2=None):
    r'''Returns [X(z, w) for z in targets] together with the data (setup times are measured by the caller).'''
    K = ComplexBallField(prec)
    M = D.M
    N2 = N2 or (2 * N1 + 1)
    vcoef = step1(D, N1, prec, N2)
    I = K(0, 1)
    PI = K(pi)
    # tau_jj = v_j(theta_j(z0)) - v_j(z0), z0 on C_j'
    tau = []
    for j in range(M):
        z0 = D.deltap[j] + D.qp[j] * K(3, 9) / 10
        z0 = D.deltap[j] + D.qp[j] * (z0 - D.deltap[j]) / (z0 - D.deltap[j]).abs()
        tau.append(eval_v(D, N1, j, vcoef[j], D.theta(j, z0), K) - eval_v(D, N1, j, vcoef[j], z0, K))
    ncoef = 2 * M * N1
    rows, rhs = [], []
    vw = [eval_v(D, N1, j, vcoef[j], w, K) for j in range(M)]
    for j in range(M):
        for z in circle_points(D.deltap[j], D.qp[j], N2, K, QQ(11)/100):
            tz = D.theta(j, z)
            E = (-2 * PI * I * (2 * (eval_v(D, N1, j, vcoef[j], z, K) - vw[j]) + tau[j])).exp()
            f1 = (tz - w) ** 2
            f2 = E * D.dtheta(j, z) * (z - w) ** 2
            B1, B2 = basis(D, N1, tz), basis(D, N1, z)
            rows.append([f1 * a - f2 * b for a, b in zip(B1, B2)])
            rhs.append(f2 - f1)
    A = Matrix(K, rows)
    x = lsq_complex(A, vector(K, rhs))

    def Xhat(z):
        return 1 + sum(c * b for c, b in zip(x, basis(D, N1, z)))
    Aconst = 1 / Xhat(w)
    return [(z - w) ** 2 * Aconst * Xhat(z) for z in targets]
