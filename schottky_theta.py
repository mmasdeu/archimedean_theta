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
Self-contained implementation (Sage) of the iterative algorithm of the paper
for the Abelian integrals Omega_n and the period matrix of a classical
Schottky group.  It only needs Sage; it does not depend on `darmonpoints`.

Conventions (see the paper):
  * `discs[i]`, i in {+-1, ..., +-N}, are open discs D_i with pairwise disjoint
    closures, infinity outside all of them, and sigma_i maps the complement of
    D_{-i} onto the closure of D_i;  sigma_{-i} = sigma_i^{-1}.
  * A_n (repelling) lies in D_{-n}, B_n (attracting) lies in D_n.
  * t_i(z) = r_i / (z - a_i) is the coordinate on the complement of D_i.
'''
from sage.all import (ComplexField, Matrix, vector, log, pi, sqrt, cos, sin)


class Disc:
    def __init__(self, center, radius):
        self.center = center
        self.radius = radius


def mobius(g, x):
    a, b, c, d = g.list()
    return (a * x + b) / (c * x + d)


def generator_from_discs(Dp, Dm):
    r'''Matrix taking the complement of Dm onto the closure of Dp.'''
    a, r, b, s = Dp.center, Dp.radius, Dm.center, Dm.radius
    K = a.parent()
    return (Matrix(K, 2, 2, [1, -a, 0, r]) ** -1 * Matrix(K, 2, 2, [0, 1, 1, 0])
            * Matrix(K, 2, 2, [1, -b, 0, s]))


def fixed_points(g):
    r'''Returns (A, B): repelling and attracting fixed points of g.'''
    a, b, c, d = g.list()
    D = sqrt((a - d) ** 2 + 4 * b * c)
    p, q = ((a - d) + D) / (2 * c), ((a - d) - D) / (2 * c)
    det = a * d - b * c
    if abs(det / (c * p + d) ** 2) < 1:  # p attracting
        return q, p
    return p, q


def coordinate(D):
    a, r = D.center, D.radius
    return Matrix(a.parent(), 2, 2, [0, r, 1, -a])


def substitution_matrix(p, sigma, ti, tj):
    r'''p x p matrix whose m-th row holds the Taylor coefficients (mod u^p)
    of phi(u)^m, where phi = t_j o sigma^{-1} o t_i^{-1}.

    With phi = (a u + b)/(c u + d), f = phi^m satisfies
    (a u + b)(c u + d) f' = m det(phi) f, which gives a three-term
    recurrence: O(p) operations per row, O(p^2) in total.'''
    K = sigma.base_ring()
    a, b, c, d = (tj * sigma.adjugate() * ti.adjugate()).list()
    det = a * d - b * c
    ac, s1, bd = a * c, a * d + b * c, b * d
    M = Matrix(K, p, p)
    M[0, 0] = 1
    f0 = K(1)
    for m in range(1, p):
        f0 = f0 * b / d
        prev, cur = K(0), f0
        M[m, 0] = cur
        for j in range(p - 1):
            nxt = ((m * det - s1 * j) * cur - ac * (j - 1) * prev) / (bd * (j + 1))
            M[m, j + 1] = nxt
            prev, cur = cur, nxt
    return M


def log1p_coefficients(alpha, p):
    r'''Coefficients of log(1 + alpha t) modulo t^p.'''
    v = [alpha.parent()(0)] * p
    pw = alpha.parent()(1)
    for k in range(1, p):
        pw *= alpha
        v[k] = (-1) ** (k + 1) * pw / k
    return v


class SchottkyTheta:
    r'''Precomputes the operator T of the paper truncated at p Taylor
    coefficients.  Usage:

        S = SchottkyTheta(discs, p)
        F = S.solve()                  # all N functions at once
        S.Omega(n, F.column(n-1), z)
        S.period_matrix(F)
    '''

    def __init__(self, discs, p, generators=None):
        r'''`generators[i]` (i = 1..N), if given, must map the complement of
        discs[-i] onto the closure of discs[i]; otherwise they are taken to be
        generator_from_discs(discs[i], discs[-i]).'''
        self.discs = discs
        self.N = max(discs)
        self.p = p
        self.q = p - 1
        self.K = discs[1].center.parent()
        K, N, q = self.K, self.N, self.q
        self.I = list(range(1, N + 1)) + [-i for i in range(1, N + 1)]
        self.idx = {i: k for k, i in enumerate(self.I)}
        self.sigma = {}
        for i in range(1, N + 1):
            if generators is None:
                self.sigma[i] = generator_from_discs(discs[i], discs[-i])
            else:
                self.sigma[i] = Matrix(K, generators[i])
            self.sigma[-i] = self.sigma[i].adjugate()
        self.t = {i: coordinate(discs[i]) for i in self.I}
        self.fixed = {i: fixed_points(self.sigma[i]) for i in range(1, N + 1)}
        n = len(self.I)
        T = Matrix(K, n * q, n * q)
        for i in self.I:
            for j in self.I:
                if j == -i:
                    continue
                M = substitution_matrix(p, self.sigma[i], self.t[i], self.t[j])
                # f = sum_{k>=1} f_k t_j^k  |->  sum_k f_k (row k of M) minus constant
                T[self.idx[i] * q:(self.idx[i] + 1) * q,
                  self.idx[j] * q:(self.idx[j] + 1) * q] = M[1:, 1:].transpose()
        self.T = T

    def seed(self, n):
        K, p, q = self.K, self.p, self.q
        A, B = self.fixed[n]
        E = [K(0)] * (len(self.I) * q)
        for i in self.I:
            if abs(i) == n:
                continue
            P, Q = mobius(self.sigma[i], B), mobius(self.sigma[i], A)
            a, r = self.discs[i].center, self.discs[i].radius
            v1 = log1p_coefficients((a - P) / r, p)
            v2 = log1p_coefficients((a - Q) / r, p)
            for m in range(1, p):
                E[self.idx[i] * q + m - 1] = v1[m] - v2[m]
        return vector(K, E)

    def seeds(self):
        cols = [self.seed(n) for n in range(1, self.N + 1)]
        return Matrix(self.K, cols).transpose()

    def solve(self):
        r'''Solve (1 - T) F = E_n for n = 1..N  (direct method).'''
        n = len(self.I) * self.q
        return (Matrix.identity(self.K, n) - self.T).solve_right(self.seeds())

    def iterate(self, n, steps):
        r'''F = sum_{k <= steps} T^k E_n  (the iteration of the paper).'''
        E = self.seed(n)
        F = E
        for _ in range(steps):
            F = E + self.T * F
        return F

    def Omega(self, n, F, z):
        A, B = self.fixed[n]
        ans = ((z - B) / (z - A)).log()
        for i in self.I:
            tt = self.discs[i].radius / (z - self.discs[i].center)
            if float(abs(tt)) > 1 + 1e-9:
                raise ValueError('z lies inside the disc D_%d; move it into the closure of the fundamental domain first' % i)
            ans += sum(F[self.idx[i] * self.q + m - 1] * tt ** m
                       for m in range(1, self.p))
        return ans

    def omega(self, n, F, z):
        r'''The differential omega_n = dOmega_n, as the coefficient of dz.'''
        A, B = self.fixed[n]
        ans = 1 / (z - B) - 1 / (z - A)
        for i in self.I:
            r = self.discs[i].radius
            tt = r / (z - self.discs[i].center)
            if float(abs(tt)) > 1 + 1e-9:
                raise ValueError('z lies inside the disc D_%d; move it into the closure of the fundamental domain first' % i)
            dF = sum(m * F[self.idx[i] * self.q + m - 1] * tt ** (m - 1)
                     for m in range(1, self.p))
            ans += dF * (-tt ** 2 / r)
        return ans

    def period_matrix(self, F, angle=0.7):
        r'''B_{mn} = Omega_n(sigma_m z) - Omega_n(z), z on the boundary of D_{-m}.'''
        N, K = self.N, self.K
        P = Matrix(K, N, N)
        for m in range(1, N + 1):
            D = self.discs[-m]
            z = D.center + D.radius * K(cos(angle), sin(angle))
            w = mobius(self.sigma[m], z)
            for n in range(1, N + 1):
                P[m - 1, n - 1] = self.Omega(n, F.column(n - 1), w) - self.Omega(n, F.column(n - 1), z)
        return P


def naive_Omega(S, n, z, L):
    r'''Truncation of the defining series to words of length <= L.
    Returns (value, number of terms).'''
    K = S.K
    A, B = S.fixed[n]
    ans, cnt = K(0), 0

    def rec(g, last, depth):
        nonlocal ans, cnt
        if abs(last) != n:
            ans += ((z - mobius(g, B)) / (z - mobius(g, A))).log()
            cnt += 1
        if depth == L:
            return
        for i in S.I:
            if last != 0 and last == -i:
                continue
            rec(g * S.sigma[i], i, depth + 1)

    rec(S.sigma[1] ** 0, 0, 0)
    return ans, cnt


def outer_disc_diameter(gp, D):
    r'''Diameter of the image of the closed disc D under the Moebius map gp
    (assumed to have its pole outside D).'''
    a, b, c, d = gp.list()
    det = a * d - b * c
    x = D.center
    return 2 * D.radius * abs(det) / (abs(c * x + d) ** 2 - abs(c) ** 2 * D.radius ** 2)


def adaptive_Omega(S, n, z, mu):
    r'''Depth-first summation which stops descending below a word sigma as soon as
    diam O(sigma) < mu, where O(sigma) is the outer disc of sigma (Bogatyrev-type
    pruning with K = 1; there is no proven error bound).'''
    K = S.K
    A, B = S.fixed[n]
    ans, cnt = K(0), 0

    def rec(g, last, parent):
        nonlocal ans, cnt
        if abs(last) != n:
            ans += ((z - mobius(g, B)) / (z - mobius(g, A))).log()
            cnt += 1
        if last != 0 and float(outer_disc_diameter(parent, S.discs[last])) < float(mu):
            return
        for i in S.I:
            if last != 0 and last == -i:
                continue
            rec(g * S.sigma[i], i, g)

    rec(S.sigma[1] ** 0, 0, None)
    return ans, cnt
