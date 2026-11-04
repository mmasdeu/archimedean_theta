# archimedean_theta

Polynomial-time computation of Poincaré theta series of classical Schottky groups: the Abelian integrals
Ω<sub>n</sub>, the holomorphic differentials ω<sub>n</sub>, the period matrix and the Schottky–Klein prime function,
to hundreds of digits.

This is the code accompanying the article

> M. Masdeu, *A polynomial-time algorithm to compute theta functions of Schottky groups*.

## The idea in one paragraph

The series Ω<sub>n</sub>(z) = Σ<sub>σ∈G/⟨σ<sub>n</sub>⟩</sub> log((z − σB<sub>n</sub>)/(z − σA<sub>n</sub>)) is a sum over a free group.
Instead of summing it term by term (which needs a number of terms exponential in the number of digits), we group the
terms according to the first letter of the word. Each group is a function holomorphic outside one of the Schottky discs,
and these 2N functions satisfy an exact linear equation F = E + 𝒯F, where 𝒯 is the composition with the generators
(a transfer operator). Expanding each function in a Taylor series in a coordinate adapted to its disc turns this into a
finite linear system, which is solved by iteration or directly. For a fixed Schottky datum the cost is O(N³P³) arithmetic
operations for P bits of precision and genus N, and the convergence is geometric.

## Quick start

```python
# run from the repository root with:  sage examples/genus2_basic.sage
import os, sys; sys.path.insert(0, os.getcwd())
from schottky_theta import *

K = ComplexField(200); i = K.gen()
discs = {1: Disc(-2-2*i, K(1)), -1: Disc(-2+2*i, K(1)),      # sigma_i maps the complement of
         2: Disc( 2-2*i, K(1)), -2: Disc( 2+2*i, K(1))}      # discs[-i] onto the closure of discs[i]
S = SchottkyTheta(discs, 60)          # 60 Taylor coefficients per disc
F = S.solve()                         # all the functions Omega_n at once
z = K(3)/10 + K(1)/10*i
S.Omega(1, F.column(0), z)            # Abelian integral Omega_1(z)
S.omega(1, F.column(0), z)            # holomorphic differential omega_1 (coefficient of dz)
S.period_matrix(F)                    # period matrix B (Riemann matrix tau = B / (2 pi i))
```

Your own generators can be passed with `SchottkyTheta(discs, p, generators={1: M1, 2: M2, ...})`, where `Mi` is a
2×2 matrix mapping the complement of `discs[-i]` onto the closure of `discs[i]`. See [INSTALL.md](INSTALL.md) for
installation and more examples.

## Contents

| Path | Contents |
|------|----------|
| `schottky_theta.py` | Self-contained implementation (Sage only): transfer matrices (three-term recurrence), iteration and direct solver, Ω<sub>n</sub>, ω<sub>n</sub>, period matrix; baselines (word-length truncation, pruned depth-first summation). |
| `theta_arch.py` | First implementation, iterating on exponentials of power series; used for the Schottky–Klein prime function (`SKprime`). Needs [darmonpoints](https://github.com/mmasdeu/darmonpoints). |
| `examples/genus2_basic.sage` | Basic usage. |
| `examples/crowdy_marshall_prime.sage` | The prime function of a triply connected domain (Crowdy–Marshall 2007), to 30 digits. |
| `examples/theta_arch_genus2.sage`, `examples/fro_family_theta_arch.sage` | The first implementation on two genus 2 examples. |
| `examples/fro/` | Example 8.1 of the paper: the family of Fairchild–Ríos Ortiz, a genus 2 curve to 109 digits and an integer relation search on its Igusa invariants; their genus 3 Riemann matrix. |
| `examples/x0_23/` | Example 8.2: a Schottky uniformization of the modular curve X<sub>0</sub>(23) to 40 digits. |
| `benchmarks/` | Scripts and outputs behind the tables and figures of Section 7 (see the table below). |

### Scripts behind the paper

| Paper | Script | Output |
|-------|--------|--------|
| §7.2 correctness checks | `benchmarks/correctness.sage` | `benchmarks/correctness.out` |
| §7.2 double coset formula (Remark 5.2) | `benchmarks/check_period_formula.sage` | printed |
| §7.2 error on the boundary circles | `benchmarks/boundary.sage` | `benchmarks/boundary.out` |
| Table 2, Figure 1 | `benchmarks/server/{mkref,newalg,pruned}.sage` | `benchmarks/server/*.out` |
| Table 3, Figure 2 | `benchmarks/server/genus.sage`, `benchmarks/genus_p64.sage` | `benchmarks/server/genus.out`, `benchmarks/genus_p64.out` |
| Figure 3 | `benchmarks/geometry.sage` | `benchmarks/geometry.out` |
| §7.6 prime function | `benchmarks/server/skprime.sage` | `benchmarks/server/skprime.out` |
| Remark after Theorem 4.4 | `benchmarks/recurrence_check.sage` | `benchmarks/recurrence_check.out` |
| Example 8.1 | `examples/fro/step*.sage` | `examples/fro/*.out` |
| Example 8.2 | `examples/x0_23/*.sage` | `examples/x0_23/*.out` |

All scripts are meant to be run **from the repository root**, e.g. `sage benchmarks/correctness.sage`.

## Citing

If you use this code, please cite the article above.

## License

Copyright (C) 2026 Marc Masdeu. The code is free software, distributed under the
[GNU General Public License, version 3](LICENSE) or (at your option) any later version. It comes with
no warranty. The text of the license is in the file `LICENSE`. `theta_arch.py` contains functions translated from
MIT-licensed Julia code of Fairchild and Ríos Ortiz; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
