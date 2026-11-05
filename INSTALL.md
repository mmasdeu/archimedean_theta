# Installation

## Requirements

* [SageMath](https://www.sagemath.org) (tested with 10.6 and 10.9). `schottky_theta.py` needs nothing else.
* For `theta_arch.py` (the first implementation, used for the Schottky–Klein prime function): the
  [darmonpoints](https://github.com/mmasdeu/darmonpoints) package and `tqdm`.
* For `examples/x0_23/check_tau.sage`: Sage's `RiemannSurface`, which needs `mpmath` (included in standard Sage
  installations).

## Installing Sage

The simplest way is through conda (see the [Sage installation guide](https://doc.sagemath.org/html/en/installation/)):

```bash
conda create -n sage -c conda-forge sage python=3.11
conda activate sage
```

Then, for the examples that use `theta_arch.py`:

```bash
sage -pip install darmonpoints tqdm
```

## Getting the code

```bash
git clone https://github.com/mmasdeu/archimedean_theta.git
cd archimedean_theta
```

No installation step is needed: the scripts add the repository root to the Python path, so they must be **run from the
repository root**. To use the modules from elsewhere, add the repository to `PYTHONPATH` (or copy
`schottky_theta.py` next to your script).

## Examples

| Command | What it does | Time |
|---------|--------------|------|
| `sage examples/genus2_basic.sage` | Ω<sub>1</sub>, ω<sub>1</sub>, the period matrix and the Riemann matrix of a genus 2 Schottky group, and a comparison with the truncated series | seconds |
| `sage examples/crowdy_marshall_prime.sage` | The square of the Schottky–Klein prime function of a triply connected domain, compared with the value of Crowdy–Marshall (needs darmonpoints) | ~15 s |
| `sage examples/collocation_cm2007.sage` | The same value with the collocation method of Crowdy–Marshall (`collocation.py`) | ~30 s |
| `sage examples/theta_arch_genus2.sage` | The first implementation against the truncated series (needs darmonpoints) | ~30 s |
| `sage examples/fro/step2.sage` | Riemann matrix of a genus 3 curve of the Fairchild–Ríos Ortiz family | ~40 s |
| `sage examples/fro/step3.sage 1/12 1/4 60` | Weierstrass values and Igusa–Clebsch invariants of the genus 2 curve S<sub>π/12,π/4</sub> (use 200 or 240 for the 109 digits of the paper, 12–21 minutes) | ~15 s |
| `sage examples/fro/step4.sage`, `sage examples/fro/step5.sage` | Integer relation searches on these invariants (use the stored `phis_200.sobj`, `phis_240.sobj`) | ~1 min |
| `sage examples/x0_23/verify_equation.sage` | Checks the model of X<sub>0</sub>(23) against the newform of level 23 | seconds |
| `sage examples/x0_23/newton.sage` | Computes the Schottky group of X<sub>0</sub>(23) (Newton's method, up to p = 200) | ~50 min |
| `sage examples/x0_23/report.sage` | Prints the parameters, fixed points and multipliers from the stored solution | ~1 min |
| `sage examples/x0_23/check_tau.sage 120` | Compares its Riemann matrix with the one computed by Sage's `RiemannSurface` | ~2 min |
| `sage benchmarks/correctness.sage` | The correctness checks of Section 7.2 of the paper | ~1 min |

Timings are for one core of a desktop computer.

### Using your own Schottky group

```python
import os, sys; sys.path.insert(0, os.getcwd())
from schottky_theta import *
K = ComplexField(300)
discs = {...}          # discs[i] = Disc(center, radius) for i = +-1, ..., +-N; closures pairwise disjoint,
                       # infinity outside all of them
gens  = {...}          # gens[i] (i = 1..N): 2x2 matrix over K mapping the complement of discs[-i]
                       # onto the closure of discs[i]
S = SchottkyTheta(discs, p=80, generators=gens)
F = S.solve()          # or, for one function: S.iterate(n, steps)
tau = S.period_matrix(F) / (2*K.pi()*K.gen())
```

Choose `p` by increasing it until two consecutive values give the same result to the digits you need; the number of
correct digits per Taylor coefficient depends on how close the discs are to each other (about 1 for well separated
discs, 0.2–0.3 when they nearly touch). Use a working precision of about 3.33·(digits) + 60 bits. The series converge
under the Schottky-type hypothesis of the paper (for example, when the group preserves a circle).
