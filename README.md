# Perfect cuboid: a census by base rectangle

A perfect cuboid is a box with integer edges `a, b, c` whose three face diagonals and space diagonal are integers. Whether one exists is an open problem. This repository takes one base rectangle at a time and records which bases can be excluded by an elliptic curve of rank 0. It is a first step, made with the tools of [magic-square-of-squares-lean](https://github.com/drew2apps/magic-square-of-squares-lean), and **the method is not new** (see "What is known").

**Status.** Started on 23 November 2026. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). Nobody has reviewed it. Nothing is in Lean yet. It rests on PARI's `ellrank` and `elltors`; the descent in two steps, added on 24 November, also on Sage and on two programs that nobody has read.

## The setting

Fix a primitive Pythagorean triple `(a, b, d) = (m² − n², 2mn, m² + n²)`, with `m > n ≥ 1` coprime of different parity: the base rectangle and its diagonal. Up to scaling, a perfect cuboid on that base is a rational `t = c ≠ 0` at which

    t² + a²,    t² + b²,    t² + d²

are all squares. With `X = t²`, such a `t` gives a rational point with `X` a nonzero square on each of four elliptic curves:

| Curve | Equation | Meaning |
|---|---|---|
| `E_ab` | `Y² = X(X + a²)(X + b²)` | the Euler bricks on the base |
| `E_ad` | `Y² = X(X + a²)(X + d²)` | |
| `E_bd` | `Y² = X(X + b²)(X + d²)` | |
| `E_2` | `Y² = X(X + a²)(X + b²)(X + d²)` | a quartic with the point `X = 0` |

(These are four of the five elliptic curves in the Jacobian of the curve of genus 5 that the three conditions define. The fifth, `Y² = (X + a²)(X + b²)(X + d²)`, always has the point `X = 0` of infinite order.) If one of the four has rank 0, its rational points are its torsion points, a finite list; if none of them has `X` a nonzero square, **there is no perfect cuboid on that base**, whatever its height.

## Result

`scripts/base_census.py`, output in `scripts/base_census_100.jsonl` (one line per base, with the rank bounds of the four curves), counts in `scripts/summary_100.txt`:

- There are 2,040 primitive bases with `m ≤ 100`.
- **1,629 are closed**: one of the four curves has rank 0 (PARI's upper bound is 0) and no torsion point has `X` a nonzero square. For `m ≤ 20` it is 83 of 86.
- 411 are open; the first are (15,8), that is the base 161 × 240 with diagonal 289, then (17,16), (19,16), (22,1).
- Each of the four curves is the only one that closes some bases (between 178 and 214 each), so all four are needed.
- For 171 of the open bases PARI does not determine some rank (lower bound below upper bound), so a rank 0 hidden behind a Tate–Shafarevich group is possible there.

A closed base is closed for every multiple of the triple, since the conditions are homogeneous. The census is by base, not by size: it says nothing about a bound on the edges, and the published searches by size say nothing about a base.

## The open bases: a descent in two steps (experiment, counted apart)

`scripts/tuple_descent.py`, `scripts/second_step.py`; output in `scripts/second_step_open_100.txt`, counts in `scripts/second_step_summary.txt`. This part is new in this repository, it is two programs that nobody has read, and its bases are counted apart from the 1,629.

**First step: three classes.** Let `t = u/v`, in lowest terms, be a solution. For `x = a, b, d` the number `u² + x²v²` is a square, and it is the norm of `u + ixv` in `ℤ[i]`. A Gaussian integer of square norm is a squarefree positive integer `δ_x` times a square of `ℚ(i)`, and an odd prime that divides `δ_x` divides `x` (a prime above it divides `u + ixv` and its conjugate, hence `2u` and `2xv`; this uses unique factorisation in `ℤ[i]`). So a solution has one of finitely many triples `(δ_a, δ_b, δ_d)`. The program keeps the triples that are possible at every prime it tests (the primes of `2abd` and all primes up to 150): it covers the projective line over `ℤ_p` by discs and finds the square class of each form on each disc; what it cannot decide counts as possible. The triples of `t = 0` and `t = ∞` always remain. On the 411 open bases between 3 and 36 triples remain, 6 in the median.

**Second step: a curve of genus 1 over `ℚ(i)`.** Put `s = it`. Then `s − x = i(u + ixv)/v` and `s + x = i(u − ixv)/v`. So for `{x, y, z} = {a, b, d}` and signs `e, f` the number `(s² − x²)(s − ey)(s − fz)` is `δ_y·δ_z` times a square of `ℚ(i)`: with four factors the unit `i` and the denominator `v` go away. The solution is therefore a point over `ℚ(i)`, with `s = it`, on

    C:  Y² = D·(s − x)(s + x)(s − ey)(s − fz),      D = the squarefree part of δ_y·δ_z,

a curve of genus 1 defined over `ℚ` with the rational point `s = x`; there are twelve such curves for each triple. Let `E` be its Jacobian, `P` the point and `P′` its conjugate. `P + P′` lies in `E(ℚ)`, and `P − P′` in the group of the twist of `E` by `−1`. If one of these two groups has rank 0, then `P′ = ±P + T` with `T` a torsion point of `E(ℚ(i))`. The conjugate of `s = it` is `−s`, so this is a polynomial equation for the coordinate of `P`, one for each `T`, and its roots in `ℚ(i)` give every possible `t`. A triple is closed when one of its twelve curves has such a rank 0 and the only rational `t` that come out are `0` and `∞`; a base is closed when every triple left by the first step is closed.

**Result.**

- **274 of the 411 open bases are closed**: every triple of the first step lies on a curve with a rank 0, and the only rational `t` found, on all curves, are `0` and `∞`. Among them (15,8) and (17,16), the first two open bases of the census; (19,16) stays open.
- 137 stay open: some triple has no curve with a rank 0. Of 3,427 triples in all, 3,145 are closed.
- With the sieve of the next section: of the 137, 123 have no `k` up to 10⁶, and 14 have nothing.

**Controls.**

- *Planted points* (`second_control.py`, `second_control.txt`). For 40 triples `(x₁, x₂, x₃)` that are not Pythagorean and for which a rational `t` with all three `t² + x_j²` square is known, the same construction is made with the classes of that `t`: 480 curves. With the right class none of them has a rank 0, as it must be (the equations have only degenerate roots, so a real point forces both ranks to be positive). With a wrong class (`3D`, `5D`, `7D`) 224, 224 and 222 of the 480 have a rank 0. So the curves and classes are the right ones, and a rank 0 is not rare by itself.
- *First step* (`first_step_15_8.txt`): 33,429 random local points and 9,342 points placed near the zeros of the forms, none outside the allowed sets; the class rules against lists of all squares modulo `2⁶, 3⁴, 5³, 7³, 13²`, no disagreement in 25,824 values.
- *First step by exhaustion* (`tuple_bruteforce.py`, `bruteforce_15_8.txt`): a second program with no code in common goes through every point modulo `p^N` with lists of squares. On the bases (15,8), (17,16), (19,16) it gives the same set as the descent at every prime up to 61. It also found that an earlier version of the descent was weaker than necessary (33 triples on (15,8) where 9 remain now); with that version the second step closed 219 bases, all of them among the 274.
- *Torsion over `ℚ(i)`*: the order given by Sage is compared with PARI's `elltors` on every curve used.

**What it rests on.** The two paragraphs above, which nobody has read; PARI's upper bound for the rank (`ellrank`), with no second program; Sage (passagemath 10.8) for the torsion over `ℚ(i)` and the roots of the equations. The first step is the program of the magic-square repository (`LAST_PAIR.md`, sections 12 and 14 there, for the rules and their proofs).

## The open bases: a sieve

`scripts/open_sieve.py`, output in `scripts/open_sieve_100.jsonl`, counts in `scripts/open_sieve_summary.txt`. For an open base, take one of the four curves with rank exactly 1 and a generator `G` (PARI finds one for 369 of the 411). Every rational point is `kG + T` with `T` a torsion point. A perfect cuboid gives a point at which `X`, `X + a²`, `X + b²`, `X + d²` are all squares; modulo a prime this depends on `k` modulo the order of `G`, so each prime gives a table of the `k` it allows, and a few dozen primes leave nothing.

- **For 369 open bases no `k` with `0 < k ≤ 1,000,000` survives.** A perfect cuboid on one of these bases would be a point `kG + T` with `|k| > 10⁶`; since the height of `kG` grows like `k²` times the height of `G` (at least 2.15), its third edge would have more than `10¹¹` digits (this conversion is a heuristic; the statement about `k` is what is computed).
- 42 open bases have no curve of rank exactly 1 with a known point, and are not treated.
- Control (`open_sieve.py --control 11,2 100000`): on the base 117 × 44 the Euler brick 44 × 117 × 240 is the point `2G` of `E_ab`. With the condition on `X + d²` switched off the sieve leaves every even `k` (the Euler bricks on that base), as it must; with it switched on, nothing.
- This is a bound for each base, not a closure. It assumes that `G` generates the group modulo torsion (PARI's point, divided by 2, 3, 5, 7 where possible; a divisor of 11 or more is not excluded).

So for the 2,040 bases with `m ≤ 100`: 1,629 have no perfect cuboid by a curve of rank 0 over `ℚ`; 274 more by the descent in two steps, if that experiment is right; of the other 137, 123 have none below an absurd size and 14 are untouched.

## What is known

- Excluding cuboids with a given face ratio by an elliptic curve of rank 0 is in the literature, for instance in *Novel required properties of, and efficient algorithms to seek, perfect cuboids* (arXiv:2401.06784), and R. Peschmann's 2026 preprints treat explicit families of fibres with elliptic curves and torsion (arXiv:2604.28072, arXiv:2605.00573). We have not compared the list of closed bases here with theirs. If this census has any value of its own, it is as a reproducible table.
- The searches by size (no perfect cuboid with small edges) are a different kind of statement.
- The descent in two steps uses standard tools (a descent on a cover given by norms from `ℤ[i]`, then curves of genus 1 over a quadratic field whose rational part or twist has rank 0). We have not searched the literature for it on this problem and do not claim that it is new.

## What could be wrong

- `ellrank` returns an upper bound from a 2-descent; the statement "rank 0" is PARI's. No second program has been run.
- The torsion points are taken from `elltors`; the program checks that it lists as many points as the order of the group.
- The reduction from a perfect cuboid on the base to the four curves is elementary and is written in the docstring of `base_census.py`; nobody else has read it.

## Next

- The 137 bases left open by the descent in two steps: products of four of the six linear factors that mix the three numbers differently, or curves of higher genus.
- A second program for the ranks used by the second step, and a reader for its two paragraphs.
- An upper bound for `k`, which would turn the sieve into a proof for a base; none is known to us for this setting.
- A second rank program for the closed bases, and `L`-values or isogeny descents for the 171 undetermined ones.
- Lean statements of the reduction.

## Reproduce

```
pip install -r requirements.txt
cd scripts
(python3 base_census.py 2 75; python3 base_census.py 76 100) > base_census_100.jsonl      # the committed file was made in these two parts, about 15 minutes
python3 summary.py
python3 open_sieve.py base_census_100.jsonl 1000000 > open_sieve_100.jsonl      # about a minute
python3 open_sieve.py --control 11,2 100000
pip install -r ../requirements-second.txt
bash second_batch.sh                       # the two-step descent on the 411 open bases, about 20 minutes on two cores; second_step_open_100.txt
python3 second_summary.py second_step_open_100.txt
python3 second_control.py 40 1
```

CI reruns the bases with `m ≤ 30`, the sieve on the first 40 open bases and the control, and, in a second job, the first step with its tests on (15,8), the exhaustion on that base, the second step on the first 30 open bases and the planted control, and compares them with the committed files.
