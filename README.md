# Perfect cuboid: a census by base rectangle

A perfect cuboid is a box with integer edges `a, b, c` whose three face diagonals and space diagonal are integers. Whether one exists is an open problem. This repository takes one base rectangle at a time and records which bases can be excluded by an elliptic curve of rank 0. It is a first step, made with the tools of [magic-square-of-squares-lean](https://github.com/drew2apps/magic-square-of-squares-lean), and **the method is not new** (see "What is known").

**Status.** Started on 23 November 2026. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`). Nobody has reviewed it. Nothing is in Lean yet. It rests on PARI's `ellrank` and `elltors`.

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

## What is known

- Excluding cuboids with a given face ratio by an elliptic curve of rank 0 is in the literature, for instance in *Novel required properties of, and efficient algorithms to seek, perfect cuboids* (arXiv:2401.06784), and R. Peschmann's 2026 preprints treat explicit families of fibres with elliptic curves and torsion (arXiv:2604.28072, arXiv:2605.00573). We have not compared the list of closed bases here with theirs. If this census has any value of its own, it is as a reproducible table.
- The searches by size (no perfect cuboid with small edges) are a different kind of statement.

## What could be wrong

- `ellrank` returns an upper bound from a 2-descent; the statement "rank 0" is PARI's. No second program has been run.
- The torsion points are taken from `elltors`; the program checks that it lists as many points as the order of the group.
- The reduction from a perfect cuboid on the base to the four curves is elementary and is written in the docstring of `base_census.py`; nobody else has read it.

## Next

- For an open base with a curve of rank exactly 1 and a known generator `G`, every rational point is `kG + T`, and the conditions can be tested modulo many primes for each `k`: a sieve that would show that a perfect cuboid on that base, if any, has `|k|` beyond a large bound (the same idea as `mw_sieve.py` in the magic-square repository).
- A second rank program for the closed bases, and `L`-values or isogeny descents for the 171 undetermined ones.
- Lean statements of the reduction.

## Reproduce

```
pip install -r requirements.txt
cd scripts
(python3 base_census.py 2 75; python3 base_census.py 76 100) > base_census_100.jsonl      # the committed file was made in these two parts, about 15 minutes
python3 summary.py
```

CI reruns the bases with `m ≤ 30` and compares them with the committed file.
