#!/usr/bin/env python3
"""The local eliminations of LAST_PAIR.md, section 12, by exhaustion (EXPERIMENT).  Shares no code with tuple_descent.py.

A system is a list of linear forms gamma_j = alpha_j u + beta_j v over Z[i] whose norms have to be squares; at a solution
(u, v) coprime integers) gamma_j = d_j z^2 with d_j a squarefree positive integer made of the primes of a known support.
Then d_j gamma_j(u, v) is the square of an element of Z[i].  For a prime p and M = p^N, (u, v) is a unit lambda times
(t, 1) or (1, t) modulo M, so for some point of the projective line over Z/M and some lambda,
    d_j lambda gamma_j(point)  is a square modulo M in Z[i],  for every j.
The program goes through every point and through lambda in a set of representatives of the units modulo squares
(1 and a non-residue; 1, 3, 5, 7 for p = 2; multiplying lambda by a square changes nothing), with a list of all squares
(x + i y)^2 modulo M, and keeps the tuples (d_j) for which this happens.  A tuple that is not kept is not the tuple of a
solution: a proof by exhaustion.  A value divisible by M counts as a square, so the result can be weaker than that of
tuple_descent.py; N is the largest with p^N below a bound.
usage: tuple_bruteforce.py magic|cuboid m,n [largest prime, default 61] [largest modulus, default 5000]"""
import sys
from itertools import product

def factor(x):
    out = set(); x = abs(x); q = 2
    while q*q <= x:
        while x % q == 0: out.add(q); x //= q
        q += 1
    if x > 1: out.add(x)
    return sorted(out)
def primes_to(b): return [q for q in range(2, b + 1) if all(q % r for r in range(2, int(q**0.5) + 1))]
def products(P):
    out = [1]
    for q in P: out += [c*q for c in out]
    return sorted(out)

def system(kind, m, n):
    if kind == 'magic':     # mv + inu;  nv - imu;  (1 + i)((mu - nv) + i(mv - nu));  (1 + i)((mu + nv) + i(mv + nu))
        forms = [((0, n), (m, 0)), ((0, -m), (n, 0)), ((m + n, m - n), (-(m + n), m - n)), ((m - n, m + n), (n - m, m + n))]
        sup = [factor(2*m*n)]*2 + [factor(2*(m*m - n*n))]*2
    else:                   # u + i x v for the three numbers x of the base
        base = (m*m - n*n, 2*m*n, m*m + n*n); forms = [((1, 0), (0, x)) for x in base]; sup = [factor(2*x) for x in base]
    return forms, sup

def kept(forms, cand, p, M):
    """The set of the product sets (one frozenset of d for each form) that occur at some point and some lambda."""
    sqG = bytearray(M*M)
    for x in range(M):
        xx = x*x
        for y in range(M): sqG[((xx - y*y) % M)*M + (2*x*y) % M] = 1
    lams = [1, 3, 5, 7] if p == 2 else [1, next(r for r in range(2, p) if pow(r, (p - 1)//2, p) == p - 1)]
    boxes = set()
    for u, v in [(t, 1) for t in range(M)] + [(1, t) for t in range(0, M, p)]:
        vals = [((a[0]*u + b[0]*v) % M, (a[1]*u + b[1]*v) % M) for a, b in forms]
        for lam in lams:
            box = tuple(frozenset(d for d in cand[j] if sqG[(d*lam*vals[j][0] % M)*M + d*lam*vals[j][1] % M]) for j in range(len(forms)))
            if all(box): boxes.add(box)
    return boxes

if __name__ == '__main__':
    kind = sys.argv[1]; m, n = (int(x) for x in sys.argv[2].split(',')); pmax = int(sys.argv[3]) if len(sys.argv) > 3 else 61
    Mmax = int(sys.argv[4]) if len(sys.argv) > 4 else 5000
    forms, sup = system(kind, m, n); cand = [products(P) for P in sup]
    alive = None; total = 1
    for c in cand: total *= len(c)
    print('%s (%d,%d): candidates for each form: %s; tuples: %d' % (kind, m, n, [len(c) for c in cand], total))
    for p in primes_to(pmax):
        N = 1
        while p**(N + 1) <= Mmax: N += 1
        boxes = kept(forms, cand, p, p**N)
        if alive is None:
            alive = set()
            for box in boxes: alive |= set(product(*box))
        else: alive = set(t for t in alive if any(all(t[j] in box[j] for j in range(len(t))) for box in boxes))
        print('prime %d, modulo %d^%d: tuples kept at every prime so far: %d%s' % (p, p, N, len(alive), ' ' + str(sorted(alive)) if len(alive) <= 4 else ''))
    print('kept at every prime: %d: %s' % (len(alive), sorted(alive) if len(alive) <= 40 else 'more than 40'))
