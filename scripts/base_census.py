#!/usr/bin/env python3
"""Perfect cuboid, one base rectangle at a time: which bases are decided by an elliptic curve of rank 0.

A perfect cuboid has edges a, b, c with a^2+b^2, a^2+c^2, b^2+c^2 and a^2+b^2+c^2 all squares.  Fix the base rectangle: a
primitive Pythagorean triple (a, b, d) = (m^2-n^2, 2mn, m^2+n^2), m > n >= 1 coprime of different parity.  Up to scaling, a
perfect cuboid on that base is a rational t = c != 0 at which the three forms
        f_a = t^2 + a^2,    f_b = t^2 + b^2,    f_d = t^2 + d^2
are all squares.  With X = t^2 this gives rational points with X a nonzero square on four elliptic curves:
        E_ab: Y^2 = X (X + a^2)(X + b^2)     (Y = t * sqrt(f_a f_b); the Euler bricks on the base)
        E_ad: Y^2 = X (X + a^2)(X + d^2)
        E_bd: Y^2 = X (X + b^2)(X + d^2)
        E_2 : Y^2 = X (X + a^2)(X + b^2)(X + d^2), a quartic with the point X = 0; its Jacobian is y^2 = (x + b^2 d^2)(x + a^2 d^2)(x + a^2 b^2), x = a^2 b^2 d^2 / X
(the fifth quotient, Y^2 = (X + a^2)(X + b^2)(X + d^2), always has the point X = 0 of infinite order and is not used).
If one of the four has rank 0 (PARI's ellrank: the upper bound is 0), its rational points are its torsion points; they are
listed, and if none has X a nonzero square, there is no perfect cuboid on this base.  Nothing here is new as a method.
usage: base_census.py M0 M1          one line of JSON per base with M0 <= m <= M1"""
import sys, math, json, itertools, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
def torsion_points(E):
    T = pari.elltors(E); gens = list(T[2]); ords = [int(o) for o in T[1]]; pts = []
    for ks in itertools.product(*[range(o) for o in ords]):
        P = pari('[0]')
        for k, g in zip(ks, gens): P = pari.elladd(E, P, pari.ellmul(E, g, k))
        if len(P) == 2: pts.append(P)
    assert len(pts) + 1 == int(T[0]); return pts
def nonzero_square(X): return X > 0 and bool(pari.issquare(X))
def decide(m, n):
    a, b, d = m*m - n*n, 2*m*n, m*m + n*n; A, B, D = a*a, b*b, d*d
    curves = [('E_ab', [0, A + B, 0, A*B, 0], lambda x: x), ('E_ad', [0, A + D, 0, A*D, 0], lambda x: x), ('E_bd', [0, B + D, 0, B*D, 0], lambda x: x),
              ('E_2', [0, B*D + A*D + A*B, 0, A*B*D*(A + B + D), (A*B*D)**2], lambda x: pari(A*B*D)/x if x != 0 else None)]
    rec = dict(base=[m, n], edges=[a, b, d], rank={}, closed_by=[])
    for name, ai, toX in curves:
        E = pari.ellinit(ai); r = pari.ellrank(E); rec['rank'][name] = [int(r[0]), int(r[1])]
        if int(r[1]) == 0:
            Xs = [toX(P[0]) for P in torsion_points(E)]; sq = [str(X) for X in Xs if X is not None and nonzero_square(X)]
            if sq: rec.setdefault('torsion_with_square_X', {})[name] = sq
            else: rec['closed_by'].append(name)
    rec['status'] = 'closed' if rec['closed_by'] else 'open'
    return rec
if __name__ == '__main__':
    M0, M1 = int(sys.argv[1]), int(sys.argv[2])
    for m in range(M0, M1 + 1):
        for n in range(1, m):
            if math.gcd(m, n) == 1 and (m - n) % 2 == 1: print(json.dumps(decide(m, n)), flush=True)
