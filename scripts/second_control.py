#!/usr/bin/env python3
"""A control for second_step.py on real points (EXPERIMENT).

Triples (x1, x2, x3) of positive integers, not Pythagorean, for which a rational t with t^2 + x_j^2 all squares IS known:
three rational numbers x = (p^2 - q^2)/(2 p q) (so that 1 + x^2 is a square), scaled to integers by their common
denominator L, with t = L.  For such a triple the argument of second_step.py applies word for word: the point has a
tuple (d_1, d_2, d_3), and on each of the twelve curves it gives a point over Q(i) with s = i t.  Whenever one of the two
ranks is 0, the function candidates() has to return t or -t.  The program counts the curves where a rank is 0 and the
curves where the planted t is missing (there must be none).  In practice the equations have only degenerate roots, so a
planted point forces both ranks to be positive and the first count is 0; to show that this is not an accident of the
family, the same count is made with a wrong class (3 D, 5 D, 7 D in place of D), where rank 0 is common.
usage: second_control.py [number of triples, default 40] [seed]"""
import sys, random
from math import gcd
from fractions import Fraction as Fr
import second_step as ss
td = ss.td

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40; rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
triples = tested = found = missing = identities = curves = 0; wrong = {3: 0, 5: 0, 7: 0}
while triples < N:
    xs = set()
    while len(xs) < 3:
        p, q = rnd.randrange(2, 9), rnd.randrange(1, 8)
        if p > q and gcd(p, q) == 1: xs.add(Fr(p*p - q*q, 2*p*q))
    xs = sorted(xs); L = 1
    for x in xs: L = L*x.denominator//gcd(L, x.denominator)
    xi = [int(x*L) for x in xs]; t = L
    assert ss.is_solution(t, xi)
    tup = [td.cls_of((t, x)) for x in xi]                     # u + i x v at (u, v) = (t, 1)
    triples += 1
    for i in range(3):
        j, k = [q for q in range(3) if q != i]
        for ej in (1, -1):
            for ek in (1, -1):
                D = td.sqf(tup[j]*tup[k]); x, r1, r2 = xi[i], ej*xi[j], ek*xi[k]; cc = [2*x, x - r1, x - r2]
                e = [D*cc[1]*cc[2], D*cc[0]*cc[2], D*cc[0]*cc[1]]; curves += 1
                for w in wrong:
                    ew = [w*v for v in e]; wrong[w] += ss.rank0(ew) or ss.rank0([-v for v in ew])
                if not (ss.rank0(e) or ss.rank0([-v for v in e])): continue
                cand = ss.candidates(D, x, r1, r2)
                if cand is None: identities += 1; continue
                tested += 1
                if t in cand or -t in cand: found += 1
                else: missing += 1; print('MISSING', xi, t, (i, ej, ek), flush=True)
print('%d planted triples, %d curves; with a rank equal to 0 and an equation that is not an identity: %d; the planted t is among the candidates: %d; missing: %d; identities: %d' % (triples, curves, tested, found, missing, identities))
print('with a wrong class, curves with a rank equal to 0: %s' % ', '.join('%d D: %d of %d' % (w, wrong[w], curves) for w in wrong))
