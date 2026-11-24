#!/usr/bin/env python3
"""A second step for a base of the perfect cuboid, after the descent with rational classes (EXPERIMENT).

Base: a Pythagorean triple (a, b, c) = (m^2 - n^2, 2 m n, m^2 + n^2).  A perfect cuboid on it gives t = u/v with
t^2 + a^2, t^2 + b^2, t^2 + c^2 squares.  First step (tuple_descent.py cuboid): u + i x v = d_x (square) in Q(i) for
x = a, b, c, with squarefree positive integers d_x, and the tuples (d_a, d_b, d_c) that pass every local test.

Second step.  Put s = i t.  Then s - x = i (u + i x v)/v and s + x = i (u - i x v)/v, so for {x, y, z} = {a, b, c} and
signs e, f the number  (s^2 - x^2)(s - e y)(s - f z)  is  d_y d_z  times a square of Q(i) (four factors: the denominator
and the unit go away).  So the solution gives a point over Q(i), with s = i t, on the curve of genus 1
    C:  Y^2 = D (s - x)(s + x)(s - e y)(s - f z),      D = the squarefree part of d_y d_z,
which is defined over Q and has rational points.  With X = k/(s - x), k = D c1 c2 c3, c1 = 2x, c2 = x - e y, c3 = x - f z,
it is  E: y^2 = (X + D c2 c3)(X + D c1 c3)(X + D c1 c2).  Let P be the point of E(Q(i)), and P' its conjugate.
If E(Q) has rank 0, P + P' is a torsion point; if the twist of E by -1 has rank 0, P - P' is one.  In both cases
P' = +-P + T with T a torsion point of E(Q(i)).  Since the conjugate of s is -s, X(P') = -k X/(2 x X + k), and the relation
becomes a polynomial equation for X(P) for each T; its roots, and the X of the torsion points, give all possible s.
A tuple is closed when some curve does this and every rational t found is t = 0 or infinity (or is not a solution);
a base is closed when every tuple left by the first step is closed.  A rational t that is a solution would be printed
as SOLUTION: it would be a perfect cuboid.

What is trusted: PARI's upper bound 0 for the rank (ellrank), Sage's torsion subgroup over Q(i) (its order is compared
with PARI's elltors over Q(i) on every curve used) and its roots of polynomials over Q(i), and the first step.  usage: second_step.py m,n [m,n ...]"""
import sys, os
from math import prod, isqrt
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tuple_descent as td
from sage.all__sagemath_schemes import *
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
K = QuadraticField(-1, 'i'); I = K.gen(); RK = PolynomialRing(K, 'X'); Xv = RK.gen()
nfK = pari('nfinit(y^2 + 1)')

def issq(q):
    q = Fr(q); return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator
def is_solution(t, xs): return t != 0 and all(issq(Fr(t)**2 + x*x) for x in xs)

_rk = {}
def rank0(e):
    key = tuple(e)
    if key not in _rk:
        a2 = sum(e); a4 = e[0]*e[1] + e[0]*e[2] + e[1]*e[2]; a6 = e[0]*e[1]*e[2]
        _rk[key] = int(pari.ellrank(pari.ellinit([0, a2, 0, a4, a6]), 1)[1]) == 0
    return _rk[key]

def candidates(D, x, r1, r2):
    """All rational t that a point with P' = +-P + T can have, on Y^2 = D (s - x)(s + x)(s - r1)(s - r2); None if some
    equation is an identity."""
    c = [2*x, x - r1, x - r2]; k = D*prod(c); e = [D*c[1]*c[2], D*c[0]*c[2], D*c[0]*c[1]]
    E = EllipticCurve(K, [0, sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], prod(e)]); a2 = K(sum(e))
    f = (Xv + e[0])*(Xv + e[1])*(Xv + e[2]); den = 2*x*Xv + k; out = set()
    def take(Xval):
        if Xval == 0: out.add('infinity'); return
        s = x + k/Xval; t = s/I
        if t in QQ: out.add(QQ(t))
    tors = E.torsion_subgroup().list()
    # a second program for the torsion over Q(i): PARI's elltors has to give the same order as Sage
    assert int(pari.elltors(pari.ellinit([0, sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], prod(e)], nfK))[0]) == len(tors), ('torsion', D, x, r1, r2)
    for T in tors:
        T = E(T)
        if T.is_zero(): eq = Xv*den + k*Xv                    # X(P') = X(P)
        else:
            XT, yT = T[0], T[1]; take(XT)
            # (y -+ yT)^2 = (X - XT)^2 (a2 + X + XT + R), R = -k X/den;  H := the right side times den
            H = (Xv - XT)**2*((a2 + Xv + XT)*den - k*Xv)
            eq = H - f*den if yT == 0 else 4*yT**2*f*den**2 - (H - (f + yT**2)*den)**2
        if eq == 0: return None
        for root, _ in eq.roots():
            if den(root) != 0: take(root)
    return out

def base(m, n):
    a, b, c = m*m - n*n, 2*m*n, m*m + n*n; xs = (a, b, c)
    forms, sup, known = td.cuboid(m, n); S, used, sizes = td.descend(forms, sup, known, 150)
    curves = [(i, j, k, ej, ek) for i in range(3) for (j, k) in [tuple(q for q in range(3) if q != i)] for ej in (1, -1) for ek in (1, -1)]
    closed = 0; left = []; sols = set(); seen = set()
    for tup in S:
        done = False
        for (i, j, k, ej, ek) in curves:
            D = td.sqf(tup[j]*tup[k]); x, r1, r2 = xs[i], ej*xs[j], ek*xs[k]; cc = [2*x, x - r1, x - r2]
            e = [D*cc[1]*cc[2], D*cc[0]*cc[2], D*cc[0]*cc[1]]
            if not (rank0(e) or rank0([-v for v in e])): continue
            cand = candidates(D, x, r1, r2)
            if cand is None: continue
            seen |= cand; bad = [t for t in cand if t != 'infinity' and is_solution(t, xs)]
            sols |= set(bad); done = True; break
        if done: closed += 1
        else: left.append(tup)
    return dict(base=(a, b, c), tuples=len(S), closed=closed, left=left, t=sorted(str(v) for v in seen), solutions=sorted(str(v) for v in sols))

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        m, n = (int(v) for v in arg.split(',')); r = base(m, n)
        status = 'SOLUTION' if r['solutions'] else ('closed' if not r['left'] else 'open')
        print('(%d,%d) base %s: tuples %d, closed %d; rational t found: %s; %s%s' % (m, n, r['base'], r['tuples'], r['closed'], ' '.join(r['t']) or 'none', status,
              '' if not r['left'] else ' (tuples without a curve: %s)' % r['left'][:4]), flush=True)
