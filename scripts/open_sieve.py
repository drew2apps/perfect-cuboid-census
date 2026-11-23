#!/usr/bin/env python3
"""The open bases of base_census.py: a sieve on the multiples of a generator.

For an open base take one of the four curves with rank exactly 1 (PARI: bounds 1 and 1, with a point G of infinite order).
Every rational point is k G + T with T a torsion point, if G generates E(Q) modulo torsion (G is divided by 2, 3, 5, 7 as
long as PARI's ellisdivisible allows; a divisor of 11 or more is not excluded).  A perfect cuboid on the base gives such a point
at which X = t^2 is a nonzero square and X + a^2, X + b^2, X + d^2 are squares.  Modulo a prime p of good reduction these four
numbers must be squares or zero (or the point must be at infinity or at X = infinity, which is allowed, to be safe).  The values
of k G + T modulo p depend on k modulo the order of G, so each prime gives a table, and k survives if every prime allows it.
The torsion points (k = 0) include the trivial points X = 0 and X = infinity, which always survive.
Result for a base: the list of (k, T) with 0 < k <= N that survive; if it is empty, a perfect cuboid on that base needs |k| > N.
usage: open_sieve.py census.jsonl N [first [last]]       (bases number first..last of the open ones)"""
import sys, json, itertools, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
pari('''cubtab(ai, G, Ts, p, A, B, D, mode, useD) = { my(Ep = ellinit(ai, p), Gp = G*Mod(1, p), o = ellorder(Ep, Gp), out = vector(#Ts), P, v, x, X, ok);
  for (j = 1, #Ts, P = if(#Ts[j] == 1, [0], Ts[j]*Mod(1, p)); v = vectorsmall(o);
    for (k = 1, o,
      if (#P == 1, ok = 1, x = P[1]; if (mode && x == 0, ok = 1, X = if(mode, A*B*D/x, x);
        ok = kronecker(lift(X), p) >= 0 && kronecker(lift(X + A), p) >= 0 && kronecker(lift(X + B), p) >= 0 && (!useD || kronecker(lift(X + D), p) >= 0)));
      v[k] = ok; P = elladd(Ep, P, Gp));
    out[j] = v);
  [o, out] }''')
pari('divq(E, P, n) = { my(Q); if (ellisdivisible(E, P, n, &Q), Q, 0) }')
PR = [int(p) for p in pari.primes(400)][1:]
def torsion(E):
    T = pari.elltors(E); gens = list(T[2]); ords = [int(o) for o in T[1]]; pts = []
    for ks in itertools.product(*[range(o) for o in ords]):
        P = pari('[0]')
        for k, g in zip(ks, gens): P = pari.elladd(E, P, pari.ellmul(E, g, k))
        pts.append(P)
    assert len(pts) == int(T[0]); return pts
def sieve(rec, N, useD=1, only=None):
    m, n = rec['base']; a, b, d = rec['edges']; A, B, D = a*a, b*b, d*d
    curves = {'E_ab': ([0, A + B, 0, A*B, 0], 0), 'E_ad': ([0, A + D, 0, A*D, 0], 0), 'E_bd': ([0, B + D, 0, B*D, 0], 0),
              'E_2': ([0, B*D + A*D + A*B, 0, A*B*D*(A + B + D), (A*B*D)**2], 1)}
    out = dict(base=[m, n], edges=[a, b, d])
    for name in ('E_ab', 'E_ad', 'E_bd', 'E_2'):
        if rec['rank'][name][1] != 1 or (only and name != only): continue
        ai, mode = curves[name]; E = pari.ellinit(ai); r = pari.ellrank(E, 1)
        if not (int(r[0]) == 1 and int(r[1]) == 1 and len(r[3]) >= 1): continue
        G = r[3][0]; Ts = torsion(E); div = []
        for l in (2, 3, 5, 7):
            again = True
            while again:
                again = False
                for T in (Ts if l == 2 else [pari('[0]')]):
                    Q = pari('divq')(E, pari.elladd(E, G, T), l)
                    if Q != 0: G = Q; div.append(l); again = True; break
        h = float(pari.ellheight(E, G)); disc = int(pari.elldisc(E)) if hasattr(pari, "elldisc") else int(pari("(E)->E.disc")(E)); den = int(pari.denominator(pari.Vec(G)))
        alive = None; used = 0
        for p in PR:
            if disc % p == 0 or den % p == 0 or (A*B*D) % p == 0 or any(len(T) == 2 and int(pari.denominator(pari.Vec(T))) % p == 0 for T in Ts): continue
            o, tabs = pari('cubtab')(ai, G, Ts, p, A, B, D, mode, useD); o = int(o); used += 1
            tabs = [np.array([int(x) for x in v], dtype=bool) for v in tabs]
            if alive is None: alive = [np.nonzero(t[np.arange(1, N + 1) % o])[0] + 1 for t in tabs]
            else: alive = [k[t[k % o]] for k, t in zip(alive, tabs)]
            if sum(len(k) for k in alive) == 0 and used >= 12: break
            if not useD and used >= 150: break
        left = sorted((int(k), j) for j, ks in enumerate(alive) for k in ks)
        out.update(curve=name, height_of_G=round(h, 4), divided_by=div, torsion=len(Ts), primes_used=used, N=N, survivors=left[:20], status='no k up to N' if not left else 'SURVIVORS')
        return out
    out['status'] = 'no curve of rank exactly 1 with a point found'; return out
if __name__ == '__main__':
    if sys.argv[1] == '--control':
        # the Euler brick 44 x 117 x 240 on the base (11,2): without the condition on X + d^2 its point must survive on E_ab
        # (X = 240^2), and with the condition nothing may be left
        import math
        m, n = (int(v) for v in sys.argv[2].split(',')); rec = dict(base=[m, n], edges=[m*m - n*n, 2*m*n, m*m + n*n], rank={k: [1, 1] for k in ('E_ab', 'E_ad', 'E_bd', 'E_2')})
        for useD in (0, 1):
            o = sieve(rec, int(sys.argv[3]), useD, 'E_ab'); print('condition on X + d^2 %s:' % ('on' if useD else 'off'), json.dumps(o), flush=True)
        sys.exit(0)
    rows = [json.loads(l) for l in open(sys.argv[1])]; opn = [r for r in rows if r['status'] == 'open']; N = int(sys.argv[2])
    lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0; hi = int(sys.argv[4]) if len(sys.argv) > 4 else len(opn)
    for r in opn[lo:hi]: print(json.dumps(sieve(r, N)), flush=True)
