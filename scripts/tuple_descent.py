#!/usr/bin/env python3
"""A descent with one rational class for each of several quadratic forms that are norms from Z[i] (EXPERIMENT).

Given linear forms gamma_j(u, v) = alpha_j u + beta_j v over Z[i] and the system  N(gamma_j(u, v)) = square  for all j,
with (u, v) coprime integers.  Each gamma_j(u, v) then has square norm, so gamma_j = d_j z_j^2 with d_j a squarefree
positive integer (an element of Q(i) of square norm is a rational number times a square; -1 is a square): the tuple
(d_1, ..., d_J) is the class of the solution.  An odd prime p divides d_j only if a prime above p divides gamma_j and its
conjugate at (u, v), hence divides the determinant of the two; so d_j is a product of primes of a finite set P_j, given by
the caller.

Local test.  A tuple is possible at p only if some (u : v) over Q_p makes every d_j gamma_j(u, v) a square in Q_p(i)
(in both factors of Q_p x Q_p if p = 1 mod 4).  P^1(Z_p) is covered by discs t0 + p^k Z_p in the charts (t, 1) and (1, t),
p | t.  On a disc gamma_j(t) = gamma_j(t0) + (a multiple of p^k), so the square class of gamma_j is constant on the disc as
soon as its valuation at t0 is below k (for p = 2: valuation at 1 + i at most 2k - 6).  A disc is dropped when a constant
class is not the class of a rational number, kept when all are constant, and cut into p discs otherwise, down to depth K,
where a class that is not yet constant counts as any class.  The coprime integers (u, v) of a solution are (t, 1) or (1, t)
times a p-adic unit, and the forms are linear, so the allowed set is closed under multiplying all classes by the class of a
unit of Z_p.  So the set of tuples allowed at p can be too large, never too small (a first version forgot the unit; the
test with random points, `selftest`, found it).  The tuples that pass at every prime are found by joining the forms one after the other.

Systems:  magic m,n   the four forms F0, F8, GL, GU of the family (classes d0, d8, dL, dU; the class of LAST_PAIR.md,
                      section 11, is their product);
          cuboid m,n  the three forms t^2 + a^2, t^2 + b^2, t^2 + c^2 for the Pythagorean triple (a, b, c) =
                      (m^2 - n^2, 2 m n, m^2 + n^2): a perfect cuboid with this face.
usage: tuple_descent.py magic|cuboid m,n [largest good prime, default 150]
       tuple_descent.py magic|cuboid m,n selftest | neartest | oracle
       tuple_descent.py batch magic|cuboid pairs.txt jobs [first [last]]     one line per ratio; the file has one ratio m,n at the
                                                                            start of each line, or ratios separated by spaces on one line"""
import sys
from math import gcd, prod
from itertools import product
from sympy import factorint, primerange

def vp(x, p):
    if x == 0: return None
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
def leg(a, p): return 0 if pow(a % p, (p - 1)//2, p) == 1 else 1
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def sqf(x): return prod(q for q, e in factorint(abs(x)).items() if e % 2)
SQ8 = sorted(set(((x*x - y*y) % 8, (2*x*y) % 8) for x in range(8) for y in range(8) if (x + y) % 2))
def key2(z):
    """Valuation at 1 + i of a nonzero Gaussian integer, and its class in Q_2(i)* modulo squares."""
    x, y = z; w = 0
    while (x + y) % 2 == 0: x, y = (x + y)//2, (y - x)//2; w += 1
    return w, (w % 2, min(((x*a - y*b) % 8, (x*b + y*a) % 8) for a, b in SQ8))
D2 = {key2((d, 0))[1]: None for d in (1, 2, 3, 6)}            # the classes of rational numbers in Q_2(i)

class Local:
    def __init__(self, forms, p, K):
        self.forms, self.p, self.K = forms, p, K
        self.kind = 2 if p == 2 else (1 if p % 4 == 1 else 3)
        if self.kind == 1:
            r = next(x for x in range(2, p) if (x*x + 1) % p == 0); M = p**(K + 4)
            for _ in range(K + 6): r = (r - (r*r + 1)*pow(2*r, -1, M)) % M
            self.iota, self.M = r, M
    def cls(self, z, k):
        """(dead, determined, label): the class of gamma = z on a disc of radius p^k, as a label of dlabel; dead if the class is
        constant and is not the class of a rational number.  When it is not determined, label is None or a hint (p = 1 mod 4)."""
        p = self.p
        if z == (0, 0): return False, False, None
        if self.kind == 2:
            w, key = key2(z)
            if w + 6 > 2*k: return False, False, None
            return (key not in D2), True, key
        if self.kind == 3:
            v = min(x for x in (vp(z[0], p), vp(z[1], p)) if x is not None)
            if v >= k: return False, False, None
            a, b = (z[0]//p**v) % p, (z[1]//p**v) % p
            return (leg(a*a + b*b, p) == 1), True, (v % 2, 0)       # a unit of Z_p[i] is a square iff its norm is a square in F_p
        out = []
        for s in (1, -1):
            w = (z[0] + s*self.iota*z[1]) % self.M; v = vp(w, p)
            if v is None: return False, False, None                 # beyond the precision of iota
            out.append(None if v >= k else (v % 2, leg(w//p**v, p)))
        # The same rational d has to work at both primes above p.  If only one component is constant on the disc, the class is
        # not determined yet, but a rational class on the disc can only be the class of that component: it is returned as a
        # hint, which allowed() uses at the largest depth in place of "every class" (LAST_PAIR.md, section 14, R4).
        if out[0] is None or out[1] is None: return False, False, (out[1] if out[0] is None else out[0])
        return (out[0] != out[1]), True, out[0]
    def dlabel(self, d):
        if self.kind == 2: return key2((d, 0))[1]
        v = vp(d, self.p)
        return (v % 2, leg(d//self.p**v, self.p) if self.kind == 1 else 0)
    def labels(self):
        if self.kind == 2: return list(D2)
        return [(0, 0), (1, 0)] if self.kind == 3 else [(0, 0), (0, 1), (1, 0), (1, 1)]
    def allowed(self):
        """A set of tuples of labels containing the labels of every tuple of classes that is possible at p."""
        p, K, forms = self.p, self.K, self.forms; out = set(); every = self.labels()
        stack = [(c, t0, 1) for c in (0, 1) for t0 in range(p) if c == 0 or t0 == 0]
        while stack:
            c, t0, k = stack.pop(); u, v = (t0, 1) if c == 0 else (1, t0)
            res = [self.cls((al[0]*u - 0 + be[0]*v, al[1]*u + be[1]*v), k) for al, be in forms]
            if any(r[0] for r in res): continue
            if all(r[1] for r in res): out.add(tuple(r[2] for r in res)); continue
            if k < K: stack += [(c, t0 + j*p**k, k + 1) for j in range(p)]; continue
            out |= set(product(*[[r[2]] if (r[1] or r[2] is not None) else every for r in res]))     # r[2] without r[1]: the hint of cls
        # The integers (u, v) are (t, 1) or (1, t) times a p-adic unit, and the forms are linear: all classes are multiplied
        # by the class of that unit.  Units of Z_p are squares in Q_p(i) when p = 3 mod 4; otherwise there is one more class.
        if self.kind == 1: out |= set(tuple((l[0], 1 - l[1]) for l in t) for t in out)
        if self.kind == 2: out |= set(tuple(self.times3(l) for l in t) for t in out)
        return out
    def times3(self, key):
        x, y = key[1]
        return (key[0], min(((3*x*a - 3*y*b) % 8, (3*x*b + 3*y*a) % 8) for a, b in SQ8))

def descend(forms, supports, known, pmax=150, K=None, extra_primes=()):
    """The tuples (d_1, ..., d_J) that pass at every prime tested; known: tuples that must pass (checked)."""
    J = len(forms); divs = []
    for P in supports:
        ds = [1]
        for q in P: ds += [x*q for x in ds]
        divs.append(ds)
    tests = sorted(set(q for P in supports for q in P) | set(primerange(2, pmax + 1)) | set(extra_primes))
    loc = []
    for p in tests:
        L = Local(forms, p, K or (12 if p == 2 else (6 if p < 50 else 3))); al = L.allowed()
        for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('a known point is rejected', p, kt)
        if len(al) < len(L.labels())**J: loc.append((p, L, al))
    S = [()]
    for j in range(J):                                      # join: the prefixes (d_1, ..., d_j) allowed at every prime
        new = []
        projs = [(L, set(a[:j + 1] for a in al), {d: L.dlabel(d) for d in divs[j]}) for p, L, al in loc]
        cache = [dict() for _ in loc]
        for s in S:
            keys = []
            for idx, (L, proj, lab) in enumerate(projs):
                if s not in cache[idx]: cache[idx][s] = tuple(L.dlabel(d) for d in s)
                keys.append(cache[idx][s])
            for d in divs[j]:
                if all(keys[idx] + (lab[d],) in proj for idx, (L, proj, lab) in enumerate(projs)): new.append(s + (d,))
        S = new
    return sorted(S), [p for p, L, al in loc], [2**len(P) for P in supports]

def magic(m, n):
    al = lambda z: gmul((1, 1), z)
    forms = [((0, n), (m, 0)), ((0, -m), (n, 0)), (al((m, -n)), al((-n, m))), (al((m, n)), al((n, m)))]   # F0, F8, GL, GU
    P0 = sorted(set([2]) | set(factorint(m)) | set(factorint(n))); PL = sorted(set([2]) | set(factorint(abs(m - n))) | set(factorint(m + n)))
    S = m*m + n*n; D = m*m - n*n; known = []
    if int(S**0.5 + 0.5)**2 == S:                           # t = 1 and t = -1
        for u, v in ((1, 1), (-1, 1)): known.append(tuple(cls_of((a[0]*u + b[0]*v, a[1]*u + b[1]*v)) for a, b in forms))
    if int((2*S)**0.5 + 0.5)**2 == 2*S:                     # t = 0 and t = infinity
        for u, v in ((0, 1), (1, 0)): known.append(tuple(cls_of((a[0]*u + b[0]*v, a[1]*u + b[1]*v)) for a, b in forms))
    return forms, [P0, P0, PL, PL], sorted(set(known))

def cuboid(m, n):
    a, b, c = m*m - n*n, 2*m*n, m*m + n*n
    forms = [((1, 0), (0, x)) for x in (a, b, c)]           # u + i x v, of norm u^2 + x^2 v^2
    sup = [sorted(set([2]) | set(factorint(x))) for x in (a, b, c)]
    known = sorted(set(tuple(cls_of((al[0]*u + be[0]*v, al[1]*u + be[1]*v)) for al, be in forms) for u, v in ((0, 1), (1, 0))))
    return forms, sup, known

def cls_of(z):
    """The squarefree positive integer d with z = d (square) in Q(i), for a Gaussian integer z of square norm."""
    x, y = z; g = gcd(x, y); x, y = x//g, y//g                # z = g (x + i y), x + i y primitive of square norm
    # a primitive Gaussian integer of square norm is a unit times a square, or (1 + i) times that with norm 2 (square): only units
    w = 0
    while (x + y) % 2 == 0: x, y = (x + y)//2, (y - x)//2; w += 1
    assert w % 2 == 0, 'the norm is not a square'
    unit = 1
    for cand, val in (((1, 0), 1), ((-1, 0), 1), ((0, 1), 2), ((0, -1), 2)):    # x + i y = unit * square: test by a square root
        zz = gmul((x, y), (cand[0], -cand[1]))                 # divide by the unit
        r = complex(zz[0], zz[1])**0.5; rr = (round(r.real), round(r.imag))
        if gmul(rr, rr) == zz: unit = val; break
    else: raise ValueError('not a unit times a square: %s' % (z,))
    return sqf(g*unit)                                      # (1 + i)^w = (2 i)^(w/2), and 2 i is a square times -1 times ... a square class: 2 and i are in the same class

def _job(a): return line(a[0], *(int(x) for x in a[1].split(',')))

def selftest(forms, primes, samples=20000, seed=1):
    """Random coprime (u, v): whenever no form is dead at p (decided from the exact values), the tuple of labels has to be in
    the allowed set.  Returns (points tested, failures)."""
    import random
    rnd = random.Random(seed); hits = fails = 0
    for p in primes:
        L = Local(forms, p, 12 if p == 2 else 6); al = L.allowed()
        for _ in range(samples):
            u, v = rnd.randrange(-10**6, 10**6), rnd.randrange(1, 10**6)
            if gcd(u, v) != 1: continue
            res = [L.cls((a[0]*u + b[0]*v, a[1]*u + b[1]*v), 10**6) for a, b in forms]
            if any(r[0] or not r[1] for r in res): continue
            hits += 1; fails += tuple(r[2] for r in res) not in al
    return hits, fails

def neartest(forms, primes, depth=6, samples=60, seed=2):
    """Points close to the zeros of the forms and to infinity, where the classes change: for p = 1 mod 4 each linear form
    has a zero t* in Q_p (in both embeddings); the points t* + p^j r, j = 1..depth, and (1, p^j r), are tested like the
    random points of selftest.  Returns (points tested, failures)."""
    import random
    rnd = random.Random(seed); hits = fails = 0
    for p in primes:
        L = Local(forms, p, 12 if p == 2 else 6); al = L.allowed(); starts = [(0, 0), (1, 0)]
        if L.kind == 1:
            M = L.M
            for al_, be in forms:
                for sgn in (1, -1):
                    a = (al_[0] + sgn*L.iota*al_[1]) % M; b = (be[0] + sgn*L.iota*be[1]) % M
                    if a % p: starts.append((0, (-b*pow(a, -1, M)) % M))
                    elif b % p: starts.append((1, (-a*pow(b, -1, M)) % M))
        for c, t0 in starts:
            for j in range(1, depth + 1):
                for _ in range(samples):
                    t = t0 + p**j*rnd.randrange(1, 10**6)
                    u, v = (t, 1) if c == 0 else (1, t)
                    if c == 1 and t % p: continue
                    res = [L.cls((a[0]*u + b[0]*v, a[1]*u + b[1]*v), 10**6) for a, b in forms]
                    if any(r[0] or not r[1] for r in res): continue
                    hits += 1; fails += tuple(r[2] for r in res) not in al
    return hits, fails

def oracle(forms, primes=((2, 6), (3, 4), (5, 3), (7, 3), (13, 2)), samples=3000, seed=3):
    """A second implementation of the local classes, by brute force.  For p^N small, the set of all squares w^2 in
    Z[i]/(p^N) is listed.  For random (u, v) and each form, the rational classes d (four representatives) with
    d gamma(u, v) a square modulo p^N are found from that list, when the valuation of gamma is small enough for the
    congruence to decide (Hensel); the program's own label has to name exactly these d.  Returns (values tested, failures)."""
    import random
    rnd = random.Random(seed); tested = fails = 0
    for p, N in primes:
        M = p**N; SQ = set(((x*x - y*y) % M, (2*x*y) % M) for x in range(M) for y in range(M))
        L = Local(forms, p, 12 if p == 2 else 6)
        nonres = next(r for r in range(2, p) if pow(r, (p - 1)//2, p) == p - 1) if p > 2 else None
        reps = [1, 2, 3, 6] if p == 2 else [1, nonres, p, p*nonres]
        for _ in range(samples):
            u, v = rnd.randrange(-10**6, 10**6), rnd.randrange(1, 10**6)
            if gcd(u, v) != 1: continue
            for a, b in forms:
                z = (a[0]*u + b[0]*v, a[1]*u + b[1]*v); nz = z[0]*z[0] + z[1]*z[1]
                if nz == 0 or vp(nz, p) > (N - 4 if p == 2 else N - 2): continue      # the congruence modulo p^N decides
                good = [d for d in reps if ((d*z[0]) % M, (d*z[1]) % M) in SQ]
                dead, det, lab = L.cls(z, 10**6); tested += 1
                mine = [] if dead else [d for d in reps if L.dlabel(d) == lab]
                fails += sorted(good) != sorted(mine)
    return tested, fails

def line(kind, m, n, pmax=150):
    forms, sup, known = (magic if kind == 'magic' else cuboid)(m, n)
    S, used, sizes = descend(forms, sup, known, pmax)
    if S != known and pmax < 1500: S, used, sizes = descend(forms, sup, known, 1500); pmax = 1500
    return '(%d,%d): %d tuples; good primes up to %d; known: %s; left: %s; %s' % (m, n, prod(sizes), pmax, known, S if len(S) <= 8 else '%d tuples' % len(S),
        'nothing is left' if not S else ('exactly the known points' if S == known else 'MORE THAN THE KNOWN POINTS'))

if __name__ == '__main__':
    if sys.argv[1] == 'batch':                                # batch magic|cuboid pairs.txt jobs [first [last]]
        from multiprocessing import Pool
        from functools import partial
        kind = sys.argv[2]; words = open(sys.argv[3]).read().split('\n')
        pairs = [w for l in words for w in (l.split() if all(',' in x for x in l.split()) else l.split()[:1]) if ',' in w and w.replace(',', '').isdigit()]
        lo = int(sys.argv[5]) if len(sys.argv) > 5 else 0; hi = int(sys.argv[6]) if len(sys.argv) > 6 else len(pairs)
        with Pool(int(sys.argv[4])) as pool:
            for out in pool.imap(_job, [(kind, a) for a in pairs[lo:hi]], chunksize=8): print(out, flush=True)
        sys.exit(0)
    if len(sys.argv) > 3 and sys.argv[3] in ('neartest', 'oracle'):
        m, n = (int(x) for x in sys.argv[2].split(',')); forms, sup, known = (magic if sys.argv[1] == 'magic' else cuboid)(m, n)
        if sys.argv[3] == 'neartest':
            pr = [2, 5, 13, 17, 29, 37, 41, 53, 61]; hits, fails = neartest(forms, pr)
            print('%s (%d,%d): points near the zeros of the forms and near infinity, at the primes %s: %d tested; outside the allowed set: %d' % (sys.argv[1], m, n, pr, hits, fails))
        else:
            tested, fails = oracle(forms)
            print('%s (%d,%d): local classes against a list of all squares modulo 2^6, 3^4, 5^3, 7^3, 13^2: %d values of forms tested; disagreements: %d' % (sys.argv[1], m, n, tested, fails))
        sys.exit(0)
    if len(sys.argv) > 3 and sys.argv[3] == 'selftest':
        m, n = (int(x) for x in sys.argv[2].split(',')); forms, sup, known = (magic if sys.argv[1] == 'magic' else cuboid)(m, n)
        pr = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 53]; hits, fails = selftest(forms, pr)
        print('%s (%d,%d): random points at the primes %s: %d at which every form has the class of a rational number; outside the allowed set: %d' % (sys.argv[1], m, n, pr, hits, fails)); sys.exit(0)
    kind = sys.argv[1]; m, n = (int(x) for x in sys.argv[2].split(',')); pmax = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    forms, sup, known = (magic if kind == 'magic' else cuboid)(m, n)
    S, used, sizes = descend(forms, sup, known, pmax)
    print('%s (%d,%d): classes for each form: %s; tuples: %d' % (kind, m, n, sizes, prod(sizes)))
    print('classes of the known points: %s' % (known or 'none'))
    print('tuples left: %d: %s' % (len(S), S[:20]))
    if kind == 'magic': print('products of the tuples left (the class d of section 11): %s' % sorted(set(sqf(prod(s)) for s in S)))
    print('nothing is left: no rational point' if not S else ('exactly the known points' if S == known else 'more than the known points'))
