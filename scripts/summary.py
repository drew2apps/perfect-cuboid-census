#!/usr/bin/env python3
"""Counts from base_census_100.jsonl.  usage: summary.py [file]"""
import sys, json, math, collections
rows = [json.loads(l) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'base_census_100.jsonl')]
M = max(r['base'][0] for r in rows)
want = [(m, n) for m in range(2, M + 1) for n in range(1, m) if math.gcd(m, n) == 1 and (m - n) % 2 == 1]
assert [tuple(r['base']) for r in rows] == want, 'the file does not list every base once, in order'
closed = [r for r in rows if r['status'] == 'closed']; opn = [r for r in rows if r['status'] == 'open']
print('primitive bases (m, n) with m <= %d: %d' % (M, len(rows)))
print('closed by a curve of rank 0 (no perfect cuboid on the base): %d; open: %d' % (len(closed), len(opn)))
for name in ('E_ab', 'E_ad', 'E_bd', 'E_2'):
    print('  %s has rank 0 and closes: %d; it is the only one that does: %d' % (name, sum(name in r['closed_by'] for r in rows), sum(r['closed_by'] == [name] for r in rows)))
print('bases where a curve of rank 0 has a torsion point with X a nonzero square: %d' % sum('torsion_with_square_X' in r for r in rows))
undet = sum(any(lo != hi for lo, hi in r['rank'].values()) for r in rows)
print('bases with a rank that PARI does not determine (lower bound below upper bound): %d, of them open: %d' % (undet, sum(any(lo != hi for lo, hi in r['rank'].values()) for r in opn)))
c = collections.Counter(tuple(sorted(hi for lo, hi in r['rank'].values())) for r in opn)
print('open bases by the upper bounds of the four ranks:', ', '.join('%s: %d' % (k, v) for k, v in sorted(c.items())))
for b in (20, 40, 60, 80, 100):
    if b <= M: print('  m <= %d: %d of %d closed' % (b, sum(r['base'][0] <= b for r in closed), sum(r['base'][0] <= b for r in rows)))
print('the first open bases:', ' '.join('(%d,%d)' % tuple(r['base']) for r in opn[:12]))
