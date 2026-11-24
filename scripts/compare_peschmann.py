#!/usr/bin/env python3
"""Compares the bases of this census with the fibers of R. Peschmann, "A torsion-intersection proof of perfect-cuboid
nonexistence on 1,072 explicit master-tuple fibers" (arXiv:2604.28072).  His fiber (m, n) fixes the ratio of two edges,
2 m n : m^2 - n^2, which is the base (m, n) here; his range, coprime m > n of different parity with m <= 100, is the same
2,040 pairs.  The list is his file paper3/data/proven_fibers.csv (github.com/renpe/euler-brick-obstructions, MIT licence;
commit 3585937 of 4 May 2026 was used); it is not copied into this repository.
usage: compare_peschmann.py path/to/proven_fibers.csv"""
import sys, csv, json, re
from collections import Counter
rows = {(int(r['m']), int(r['n'])): r for r in csv.DictReader(open(sys.argv[1]))}
cen = {tuple(json.loads(l)['base']): json.loads(l) for l in open('base_census_100.jsonl')}
st = {}
for l in open('second_step_open_100.txt'):
    g = re.match(r'\((\d+),(\d+)\).*; (closed|open|SOLUTION)', l); st[(int(g.group(1)), int(g.group(2)))] = g.group(3)
status = {k: 'census' if r['status'].startswith('closed') else ('two-step' if st.get(k) == 'closed' else 'open') for k, r in cen.items()}
print('his fibers: %d; bases here: %d; his fibers that are bases here: %d' % (len(rows), len(cen), len(set(rows) & set(cen))))
c = Counter((k in rows, status[k]) for k in cen)
print('closed by the census here (rank 0 of one of four curves, PARI): %d, of which in his list: %d' % (c[(True, 'census')] + c[(False, 'census')], c[(True, 'census')]))
print('closed by the two-step descent here: %d, of which in his list: %d' % (c[(True, 'two-step')] + c[(False, 'two-step')], c[(True, 'two-step')]))
print('open here: %d, of which in his list: %d: %s' % (c[(True, 'open')] + c[(False, 'open')], c[(True, 'open')], ' '.join('(%d,%d)' % k for k in sorted(k for k in rows if status[k] == 'open'))))
miss = Counter(rows[k]['rank_resolution'] for k in rows if status[k] != 'census')
print('his fibers not closed by the census here, by the way he found the rank: %s' % dict(miss))
print('closed by him or here: %d; open in both: %d' % (sum(1 for k in cen if k in rows or status[k] != 'open'), sum(1 for k in cen if k not in rows and status[k] == 'open')))
