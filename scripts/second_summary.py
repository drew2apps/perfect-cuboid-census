#!/usr/bin/env python3
"""Counts for the output of second_step.py.  usage: second_summary.py second_step_open_100.txt"""
import sys, re
n = closed = opened = sol = tup = tupc = 0; ts = set(); firstopen = []
for l in open(sys.argv[1]):
    g = re.match(r'\((\d+),(\d+)\) base \((\d+), (\d+), (\d+)\): tuples (\d+), closed (\d+); rational t found: (.*); (closed|open|SOLUTION)', l)
    if not g: continue
    n += 1; tup += int(g.group(6)); tupc += int(g.group(7)); ts |= set(g.group(8).split())
    st = g.group(9); closed += st == 'closed'; opened += st == 'open'; sol += st == 'SOLUTION'
    if st == 'open' and len(firstopen) < 6: firstopen.append('(%s,%s)' % (g.group(1), g.group(2)))
print('bases: %d' % n)
print('closed by the second step (every tuple of the first step is closed, no solution): %d' % closed)
print('still open (some tuple has no curve with a rank equal to 0): %d; the first: %s' % (opened, ' '.join(firstopen)))
print('with a solution: %d' % sol)
print('tuples left by the first step: %d; closed: %d' % (tup, tupc))
print('rational t found on the curves: %s' % ' '.join(sorted(ts)))
