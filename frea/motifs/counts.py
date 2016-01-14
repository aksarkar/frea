"""Aggregate motif counts by group and region

Usage: python counts.py GROUPS

GROUPS is space-separated (motif, group) pairs, one per line. Expects
space-separated lines (chromosome, start, end, motif) on stdin. Prints (group,
count) pairs on stdout.

Position information is used to de-dupe per region; use 1 bp regions to count
SNPs.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import operator
import sys

with open(sys.argv[1]) as f:
    motif_to_group = dict(line.split() for line in f)

counter = collections.Counter()
data = (line.split() for line in sys.stdin)
for k, g in itertools.groupby(data, key=operator.itemgetter(0, 1, 2)):
    group = set(motif_to_group[motif] for _, _, _, motif in g)
    counter.update(group)
for k in counter:
    print(k, counter[k])
