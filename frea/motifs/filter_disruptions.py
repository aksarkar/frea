"""Filter motif disruptions by information content

Usage: python filter_disruptions.py COUNTS PWMS

COUNTS contains space-separated records (master regulator name, cofactor name,
offset, count).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import operator
import math
import re
import sys

with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    types = [str, str, int, int]
    parsed = ([f(x) for f, x in zip(types, row)] for row in data)
    cofactors = {k: list(g) for k, g in itertools.groupby(parsed, key=operator.itemgetter(0))}
    lookup = collections.defaultdict(set)
    for _, pwm, offset, _ in itertools.chain.from_iterable(cofactors.values()):
        lookup[pwm].add(offset)

with open(sys.argv[2]) as f:
    is_pwm_row = lambda line: line[0] != '>'
    keep = {}
    for line in f:
        if line.startswith('>'):
            name = line.split()[0][1:]
            if name not in lookup:
                itertools.dropwhile(is_pwm_row, f)
                continue
            pwm = [[float(x) for x in line.split()[1:]] for line in itertools.takewhile(is_pwm_row, f)]
            for offset in lookup[name]:
                info = 2 - sum(-x * math.log(x, 2) for x in pwm[offset] if x > 0)
                if info > 1:
                    keep[(name, offset)] = info
for seed in cofactors:
    for _, name, offset, count in cofactors[seed]:
        if (name, offset) in keep:
            print(seed, name, offset, count, '{:.3f}'.format(keep[(name, offset)]))
