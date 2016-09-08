""" Collect pruned output and output a single BED file

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import operator
import sys

from ..algorithms import hashjoin
from ..formats import *

def _chromosome(c):
    with gzip.open('{}.txt.gz'.format(c), 'rt') as f, gzip.open('{}-pruned.gz'.format(c), 'rt') as g:
        original = parse(ucsc_bed_format, (line.split() for line in f))
        pruned = (line.split() for line in g)
        for o, p in hashjoin(original, pruned, key1=operator.itemgetter(3), key2=operator.itemgetter(1)):
            chr_, start, end, _, score = o
            name, *_ = p
            yield chr_, int(start), int(end), name, score

if __name__ == '__main__':
    for c in sorted([chromosome(x) for x in range(1, 23)]):
        for row in sorted(_chromosome(c)):
            print(*row, sep='\t')
