""" Sample causal variants to feed into simulation

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import contextlib
import itertools
import operator
import gzip
import random
import signal
import sys

def _intersect(snps, regions):
    snps = iter(snps)
    regions = iter(regions)
    a = next(snps)
    b = next(regions)    
    while True:
        if a[3] > b[2]:
            b = next(regions)
        else:
            if a[3] >= b[1]:
                yield a
            a = next(snps)

def intersect(snps, regions):
    regions_by_chromosome = {k: list(g) for k, g in itertools.groupby(regions, key=operator.itemgetter(0))}
    for k, g in itertools.groupby(snps, key=operator.itemgetter(0)):
        for snp in _intersect(g, regions_by_chromosome['chr{}'.format(k)]):
            yield snp

def sample_annotated(snps, regions, **kwargs):
    annotated = intersect(snps, regions)
    for snp in sample_uniform(annotated):
        yield snp

def sample_uniform(data, n_causal=1000, window_size=1e6, n_per_window=1):
    windows = [list(g) for _, g in itertools.groupby(data, key=lambda x: (x[0], int(x[3] / window_size)))]
    n = 0
    for w in random.sample(windows, n_causal):
        for _, snp, *_ in random.sample(w, n_per_window):
            yield snp

@contextlib.contextmanager
def bim(filename):
    with open(filename) as f:
        types = [int, str, int, int]
        data = (line.split() for line in f)
        try:
            yield ([g(x) for g, x in zip(types, row)] for row in data)
        except Exception as e:
            raise e

@contextlib.contextmanager
def ucsc_bed(filename):
    with gzip.open(filename, 'rt') as f:
        types = [str, int, int, str, float]
        data = (line.split() for line in f)
        try:
            yield ([g(x) for g, x in zip(types, row)] for row in data)
        except Exception as e:
            raise e

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    random.seed(int(sys.argv[1]))
    with bim(sys.argv[2]) as snps, ucsc_bed(sys.argv[3]) as regions:
        for x in sample_annotated(snps, regions):
            print(x)

