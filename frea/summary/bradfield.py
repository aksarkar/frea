import itertools
import gzip
import math
import operator
import os
import os.path
import re
import sys

import scipy.stats

from ..algorithms import join

def bradfield(file):
    types = [str, int, int, float]
    next(file)
    data = (line.split('\t') for line in file)
    for row in data:
        if not row[0]:
            row[0] = '.'
        yield [f(x) for f, x in zip(types, row)]

def snptest(file):
    with gzip.open(file, 'rt') as f:
        next(f)
        data = (line.split() for line in f)
        for row in data:
            yield row[1], int(row[3]), row[4], row[5], float(row[34])

def wtccc(chromosome, snptest_dir):
    patt = re.compile(r'{}\..*\.txt\.gz'.format(chromosome))
    infiles = sorted(os.path.join(root, f) for root, _, files in os.walk(snptest_dir)
                     for f in files if patt.match(f))
    for f in infiles:
        for row in snptest(f):
            yield row

def zscores(chromosome, snps, snptest_dir):
    for pair in join(snps, wtccc(chromosome, snptest_dir), key1=operator.itemgetter(2, 0), key2=operator.itemgetter(1, 0)):
        name, _, pos, p = pair[0]
        _, _, a0, a1, odds = pair[1]
        z = math.sqrt(scipy.stats.chi2(1).isf(p))
        if odds < 1:
            z = -z
        yield name, pos, a0, a1, z

if __name__ == '__main__':
    summary_file = sys.argv[1]
    snptest_dir = sys.argv[2]
    with gzip.open(summary_file, 'rt') as f:
        summary = sorted(bradfield(f), key=operator.itemgetter(1, 2))
        for k, g in itertools.groupby(summary, key=operator.itemgetter(1)):
            with open('{}.txt'.format(k), 'w') as outfile:
                print('id', 'pos', 'ref', 'alt', 'z', file=outfile)
                for row in zscores(k, g, snptest_dir):
                    print(*row, file=outfile)
