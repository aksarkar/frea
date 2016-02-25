"""Compute p-values for functional enrichment against matched sets

Usage: python test.py TABLE THRESH NTRIALS

Expects BED file of SNPs with score = negative log p-value on stdin. TABLE is
gzipped space-separated (rsid, key1[, key2, ...]) tuples.

This implementation resamples overlaps and counts how many resampled SNPs
exceed THRESH over resampled sets.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import contextlib
import gzip
import heapq
import itertools
import operator
import pdb
import random
import sys

import scipy.stats

from .algorithms import merge, moments
from .formats import get_pouyak_name, chromosome

I = itertools.chain.from_iterable

def extract_maf(n):
    """Extract MAF from OXSTATS legend file.

    This function is an altnerative entry point, intended to be called using
    python -c. (This can't be installed as an entry point with setuptools
    because the argument is required to uniformly format SNP names.)

    Expects legend file on stdin.

    """
    data = (line.split() for line in sys.stdin)
    next(data)
    for row in data:
        if row[4] == 'SNP' and float(row[-1]) > 0.01:
            row[1] = int(row[1])
            print(get_pouyak_name(chromosome(n), *row[:4]), row[-1])

def build_table():
    """Build lookup table for enrichment test.

    This is installed as an entry point (frea-build-table).

    """
    with contextlib.ExitStack() as stack:
        files = [stack.enter_context(gzip.open(f, 'rt')) for f in sys.argv[1:]]
        iters = [(line.split() for line in f) for f in files]
        for it in iters:
            next(it)
        key=operator.itemgetter(0)
        table = []
        for k, g in itertools.groupby(merge(*iters, key=key), key=key):
            props = list(g)
            if len(props) == len(iters):
                table.append([k] + [x[1] for x in props])
        for snp in sorted(table, key=key):
            print(*snp)

def build_bins(snps, table):
    """Load lookup table and partition it according to annotated SNPs.

    Return dictionaries from SNP properties to SNPs.

    """
    overlap_bins = collections.defaultdict(list)
    nonoverlap_bins = collections.defaultdict(list)
    for rsid, *key in table:
        k = tuple(key)
        if snps[rsid][1]:
            overlap_bins[k].append(rsid)
        else:
            nonoverlap_bins[k].append(rsid)
    return overlap_bins, nonoverlap_bins

def match(overlap_counts, nonoverlap_bins):
    """Resample overlap_counts SNPs from nonoverlap_bins"""
    return I((random.choice(nonoverlap_bins[k]) for _ in range(v)) for k, v in
             overlap_counts.items() if k in nonoverlap_bins)

def num_top_snps(rsids, thresh):
    """Return the number of SNPs with p < thresh"""
    result = 0
    for s in rsids:
        if snps[s][0] > thresh:
            result += 1
    return result

def permutation_test(overlap_bins, nonoverlap_bins, thresh, ntrials):
    """Perform permutation test for enrichment.

    This implementation tests the count of annotated SNPs with p < thresh
    against the count of resampled SNPs with p < thresh. Resampled SNPs are
    taken from outside of the annotation of interest and matched on annotated
    SNP properties.

    Returns a tuple for long-form output:
    - threshold
    - original count
    - mean null count
    - variance null count
    - exact p value
    - Anderon-Darling statistic for normality of the null distribution
    - Anderson-Darlign .05 critical value

    """
    X = num_top_snps(I(overlap_bins.values()), thresh)
    if X == 0:
        return 0, 0, 0, 0, 1, 0, 0
    overlap_counts = {k: len(overlap_bins[k]) for k in overlap_bins}
    Y = [num_top_snps(match(overlap_counts, nonoverlap_bins), thresh) for _ in range(ntrials)]
    mean, variance = moments(Y)
    anderson, critical_values, _ = scipy.stats.anderson(Y)
    exact_p = (1 + len([y for y in Y if y >= X])) / (1 + ntrials)
    return thresh, X, mean, variance, exact_p, anderson, critical_values[2]

if __name__ == '__main__':
    random.seed(0)
    snps_raw = (line.split() for line in sys.stdin)
    snps = {row[3]: (float(row[4]), int(row[5])) for row in snps_raw}

    with gzip.open(sys.argv[1], 'rt') as f:
        data = (line.split() for line in f)
        filter_ = (d for d in data if d[0] in snps)
        overlap_bins, nonoverlap_bins = build_bins(snps, filter_)

    thresh = float(sys.argv[2])
    ntrials = int(sys.argv[3])
    result = permutation_test(overlap_bins, nonoverlap_bins, thresh, ntrials)
    print(' '.join('{:.3g}'.format(x) for x in result))
