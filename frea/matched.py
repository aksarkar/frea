"""Compute p-values for functional enrichment against matched sets

Usage: python test.py TABLE THRESH NTRIALS

Expects BED file of SNPs with score = negative log p-value on stdin. TABLE is
gzipped space-separated (rsid, key1[, key2, ...]) tuples.

This implementation resamples overlaps and counts how many resampled SNPs
exceed THRESH over resampled sets.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import gzip
import itertools
import random
import sys

import scipy.stats

I = itertools.chain.from_iterable

def build_bins(snps, table):
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
    return I((random.choice(nonoverlap_bins[k]) for _ in range(v)) for k, v in
             overlap_counts.items() if k in nonoverlap_bins)

def num_top_snps(rsids, thresh):
    result = 0
    for s in rsids:
        if snps[s][0] > thresh:
            result += 1
    return result

def moments(xs):
    running_mean = xs.pop(0)
    running_variance = 0
    for i, x in enumerate(xs):
        new_mean = running_mean + (x - running_mean) / (i + 2)
        running_variance += (x - running_mean) * (x - new_mean)
        running_mean = new_mean
    running_variance /= ntrials
    return running_mean, running_variance

def permutation_test(overlap_bins, nonoverlap_bins, thresh, ntrials):
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
