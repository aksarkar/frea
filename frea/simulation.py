"""Simulate phenotypes with different genetic architectures.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import gzip
import random
import sys

import numpy
import numpy.random as R

from .formats import oxstats, oxstats_gen_to_dosage

def sample_uniform(data, p_causal=0.5, window_size=1e6, n_per_window=1):
    for _, g in itertools.groupby(data, key=lambda x: (x[0], int(int(x[2]) / window_size))):
        if random.random() < p_causal:
            w = list(g)
            for snp in random.sample(w, n_per_window):
                yield numpy.array(oxstats_gen_to_dosage(snp[5:]))

def simulate_null(samples, probs, pve=0.5, **kwargs):
    """Simulate phenotypes under the null (no enrichment)

    Under the null, variants are sampled uniformly at random, independent of
    annotation. We make a stronger assumption that causal variants are in
    approximate linkage equilibrium, and sample causal variants uniformly
    within 1MB windows.

    """
    y = numpy.zeros(len(samples) - 2)
    for dose in sample_uniform(probs, **kwargs):
        dose -= dose.mean()
        y += dose * R.normal()
    y += R.normal(scale=(1 / pve - 1) * y.std())
    y -= y.mean()
    y /= y.std()
    print(' '.join(samples[0]), 'pheno')
    print(' '.join(samples[1]), 'P')
    for i, y_i in enumerate(y):
        print(' '.join(samples[i + 2]), y_i)

if __name__ == '__main__':
    random.seed(sys.argv[1])
    if sys.argv[2] == '-':
        sys.argv[2] = None
    with oxstats(sys.argv[3], sys.argv[2]) as (samples, probs):
        simulate_null(samples, probs)
