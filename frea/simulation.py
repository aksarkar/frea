"""Simulate phenotypes with different genetic architectures.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import gzip
import os.path
import random
import sys

import numpy
import numpy.random as R

from .formats import *

def windows(key, ref_dir='.', p_causal=0.5, window_size=1e6, n_per_window=1, **kwargs):
    J = os.path.join
    for chr_ in [str(i) for i in range(1, 23)]:
        sample_file = J(ref_dir, 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample')
        legend_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(chr_))
        haps_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'.format(chr_))
        with oxstats_haplotypes(sample_file, legend_file, haps_file, **kwargs) as data:
            for _, g in itertools.groupby(data, key=key):
                yield g

def sample_gaussian(n, pve=0.5, ref_dir='.', p_causal=0.5, window_size=1e6,
                    n_per_window=1, **kwargs):
    """Generate Gaussian phenotypes.

    This implementation recombines reference haplotypes to generate arbitrary
    size samples. We make two passes through the reference data:

    1. Sample causal variants, effect sizes, and recombinations (at fixed
       intervals); reconstruct locally to accumulate liabilities and genetic
       variance

    2. Add Gaussian noise; reconstruct the sample genotypes and compute
       marginal association statistics

    n - target number of samples
    pve - target proportion of variance explained
    ref_dir - base directory of 1KG reference haplotypes
    kwargs - keyword arguments to oxstats_haplotypes

    """
    key = lambda x: int(int(x[0][1]) / window_size)
    y = numpy.zeros(n)
    mosaic = []
    ## Pass 1
    for window in windows(key, ref_dir=ref_dir, **kwargs):
        window = numpy.array([h for _, h in window])
        p, k = window.shape
        mosaic.append(R.randint(0, k, size=2 * n))
        if R.rand() < p_causal:
            causal = R.choice(numpy.arange(p), size=n_per_window, replace=False)
            effects = R.normal(size=n_per_window)
            w = window.T[mosaic[-1], causal]  # 2n x n_per_window
            x = w[::2] + w[1::2]
            y += x.reshape(-1, n_per_window).dot(effects)
    y += R.normal(scale=(1 / pve - 1) * y.std(), size=n)
    y -= y.mean()
    y -= y.std()
    ## Pass 2
    numpy.seterr(all='warn')
    for ancestors, window in zip(mosaic, windows(key, ref_dir=ref_dir, **kwargs)):
        window = list(window)
        w = numpy.array([h for _, h in window]).T[ancestors]  # 2n x p
        x = (w[::2] + w[1::2]).astype('float32')
        x -= x.mean(axis=0)
        var = numpy.diag(x.T.dot(x)) + 1e-8  # Needed for monomorphic SNPs
        b = y.T.dot(x).T / var
        rss = numpy.inner(y - x.dot(b), y - x.dot(b))
        se = rss / var
        stat = numpy.square(b / se)
        logp = -scipy.stats.chi2(1).logsf(stat)
        for (l, _), p in zip(window, logp):
            print(' '.join(l), p)

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
