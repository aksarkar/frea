"""Simulate phenotypes with different genetic architectures.

Usage: python -m frea.simulation SEED N LEGEND HAPS SAMPLE

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import argparse
import collections
import itertools
import gzip
import pickle
import random
import sys

import numpy
import numpy.random as R

from .formats import *

def blocks(hotspot_file, **kwargs):
    """Yield LD blocks defined by recombination hotspots

    Blocks are pairs of iterables (instances of itertools.takewhile), allowing
    sequential access to OXSTATS reference haplotypes, and recombination rates
    at the block boundary (next recombination). Hotspots are taken as the
    center of contiguous 10kb regions with significantly different
    recombination rate from flanking regions. Recombination rates are taken as
    the mean of standardized recombination rates estimated in each 10kb window,
    in units of percentage of haplotypes broken by a recombination event.

    hotspot_file - deCODE recombination hotspot file (single chromosome)
    kwargs - arguments to oxstats_haplotypes

    """

    with decode_recombination_hotspots(hotspot_file) as hotspots, \
         oxstats_haplotypes(**kwargs) as data:
        for pos, rate in hotspots:
            block = itertools.takewhile(lambda x: int(x[0][1]) < pos, data)
            yield list(zip(*block)), rate
        yield list(zip(*data)), 0

def _reconstruct(mosaic, haplotypes, center=True):
    """Return centered genotypes (n x p) according to ancestor pointers

    haplotypes - list of list of haplotypes (p x k)
    mosaic - list of ancestor pointers (2n x 1; values in 1..k)

    """
    n = mosaic.shape[0] / 2
    w = numpy.array(haplotypes, dtype='int8').T[mosaic]  # 2n x p
    x = w.reshape(n, -1, 2).sum(axis=2)
    if center:
        x = x.astype('float32')
        x -= x.mean(axis=0)
    return x

def sample_events(seed, n, hotspot_file, p_causal=0.5, n_per_window=1,
                  **kwargs):
    """Return recombination events and estimated liabilities

    The goal is to recombine reference haplotypes to generate genotypes and
    phenotypes for arbitrary size samples in two passes through the reference
    data.

    In the first pass, sample causal variants, effect sizes, and recombinations
    (at recombination hotspots). Recombinations are stored as an initial
    configuration plus sparse updates. Causal variants and effect sizes are
    stored as genetic values (since we won't need the actual effects later).

    mosaic - initial ancestor pointers (2n x 1; values in 1..k)
    p_causal - probability each block contains a causal variant
    n_per_window - number of causal variants per block
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    kwargs - arguments to oxstats_haplotypes

    """
    mosaic = None
    events = []
    y = numpy.zeros(n)
    for (_, haplotypes), rate in blocks(hotspot_file, **kwargs):
        if len(haplotypes) == 0:
            break
        haplotypes = numpy.array(haplotypes, dtype='int8')
        p, k = haplotypes.shape
        if mosaic is None:
            R.seed(seed)
            mosaic = numpy.arange(2 * n)
        if R.rand() < p_causal:
            causal = R.choice(numpy.arange(p), size=n_per_window, replace=False)
            effects = R.normal(size=n_per_window)
            x = _reconstruct(mosaic, haplotypes)
            y += x[:,causal].dot(effects)
        events.append([])
    return y, events

_logsf = scipy.stats.chi2(1).logsf

def compute_marginal_stats(x, y):
    """Compute marginal association statistics for SNPs x against phenotype y

    Assume Gaussian phenotype with mean 0, and all SNPs centered, and perform
    univariate linear regressions in parallel for each SNP.

    """
    n, p = x.shape
    var = numpy.diag(x.T.dot(x)) + 1e-8  # Needed for monomorphic SNPs
    b = y.T.dot(x).T / var
    s = ((y ** 2).sum() - b ** 2 * var) / (n - 1)
    se = numpy.sqrt(s / var)
    stat = numpy.square(b / se)
    logp = -_logsf(stat)
    return b, se, logp

def marginal_association(seed, y, events, hotspot_file, pve=0.5, eps=1e-8,
                         **kwargs):
    """Output marginal association statistics for a Gaussian phenotype

    In the second pass, reconstruct the sample genotypes locally and compute
    marginal association statistics in parallel using matrix operations. Output
    summary statistics in BED format, with -log10(p) as the score and a
    placeholder (',') for chromosome. The placeholder is needed to avoid
    dependency on file name formatting, which changes with each 1KG release

    Assume the same random seed is used as in the first pass (to get the same
    starting configuration).

    mosaic - initial ancestor pointers (2n x 1; values in 1..k)
    events - list of (haplotype index, new ancestor) updates to mosaic
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    pve - target proportion of variance explained
    eps - tolerance (for SNP variance at monomorphic SNPs)
    kwargs - arguments to oxstats_haplotypes

    """
    n = y.shape[0]
    mosaic = None
    numpy.seterr(all='warn')
    for ((legend, haplotypes), _), event in zip(blocks(hotspot_file, **kwargs), events):
        if mosaic is None:
            k = len(haplotypes[0])
            R.seed(seed)
            mosaic = numpy.arange(2 * n)
        x = _reconstruct(mosaic, haplotypes)
        b, se, logp = compute_marginal_stats(x, y)
        for l, b_j, s_j, p in zip(legend, b, se, logp):
            name, pos, a0, a1, *_ = l
            print(name, pos, a0, a1, b_j, s_j, p)

def output_genotypes(seed, y, events, chromosome, hotspot_file,
                     file=sys.stdout, **kwargs):
    """Output reconstructed genotypes in VCF format

    Assume random seed is the same used to sample recombination events (to
    ensure same starting mosaic)

    seed - random seed (integer)
    events - recombination events
    chromosome - chromosome code
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    kwargs - arguments to oxstats_haplotypes

    """
    n = y.shape[0]
    delim='\t'
    print("##fileformat=VCFv4.1", file=file)
    print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
          'FORMAT', delim.join('GEN{}'.format(i) for i in range(n)), sep=delim,
          file=file)
    mosaic = None
    coding = ('0/0', '0/1', '1/1')
    for ((legend, haplotypes), _), event in zip(blocks(hotspot_file, **kwargs), events):
        if mosaic is None:
            k = len(haplotypes[0])
            R.seed(seed)
            mosaic = numpy.arange(2 * n)
        x = _reconstruct(mosaic, haplotypes, center=False).T
        for l, x_j in zip(legend, x):
            name, pos, a0, a1, *_ = l
            print(chromosome, pos, name, a0, a1, '.', '.', '.', 'GT',
                  delim.join(coding[x_ij] for x_ij in x_j), sep=delim,
                  file=file)

def output_phenotype(y, file=sys.stdout):
    """Output phenotype in plink format"""
    for i, y_i in enumerate(y):
        print('0', 'GEN{}'.format(i), y_i, file=file)

def combine_genetic_values(seed, files, pve=0.5):
    """Add Gaussian noise to generate a phenotype with target PVE

    seed - random seed
    files - pickled (y, events) tuples
    pve - target PVE

    """
    R.seed(seed)
    genetic_values = []
    for f_ in files:
        with open(f_, 'rb') as f:
            y, _ = pickle.load(f)
            genetic_values.append(y)
    y = numpy.array(genetic_values).sum(axis=0)
    # Use the realized genetic variance to set the noise scale. In other
    # simulations, we use the population value of the genetic variance (based
    # on MAF and effect size) instead.
    y += R.normal(scale=numpy.sqrt((1 / pve - 1) * y.var()), size=y.shape)
    y -= y.mean()
    y /= y.std()
    return y

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
    y = numpy.zeros(len(samples))
    for dose in sample_uniform(probs, **kwargs):
        dose -= dose.mean()
        y += dose * R.normal(size=y.shape[0])
    y += R.normal(scale=numpy.sqrt((1 / pve - 1) * y.var()))
    y -= y.mean()
    y /= y.std()
    for i, y_i in enumerate(y):
        print(' '.join(samples[i]), y_i)

if __name__ == '__main__':
    with contextlib.ExitStack() as stack:
        data = [stack.enter_context(oxstats_genotypes(*args)) for args in kwise(sys.argv[1:], 2)]
        samples = list(itertools.chain.from_iterable(s for _, _, s, _ in data))
        merged = merge_oxstats([d for _, _, _, d in data])
        h1, h2, *_ = data[0]
        print(' '.join(h1), 'pheno')
        print(' '.join(h2), 'P')
        simulate_null(samples, merged)
