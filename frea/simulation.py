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

from .formats import *

def sample_chromosome_mosaic(n, sample_file, legend_file, haps_file, p_causal=0.5,
                             window_size=1e6, n_per_window=1, **kwargs):
    """Sample pointers to ancestors at fixed recombination points

    n - number of descendant samples
    sample_file - reference sample information
    legend_file - reference variant information
    haps_file - reference haplotypes
    p_causal - probability each window contains a causal variant
    window_size - window between fixed recombination points (and thinning causal variants)
    n_per_window - number of causal variants per window
    **kwargs - arguments to frea.formats.oxstats_haplotypes

    """
    with oxstats_haplotypes(sample_file, legend_file, haps_file, **kwargs) as data:
        for _, g in itertools.groupby(data, key=lambda x: int(int(x[0][1]) / window_size)):
            window = list(g)
            k = len(window[0][1])
            if random.random() < p_causal:
                w = [l[:5] for l, _ in window]
                causal = random.sample(w, n_per_window)
                effects = R.normal(size=len(n_per_window))
            else:
                causal = []
                effects = []
            yield R.randint(0, k, size=2 * n), causal, effects

def sample_mosaic(n, ref_dir='.', **kwargs):
    """Yield a mosaic of reference haplotypes.

    n - number of descendent samples
    kwargs - keyword arguments to parse_oxstats_haplotypes

    """
    J = os.path.join
    for chr_ in [str(i) for i in range(23)]:
        sample_file = J(ref_dir, 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample')
        legend_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(chr_))
        haps_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'.format(chr_))
        yield from sample_chromosome_mosaic(sample_file, legend_file, haps_file, **kwargs)

def reconstruct_chromosome_mosaic(ancestors, sample_file, legend_file,
                                  haps_file, window_size=1e6, **kwargs):
    """Reconstruct genotypes from recorded recombinations

    ancestors - pointer to haplotypes to copy between recombinations
    sample_file - reference sample information
    legend_file - reference variant information
    haps_file - reference haplotypes
    window_size - window between recombination points
    kwargs - keyword arguments to parse_oxstats_haplotypes

    """
    with oxstats_haplotypes(sample_file, legend_file, haps_file, **kwargs) as data:
        for _, g in itertools.groupby(data, key=lambda x: int(int(x[0][1]) / window_size)):
            for _, h in g:
                haplotypes = [h[i] for i in ancestors]
                yield [sum(h) for h in kwise(haplotypes, 2)]

def reconstruct_mosaic(ancestors, ref_dir='.', **kwargs):
    """Yield reconstructed genotypes from recorded recombinations

    ancestors - pointer to haplotypes to copy between recombinations
    kwargs - keyword arguments to parse_oxstats_haplotypes

    """
    J = os.path.join
    for chr_ in [str(i) for i in range(23)]:
        sample_file = J(ref_dir, 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample')
        legend_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(chr_))
        haps_file = J(ref_dir, 'ALL.chr{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'.format(chr_))
        yield from reconstruct_chromosome_mosaic(ancestors, sample_file, legend_file, haps_file, **kwargs)

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
