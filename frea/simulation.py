"""Simulate phenotypes with different genetic architectures.

To generate realistic genotypes, we adapt the approach of Loh et al., Nat Genet
2015 (Supplementary Figure 4) to sample mosaics of founder individuals. We use
Thousand Genomes phased reference haplotypes as the founders and sample
recombinations at known recombination hotspots.

The simulation proceeds in two passes over the reference data:

1. Sample recombination events and genetic values
2. Compute marginal association statistics

The key design decision is to only store the sparse updates to the mosaic
induced by recombination events, minimizing the storage requirement (which in
turn improves the performance on NFS).

To control the effect of population structure on marginal summary statistics,
we restrict founder haplotypes to a single cohort of 1KG. However, it is
possible (Brand, Lin Alg Appl 2006) to compute a rank-k truncated SVD in one
additional pass through the data.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import random
import sys

import numpy
import numpy.random as R

from .formats import *

def blocks(hotspot_file, **kwargs):
    """Yield LD blocks defined by recombination hotspots

    Blocks are a list of OXSTATS reference haplotypes and recombination rates
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
    """Return genotypes (n x p) according to ancestor pointers

    If genotypes are centered, casts to float32; otherwise, casts to int8.

    haplotypes - list of list of haplotypes (p x k)
    mosaic - list of ancestor pointers (2n x 1; values in 1..k)
    center - True if returned genotypes should be centered

    """
    n = mosaic.shape[0] / 2
    w = numpy.array(haplotypes, dtype='int8').T[mosaic]  # 2n x p
    x = w.reshape(n, -1, 2).sum(axis=2)
    if center:
        x = x.astype('float32')
        x -= x.mean(axis=0)
    return x

def sample_events(seed, n, hotspot_file, p_causal=0.5, n_per_block=1,
                  **kwargs):
    """Return recombination events and estimated liabilities

    The goal is to recombine reference haplotypes to generate genotypes and
    phenotypes for arbitrary size samples in two passes through the reference
    data.

    In the first pass, sample causal variants, effect sizes, and recombinations
    (at recombination hotspots). Recombinations are stored as an initial
    configuration plus sparse updates. Causal variants and effect sizes are
    stored as genetic values (since we won't need the actual effects later).

    seed - random seed
    n - target sample size
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    p_causal - probability each block contains a causal variant
    n_per_block - number of causal variants per block
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
            mosaic = R.randint(0, k, size=2 * n)
        if R.rand() < p_causal:
            causal = R.choice(numpy.arange(p), size=n_per_block, replace=False)
            effects = R.normal(size=n_per_block)
            # (k x n_per_block) x (n_per_block x 1) => (2n x 1) => (n x 1)
            value = haplotypes[causal].T.dot(effects)[mosaic].reshape(-1, 2).sum(axis=1)
            value -= value.mean()
            y += value
        # TODO: without rejection sampling, the sampled rate of recombinations
        # will be biased towards zero because we sample ancestors with
        # replacement (so a sampled event might not actually switch
        # haplotypes). But accurately simulating recombination is not the goal
        # of this simulation, so just verify that the summary statistics are
        # reasonable
        hits = numpy.where(R.uniform(size=2 * n) < rate / 100)[0]
        ancestors = R.randint(0, k, size=hits.shape)
        events.append((hits, ancestors))
        mosaic[hits] = ancestors
    return y, events

_sf = scipy.stats.chi2(1).sf

def compute_marginal_stats(x, y):
    """Compute marginal association statistics for SNPs x against phenotype y

    Assume Gaussian phenotype with mean 0, and all SNPs centered, and perform
    univariate linear regressions in parallel for each SNP.

    x - dosage matrix (n x p)
    y - phenotype vector (n x 1)

    """
    n, p = x.shape
    var = numpy.diag(x.T.dot(x)) + 1e-8  # Needed for monomorphic SNPs
    b = y.T.dot(x).T / var
    s = ((y ** 2).sum() - b ** 2 * var) / (n - 1)
    se = numpy.sqrt(s / var)
    stat = numpy.square(b / se)
    logp = -numpy.log10(_sf(stat))
    return b, se, logp

def marginal_association(seed, y, events, hotspot_file, pve=0.5, eps=1e-8,
                         chunk_size=1000, **kwargs):
    """Output marginal association statistics for a Gaussian phenotype

    In the second pass, reconstruct the sample genotypes locally and compute
    marginal association statistics in parallel using matrix operations. Output
    summary statistics in BED format, with -log10(p) as the score and a
    placeholder (',') for chromosome. The placeholder is needed to avoid
    dependency on file name formatting, which changes with each 1KG release

    Assume the same random seed is used as in the first pass (to get the same
    starting configuration).

    seed - random seed
    y - phenotype vector (n x 1)
    events - list of (haplotype index, new ancestor) updates to mosaic
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    pve - target proportion of variance explained
    eps - tolerance (for SNP variance at monomorphic SNPs)
    chunk_size - number of models to fit in parallel
    kwargs - arguments to oxstats_haplotypes

    """
    n = y.shape[0]
    mosaic = None
    numpy.seterr(all='warn')
    for ((legend, haplotypes), _), (hits, ancestors) in zip(blocks(hotspot_file, **kwargs), events):
        if mosaic is None:
            k = len(haplotypes[0])
            R.seed(seed)
            mosaic = R.randint(0, k, size=2 * n)
        for k, g in itertools.groupby(enumerate(haplotypes),
                                      key=lambda x: int(x[0] // chunk_size)):
            x = _reconstruct(mosaic, [h for _, h in g])
            beta, se, logp = compute_marginal_stats(x, y)
            for l, p in zip(legend, logp):
                name, pos, a0, a1, *_ = l
                print(output_ucsc_bed(',', p, name, int(pos), a0, a1))
        mosaic[hits] = ancestors

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
    for ((legend, haplotypes), _), (hits, ancestors) in zip(blocks(hotspot_file, **kwargs), events):
        if mosaic is None:
            k = len(haplotypes[0])
            R.seed(seed)
            mosaic = R.randint(0, k, size=2 * n)
        x = _reconstruct(mosaic, haplotypes, center=False).T
        for l, x_j in zip(legend, x):
            name, pos, a0, a1, *_ = l
            print(chromosome, pos, name, a0, a1, '.', '.', '.', 'GT',
                  delim.join(coding[x_ij] for x_ij in x_j), sep=delim,
                  file=file)
        mosaic[hits] = ancestors

def output_phenotype(y, file=sys.stdout):
    """Output phenotype in plink format"""
    for i, y_i in enumerate(y):
        print('0', 'GEN{}'.format(i), y_i, file=file)

def combine_genetic_values(seed, values, pve=0.5):
    """Add Gaussian noise to generate a phenotype with target PVE

    seed - random seed
    values - m x n array of genetic values (m chromosomes)
    pve - target PVE

    """
    R.seed(seed)
    y = values.sum(axis=0)
    # Use the realized genetic variance to set the noise scale. In other
    # simulations, we use the population value of the genetic variance (based
    # on MAF and effect size) instead.
    y += R.normal(scale=numpy.sqrt((1 / pve - 1) * y.var()), size=y.shape)
    y -= y.mean()
    y /= y.std()
    return y
