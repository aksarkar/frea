"""Simulate phenotypes with different genetic architectures.

To generate realistic genotypes, we adapt the approach of Loh et al., Nat Genet
2015 (Supplementary Figure 4) to sample mosaics of founder individuals. We use
Thousand Genomes phased reference haplotypes as the founders and sample
recombinations at known recombination hotspots.

The simulation proceeds in three passes over the reference data:

1. Sample recombination events and genetic values
2. Compute left singular vectors of the reconstructed genotypes
3. Compute marginal association statistics

The key design decision is to only store the sparse updates to the mosaic
induced by recombination events, minimizing the storage requirement (which in
turn improves the performance on NFS).

To control the effect of population structure on marginal summary statistics,
we restrict founder haplotypes to a single cohort of 1KG. We additionally
implement streaming truncated singular value decomposition using rank one
updates through the reconstructed genotypes, optionally subsampling within LD
blocks.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import random
import re
import sys

import numpy
import numpy.random as R

from .formats import *

numpy.seterr(all='warn')

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

def _reconstruct(mosaic, haplotypes, center=True):
    """Return genotypes (n x p) according to ancestor pointers

    If genotypes are centered, casts to float32; otherwise, casts to int8.

    haplotypes - list of list of haplotypes (p x k)
    mosaic - list of ancestor pointers (2n x 1; values in 1..k)
    center - True if returned genotypes should be centered

    """
    n = mosaic.shape[0] // 2
    w = numpy.array(haplotypes, dtype='int8').T[mosaic]  # 2n x p
    x = w.reshape(n, -1, 2).sum(axis=2)
    if center:
        x = x.astype('float32')
        x -= x.mean(axis=0)
    return x

def reconstruction_pass(seed, n, events, hotspot_file, **kwargs):
    """Yield updated mosaic, legend, and haplotype entries

    seed - random seed
    n - target sample size
    events - list of (haplotype index, new ancestor) updates to mosaic
    hotspot_file - deCODE recombination hotspot file (single chromosome)
    kwargs - arguments to oxstats_haplotypes

    """
    mosaic = None
    for ((legend, haplotypes), _), (hits, ancestors) in zip(blocks(hotspot_file, **kwargs), events):
        if mosaic is None:
            R.seed(seed)
            mosaic = R.randint(0, len(haplotypes[0]), size=2 * n)
        else:
            mosaic[hits] = ancestors
        yield mosaic, legend, haplotypes

def thin_svd(k=20, n_per_block=None, partial_svd=None, debug=False, **kwargs):
    """Return the partial rank-k truncated SVD of the genotypes

    We recover the top k principal components truncated SVD in one pass through
    the data using rank one updates to a partial solution (Brand, Lin Alg Appl
    2006). This method is required when the genotype matrix cannot fit in
    memory and eigendecomposition of the kernel matrix is impossible. We
    discard the right singular vectors (SNP loadings).

    k - number of principal components to return
    n_per_block - number of variants to subsample in SVD (default: don't subsample)
    partial_svd - current factorization (U, S)
    kwargs - arguments to reconstruction_pass

    """
    if partial_svd is None:
        U = numpy.matrix(numpy.zeros((kwargs['n'], k)))
        S = numpy.matrix(numpy.zeros((k, k)))
    else:
        U, S = partial_svd
    for mosaic, _, haplotypes in reconstruction_pass(**kwargs):
        haplotypes = numpy.array(haplotypes, dtype='int8')
        p, _ = haplotypes.shape
        if n_per_block is None or n_per_block >= p:
            subsample = numpy.arange(p)
        else:
            subsample = R.choice(numpy.arange(p), size=n_per_block, replace=False)
        for h in haplotypes[subsample]:
            x = _reconstruct(mosaic, h)
            # Eq. (6)
            m = U.T * x
            P = x - U * m
            Ra = numpy.linalg.norm(p) + 1e-8
            P /= Ra
            # Eq. (9)
            K = numpy.bmat([[S, m], [numpy.zeros((1, k)), numpy.matrix(Ra)]])
            UK, SK, VK = numpy.linalg.svd(K, full_matrices=False)
            # Eq. (5); n x (k + 1) * (k + 1) x k
            U = numpy.bmat([U, P]) * UK[:, :k]
            S = numpy.matrix(numpy.diag(SK[:k]))
    return U, S

_sf = scipy.stats.chi2(1).sf

def marginal_association(y, covars=None, eps=1e-8, chunk_size=1000, **kwargs):
    """Yield marginal association statistics for a Gaussian phenotype

    Reconstruct the sample genotypes locally and fit linear regression models
    in parallel using matrix operations (Sikorska et al., BMC Bioinformatics
    2013). We do not include an intercept in the model because we assume
    genotypes and phenotypes are centered.

    y - phenotype vector (n x 1)
    covars - continuous covariates (n x k)
    eps - tolerance (for SNP variance at monomorphic SNPs)
    chunk_size - number of models to fit in parallel
    kwargs - arguments to reconstruction_pass

    """
    if covars is not None:
        # Regress out covariates
        covars -= numpy.array(covars).mean(axis=0)
        C = numpy.linalg.pinv(covars)
        y -= numpy.einsum('ij,jk,k->i', covars, C, y)
    for mosaic, legend, haplotypes in reconstruction_pass(**kwargs):
        for k, g in itertools.groupby(enumerate(zip(legend, haplotypes)), key=lambda x: int(x[0] // chunk_size)):
            legend, haplotypes = zip(*[x[1] for x in g])
            x = _reconstruct(mosaic, haplotypes)
            if covars is not None:
                x -= covars.dot(C.dot(x))
            n, p = x.shape
            var = numpy.diag(x.T.dot(x)) + 1e-8  # Needed for monomorphic SNPs
            beta = y.T.dot(x).T / var
            df = n - 1
            if covars is not None:
                df -= covars.shape[1]
            s = ((y ** 2).sum() - beta ** 2 * var) / df
            se = numpy.sqrt(s / var)
            stat = numpy.square(beta / se)
            logp = -numpy.log10(_sf(stat))
            for l, b, s, p in zip(legend, beta, se, logp):
                name, pos, a0, a1, *_ = l
                yield name, int(pos), a0, a1, b, s, p

def output_genotypes(chromosome, file=sys.stdout, **kwargs):
    """Output reconstructed genotypes in VCF format

    Assume random seed is the same used to sample recombination events (to
    ensure same starting mosaic)

    chromosome - chromosome code
    n - number of samples
    kwargs - arguments to reconstruction_pass

    """
    delim='\t'
    print("##fileformat=VCFv4.1", file=file)
    print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
          'FORMAT', delim.join('GEN{}'.format(i) for i in range(kwargs['n'])), sep=delim,
          file=file)
    coding = ('0/0', '0/1', '1/1')
    for mosaic, legend, haplotypes in reconstruction_pass(**kwargs):
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
    return y

def output_oxstats_pheno():
    parser = argparse.ArgumentParser()
    parser.add_argument('state', help='Output of frea-combine-pheno')
    parser.add_argument('sample_file', nargs='+', help='Individual cohort sample files')
    args = parser.parse_args()
    with open(args.state, 'rb') as f:
        y = pickle.load(f)
    it = iter(y)
    for filename in args.sample_file:
        with open(filename) as f, open(re.sub('.sample', '.pheno', filename), 'w') as g:
            data = (line.strip().split() for line in f)
            print(' '.join(next(data)), 'pheno', file=g)
            print(' '.join(next(data)), 'P', file=g)
            for row in data:
                print(' '.join(row[:4]), next(it), file=g)
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
