"""Parsers for common formats

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import contextlib
import functools
import gzip
import itertools
import operator
import math
import sys

import scipy.stats

from .algorithms import join, kwise, moments

ucsc_bed_format = [str, int, int, str, float]
impg_format = [str, int, str, str, float, float]
legend_format = [str, int, str, str, str]
plink_format = [lambda x: 'chr{}'.format(x), str, int, int, str, str]

def get_pouyak_name(chr_, name, pos, a0, a1):
    if len(a0) <= len(a1):
        start = pos
        end = pos
        if len(a0) == len(a1):
            delta = a1
        else:  # insertion
            delta = a1[len(a0):]
    else:  # deletion
        start = pos + len(a1)
        end = pos + len(a0) - len(a1)
        delta = ""
    return '|'.join([name, chr_, str(end), str(end), delta])

def get_plink_alleles(a0, a1):
    if len(a0) > 1:
        a0 = 'I{}'.format(len(a0))
        a1 = 'D'
    elif len(a1) > 1:
        a0 = 'D'
        a1 = 'I{}'.format(len(a1))
    return a0, a1

def parse(types, iterable):
    for row in iterable:
        yield [f(x) for f, x in zip(types, row)]

def output_ucsc_bed(chr_, score, name, pos, a0, a1):
    return '\t'.join([chr_, str(pos - 1), str(pos + len(a0) - 1), get_pouyak_name(chr_, name, pos, a0, a1), str(score)])

def output_impg_zscores(chr_, score, name, pos, a0, a1):
    return ' '.join([name, str(pos), a0, a1, str(score)])

def output_plink_bim(chr_, score, name, pos, a0, a1):
    return '\t'.join([chr_, get_pouyak_name(chr_, name, pos, a0, a1), "0", str(pos), a0, a1])

_isf = scipy.stats.chi2(1).isf
_logsf = scipy.stats.chi2(1).logsf

def zscore(p, odds):
    z = math.sqrt(_isf(float(p)))
    if float(odds) < 1:
        z *= -1
    return z

def z_to_p(z):
    return -_logsf(z * z)

def logp(p):
    return -math.log10(float(p))

def chromosome(n):
    return 'chr{}'.format(n)

def summary_stats(filename, input_key, skip_header=True):
    """Parse summary statistics

    input_key - function that returns (chromosome, position, id, z)
    """
    autosomes = set('chr{}'.format(x) for x in range(1, 23))
    with gzip.open(filename, 'rt') as f:
        if skip_header:
            next(f)
        data = (line.split() for line in f)
        for row in data:
            parsed = input_key(row)
            if parsed[0] in autosomes:
                yield parsed

def oxstats_gen_to_dosage(probs):
    probs = [float(p) for p in probs]
    return [2 * c + b for (a, b, c) in kwise(probs, 3)]

def parse_oxstats(data):
    entries = (line.split() for line in data)
    for row in entries:
        row[2] = int(row[2])
        yield row

@contextlib.contextmanager
def oxstats_genotypes(sample_file, gen_file):
    """Return the list of samples and a generator which yields genotype
probabilities.

    Expects data in OXSTATS gen format (not bgen). This implementation does
    not allow random access to avoid memory issues.

    """
    with open(sample_file) as f:
        data = (line.split() for line in f)
        h1 = next(data)
        h2 = next(data)
        samples = list(data)
    with gzip.open(gen_file, 'rt') as f:
        yield h1, h2, samples, parse_oxstats(f)

def _merge_oxstats(seq1, seq2):
    """Merge lines of OXSTATS genotypes, matching on SNP attributes"""
    for a, b in join(seq1, seq2, key1=operator.itemgetter(2)):
        if a[:5] == b[:5]:
            yield a + b[5:]

def merge_oxstats(iterables):
    """Yield merged genotypes from iterables

    iterables -parsed OXSTATS data

    """
    for row in functools.reduce(_merge_oxstats, iterables):
        yield row

def parse_oxstats_haps(samples, legend_iterable, haps_iterable, group='EUR',
                       snps_only=True, min_maf=0.01):
    """Return OXSTATS haplotypes which meet the specific criteria

    samples - list of (sample, population, group, sex) records
    legend_iterable - file-like object
    haps_iterable - file-like object
    group - sample group
    snps_only - only keep SNPs
    min_maf - only keep variants with MAF above threshold

    """
    keep = [x[2] == group for x in samples]
    legend_entries = (line.split() for line in legend_iterable)
    header = next(legend_entries)
    maf_col = header.index('{}.maf'.format(group.lower()))
    haps_entries = (line.split() for line in haps_iterable)
    for l, h in zip(legend_entries, haps_entries):
        if (not snps_only or l[4] == 'SNP') and float(l[maf_col]) > min_maf:
            filtered_haps = [pair for pair, keep_ in zip(kwise(h, 2), keep) if keep_]
            ignore_homologs = itertools.chain.from_iterable(filtered_haps)
            parsed_haps = [int(x) for x in ignore_homologs]
            yield l[:4], parsed_haps

@contextlib.contextmanager
def oxstats_haplotypes(sample_file, legend_file, haps_file, **kwargs):
    """Return the list of samples and a generator which yields haplotypes.

    By default, reads from stdin. Expects data in OXSTATS haplotypes

    This implementation does not allow random access to avoid memory issues.

    kwargs - keyword arguments to parse_oxstats_haps

    """
    with open(sample_file) as f:
        next(f)
        samples = [line.split() for line in f]
    with gzip.open(legend_file, 'rt') as f, gzip.open(haps_file, 'rt') as g:
        yield parse_oxstats_haps(samples, f, g)

@contextlib.contextmanager
def decode_recombination_hotspots(hotspot_file):
    """Return a generator which yields recombination hotspots.

    deCODE estimated recombination rates in 10kb windows and defined hotspots
    based on significant difference against flanking regions. We combine
    contiguous 10kb windows, take the midpoint as the recombination point, and
    take the mean standardized recombination rate as the recombination rate
    (percentage of meioses which recombine).

    """
    def merge(parsed):
        for _, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
            g = list(g)
            mean, _ = moments(x[-1] for x in g)
            midpoint = (g[0][2] + g[-1][3]) // 2
            yield midpoint, mean

    with gzip.open(hotspot_file, 'rt') as f:
        data = (line.split() for line in f)
        parsed = parse([str, str, int, int, float], data)
        yield merge(parsed)
