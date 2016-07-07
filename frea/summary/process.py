import math
import operator

import scipy.stats

from .lookup import *
from ..formats import *

def oxstats_snptest_to_ucsc_bed(filename, chr_, min_info=0.8):
    """Convert SNPTEST output to BED format

    Assume SNPTEST was run on 1KG-imputed genotypes and write out BED entries
    with -log10(p) as the score.

    """
    with gzip.open(filename, 'rt') as f:
        drop_comments = (line for line in f if not line.startswith('#'))
        data = (line.split() for line in drop_comments)
        header = {x: i for i, x in enumerate(next(data))}
        drop_missing = (row for row in data if row[header['frequentist_add_pvalue']] != 'NA')
        chr_ = chromosome(int(chr_))
        for row in drop_missing:
            pos = int(row[header['position']])
            rsid = row[header['rsid']]
            a0 = row[header['alleleA']]
            a1 = row[header['alleleB']]
            logp = -math.log10(float(row[header['frequentist_add_pvalue']]))
            print(chr_, pos, pos + len(a0),
                  get_pouyak_name(chr_, rsid, pos, a0, a1),
                  logp, sep='\t')

def ucsc_to_impg(ref_dir='.'):
    """Convert lifted over z-scores to ImpG input format"""
    data = (line.split('\t') for line in sys.stdin)
    parsed = parse(ucsc_bed_format, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        with open('{}.txt'.format(k[3:]), 'w') as f:
            print('id', 'pos', 'ref', 'alt', 'z', file=f)
            output_ = functools.partial(output, output_fn=output_impg_zscores, score_fn=float, file=f)
            lookup(g, output_fn=output_, ref_dir=ref_dir)

def _fix_names(fmt, output_fn, **kwargs):
    data = (line.strip().split('\t') for line in sys.stdin)
    parsed = parse(fmt, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        lookup(g, output_fn, **kwargs)

def fix_ucsc_bed_names(**kwargs):
    """Fix the names of a BED file"""
    if 'fmt' not in kwargs:
        kwargs['fmt'] = ucsc_bed_format
    if 'output_fn' not in kwargs:
        kwargs['output_fn'] = functools.partial(output, score_fn=str)
    _fix_names(**kwargs)

def fix_plink_bim_names(**kwargs):
    def output(bim, ref):
        chr_, affy_name = bim[:2]
        name, pos, a0, a1, _ = ref
        print(affy_name, get_pouyak_name(chr_, name, pos, a0, a1))
    _fix_names(plink_format, output, input_join_key=operator.itemgetter(3), **kwargs)

def _reformat(filename, key):
    for row in summary_stats(filename, key):
        print(row[0], row[1], row[1] + 1, row[2], row[3], sep='\t')

def diagram(filename):
    def key(row):
        return chromosome(row[1]), int(row[2]), int(row[2]) + 1, row[0], float(row[5]), row[3], row[4]
    for row in summary_stats(filename, key):
        # Use the given end for indels
        print(*row)

def cardiogram(filename):
    def key(row):
        chromosome, pos = row[1].split(':')
        return chromosome, int(pos), row[0], zscore(float(row[5]), math.exp(float(row[7])))
    _reformat(filename, key)

def ra(filename):
    _reformat(filename, lambda x: (chromosome(x[0]), int(x[2]), x[1], x[-4]))

def ra_subcohorts(filename):
    _reformat(filename, lambda x: (chromosome(x[0]), int(x[2]), x[1], '|'.join(x[11:17] + x[-4:-3])))

def iibdgc(filename, ref_dir):
    def key(row):
        return chromosome(row[0]), int(row[2]), float(row[10])
    lookup(sorted(summary_stats(filename, key), key=operator.itemgetter(0, 1)), ref_dir=ref_dir)

def bip(filename):
    _reformat(filename, lambda x: (chromosome(x[1]), int(x[2]), x[0], zscore(float(x[7]), float(x[5]))))

def scz(filename, ref_dir):
    def key(row):
        return row[0], int(row[4]), float(row[8])
    lookup(sorted(summary_stats(filename, key), key=operator.itemgetter(0, 1)), ref_dir=ref_dir)
