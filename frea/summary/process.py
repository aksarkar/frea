import math
import operator

import scipy.stats

from .lookup import *
from ..formats import *

def ucsc_to_impg():
    """Convert lifted over z-scores to ImpG input format"""
    data = (line.split('\t') for line in sys.stdin)
    parsed = parse(ucsc_bed_format, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        with open('{}.txt'.format(k[3:]), 'w') as f:
            print('id', 'pos', 'ref', 'alt', 'z', file=f)
            output_ = functools.partial(output, output_fn=output_impg_zscores, score_fn=float, file=f)
            lookup(g, output_fn=output_)

def _fix_names(fmt, output_fn, **kwargs):
    data = (line.split('\t') for line in sys.stdin)
    parsed = parse(fmt, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        lookup(g, output_fn, **kwargs)

def fix_ucsc_bed_names():
    """Fix the names of a BED file"""
    _fix_names(ucsc_bed_format, output_ucsc_bed)

def fix_plink_bim_names():
    def output(bim, ref):
        chr_, affy_name = bim[:2]
        name, pos, a0, a1, _ = ref
        print(affy_name, get_pouyak_name(chr_, name, pos, a0, a1))
    _fix_names(plink_format, output, input_join_key=operator.itemgetter(3))

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
