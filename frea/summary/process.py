import math
import operator

import scipy.stats

from .lookup import *

def process_sim():
    """Process plink linear regression output"""
    fmt = [lambda x: "chr{}".format(x), str, int, str, str, int, float, float, float]
    filter_garbage = (line for line in sys.stdin if not line.startswith(" CHR"))
    lookup(parse(filter_garbage, fmt), output, input_join_key=operator.itemgetter(2))

def ucsc_to_impg():
    """Convert lifted over z-scores to ImpG input format"""
    data = (line.split('\t') for line in sys.stdin)
    parsed = parse(ucsc_bed_format, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        with open('{}.txt'.format(k[3:]), 'w') as f:
            print('id', 'pos', 'ref', 'alt', 'z', file=f)
            output_ = functools.partial(output, output_fn=output_impg_zscores, score_fn=float, file=f)
            lookup(g, output_fn=output_)

def fix_names():
    """Fix the names of a BED file"""
    data = (line.split('\t') for line in sys.stdin)
    parsed = parse(ucsc_bed_format, data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        lookup(g)

_isf = scipy.stats.chi2(1).isf

def zscore(p, odds):
    z = math.sqrt(_isf(float(p)))
    if float(odds) < 1:
        z *= -1
    return z

def logp(p):
    return -math.log10(float(p))

def summary_stats(filename, input_key, skip_header=True):
    """Parse summary statistics

    input_key - function that returns (chromosome, position, id, z)
    """
    with gzip.open(filename, 'rt') as f:
        if skip_header:
            next(f)
        data = (line.split() for line in f)
        for row in data:
            yield input_key(row)

def chromosome(n):
    return 'chr{}'.format(n)

def _helper(filename, key):
    for row in summary_stats(filename, key):
        print(row[0], row[1], row[1] + 1, row[2], row[3], sep='\t')

def diagram():
    def key(row):
        return chromosome(row[1]), int(row[2]), int(row[2]) + 1, row[0], float(row[5]), row[3], row[4]
    for row in summary_stats('/broad/compbio/aksarkar/data/gwas-summary-stats/diagram/DIAGRAMv3.2012DEC17.txt.gz', key):
        print(*row)

def cardiogram():
    def key(row):
        chromosome, pos = row[1].split(':')
        return chromosome, int(pos), row[0], zscore(float(row[5]), math.exp(float(row[7])))
    _helper('/broad/compbio/aksarkar/data/gwas-summary-stats/cardiogram/cardiogram_gwas_results/CARDIoGRAM_GWAS_RESULTS.txt.gz',
            key)

def ra():
    _helper('/broad/compbio/aksarkar/data/gwas-summary-stats/ra/ra.txt.gz',
            lambda x: (chromosome(x[0]), int(x[2]), x[1], x[-4]))

def iibdgc():
    key = lambda x: (chromosome(x[0]), int(x[2]), float(x[10]))
    def parse(iterable):
        for row in iterable:
            yield key(row)
    with gzip.open('/broad/compbio/aksarkar/data/gwas-summary-stats/iibdgc/iibdgc-trans-ancestry-summary-stats/EUR.CD.gwas.assoc.gz', 'rt') as f:
        next(f)
        data = (line.split() for line in f)
        lookup(sorted(parse(data), key=operator.itemgetter(0, 1)), output_fn=output)

def pgc():
    _helper('/broad/compbio/aksarkar/data/gwas-summary-stats/pgc/pgc.bip.2012-04/pgc.bip.full.2012-04.txt.gz',
            lambda x: (chromosome(x[1]), int(x[2]), x[0], zscore(float(x[7]), float(x[5]))))
