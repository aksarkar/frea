"""Convert ImpG output to UCSC BED

Usage: python -m frea.process.process_impg CHROMOSOME

Assumes current working directory contains ImpG formatted output (including
r2pred) files with ImpG naming scheme (chrX.start-end.txt).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import functools
import os
import re
import signal
import sys

import scipy.stats

from ..algorithms import join
from ..formats import parse, impg_format
from .lookup import *

def impg_window(filename, min_r2pred=0.6):
    patt = re.compile(r'\d+\.(?P<start>\d+)-(?P<end>\d+)')
    m = patt.match(filename)
    target = '{}.txt'.format(m.group(0))
    start, end = m.groups()
    start = int(start)
    end = int(end) - 250000
    if start > 1:
        start += 250000
    if not os.path.exists(target):
        raise StopIteration
    with open(target) as f:
        data = (line.split() for line in f)
        next(data)
        parsed = parse(impg_format, data)
        for row in parsed:
            if row[1] >= start and row[1] < end and row[5] >= min_r2pred:
                yield row

def impg_chromosome(chromosome, min_r2pred=0.6):
    with open('../maps/{}.names'.format(chromosome)) as f:
        for filename in f:
            for row in impg_window(filename, min_r2pred):
                yield row

def process_impg(chromosome, ref_dir='.'):
    """Process per-window ImpG output"""
    sf = scipy.stats.chi2(1).sf
    k = 'chr{}'.format(chromosome)
    output_ = functools.partial(output, key=lambda x: (k, sf(math.pow(x[4], 2))))
    lookup(impg_chromosome(chromosome), output_fn=output_,
           input_sort_key=lambda x: k, input_join_key=operator.itemgetter(1, 2, 3),
           ref_join_key=operator.itemgetter(1, 2, 3), ref_dir=ref_dir)

def process_impz(chromosome, ref_dir='.'):
    """Process per-chromosome ImpG output"""
    sf = scipy.stats.chi2(1).sf
    k = 'chr{}'.format(chromosome)
    output_ = functools.partial(output, key=lambda x: (k, sf(math.pow(x[4], 2))))
    data = (line.split() for line in sys.stdin)
    next(data)
    parsed = parse(impg_format, data)
    filtered = (row for row in parsed if row[5] >= 0.6)
    lookup(filtered, output_fn=output_,
           input_sort_key=lambda x: k, input_join_key=operator.itemgetter(1, 2, 3),
           ref_join_key=operator.itemgetter(1, 2, 3), ref_dir=ref_dir)

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if len(sys.argv) > 2:
        ref_dir = sys.argv[2]
    else:
        ref_dir = '.'
    process_impg(sys.argv[1], ref_dir)
