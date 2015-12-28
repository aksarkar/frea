import contextlib
import functools
import gzip
import math
import os
import signal
import sys

import scipy.stats

from .join import lookup
from .process import chromosome
from ..formats import output_impg_zscores

_isf = scipy.stats.chi2(1).isf
_z = lambda x: math.sqrt(_isf(x))

def output_check_alleles(bed_entry, ref_entry, file=sys.stdout):
    chr_, _, p, b0, b1 = bed_entry
    name, pos, a0, a1, _ = ref_entry
    z = _z(p)
    if (a0, a1) == (b1, b0):
        z *= -1
        b0, b1 = a0, a1
    if (a0, a1) == (b0, b1):
        print(output_impg_zscores(chr_, z, name, pos, a0, a1), file=file)
    else:
        print(chr_, pos, name, a0, a1, b0, b1, file=sys.stderr)

def _key(row):
    return row[0], int(row[1]), float(row[4]), row[5], row[6]

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    data = (line.split() for line in sys.stdin)
    parsed = (_key(row) for row in data)
    for k, g in itertools.groupby(parsed, key=operator.itemgetter(0)):
        with open('{}.txt'.format(k[3:]), 'w') as outfile:
            output_fn = functools.partial(output_check_alleles, file=outfile)
            print('id', 'pos', 'ref', 'alt', 'z', file=outfile)
            lookup(g, output_fn)
