"""Filter IMPUTE2 output based on the information metric

We cannot use the imputation output directly because we impute each cohort
independently, then filter after computing quality on cases and controls
jointly. This implementation is faster than SNPTEST and avoids an apparent bug
in its gzip-handling code.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import contextlib
import functools
import operator
import gzip
import re
import signal
import sys

from .algorithms import join, kwise
from .formats import chromosome, get_pouyak_name

def _merge_oxstats(seq1, seq2):
    for a, b in join(seq1, seq2, key1=operator.itemgetter(2)):
        if a[:5] == b[:5]:
            yield a + b[5:]

def merge_oxstats(iterables):
    for row in functools.reduce(_merge_oxstats, iterables):
        yield row

def parse_oxstats(data):
    entries = (line.split() for line in data)
    for row in entries:
        row[2] = int(row[2])
        yield row

def info(p):
    """Compute the ratio between observed and complete information.

    This implementation follows the description in
    https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.v2.pdf

    """
    p = list(p)
    e = [a + 2 * b for _, a, b in kwise(p, 3)]
    f = [a + 4 * b for _, a, b in kwise(p, 3)]
    theta_hat = sum(e) / float(2 * len(e))
    info = 1
    if theta_hat > 0 and theta_hat < 1:
        info -= sum(f_ - e_ * e_ for e_, f_ in zip(e, f)) / (2 * len(e) * theta_hat * (1 - theta_hat))
    return info

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    chrom = chromosome(int(re.search(r'(\d+)', sys.argv[1]).group()))
    with contextlib.ExitStack() as stack:
        files = [stack.enter_context(gzip.open(f, 'rt')) for f in sys.argv[1:]]
        parsed = [parse_oxstats(f) for f in files]
        merged = merge_oxstats(parsed)
        for row in merged:
            print(get_pouyak_name(chrom, *row[1:5]), info(float(x) for x in row[5:]))
