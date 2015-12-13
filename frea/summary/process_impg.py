import functools
import re
import signal

import scipy.stats

from .join import *

def impg_window(filename, min_r2pred=0.6):
    patt = re.compile(r'\d+\.(?P<start>\d+)-(?P<end>\d+)')
    m = patt.match(filename)
    start, end = m.groups()
    start = int(start)
    end = int(end) - 250000
    if start > 1:
        start += 250000
    with open('{}.txt'.format(m.group(0))) as f:
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

def process_impg(chromosome):
    """Process ImpG output"""
    sf = scipy.stats.chi2(1).sf
    k = 'chr{}'.format(chromosome)
    output_ = functools.partial(output, key=lambda x: (k, sf(math.pow(x[4], 2))))
    lookup(impg_chromosome(chromosome), output_fn=output_,
           input_sort_key=lambda x: k, input_join_key=operator.itemgetter(1, 2, 3),
           ref_join_key=operator.itemgetter(1, 2, 3))

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    process_impg(sys.argv[1])
