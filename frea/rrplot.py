"""Bin binary annotations for RR plots

Usage: python bin.py PHENOTYPE FEATURE CELLTYPE BINSIZE [CUMULATIVE]

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import itertools
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

phenotype = sys.argv[1]
feature = sys.argv[2]
celltype = sys.argv[3]
binsize = int(sys.argv[4])
cumulative = len(sys.argv) > 5
data = [float(x) for x in sys.stdin]
raw_exp = itertools.repeat(binsize / (1 + len(data)) * sum(data))
bins = (data[i:i+binsize] for i in range(0, len(data), binsize))
raw_counts = (sum(b) for b in bins)
if cumulative:
    counts = itertools.accumulate(raw_counts)
    exp = itertools.accumulate(raw_exp)
else:
    counts = raw_counts
    exp = raw_exp
for i, (c, e) in enumerate(zip(counts, exp)):
    print((i + 1) * binsize, phenotype, celltype, feature, c, e, sep=',')
