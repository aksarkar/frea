"""Score LD blocks for functional enrichment

Expects output of bedtools intersect -c on stdin. Writes out (-log(p), score)
on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import operator
import signal
import sys

if __name__ == '__main__':
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    data = (line.split() for line in sys.stdin)
    for k, g in itertools.groupby(data, key=operator.itemgetter(3)):
        snps = list(g)
        logp = max(s[4] for s in snps)
        score = sum(s[5] for s in snps) / len(snps)
        print(logp, score)
