""" Prune SNPs according to LD

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import collections
import gzip
import operator
import sys

with gzip.open(sys.argv[1], 'rt') as f:
    types = [str, int, int, str, float]
    scores = (line.split() for line in f)
    parsed = ([f(x) for f, x in zip(types, snp)] for snp in scores)
    snps = {name: score for _, _, _, name, score in parsed}

with gzip.open(sys.argv[2], 'rt') as f:
    thresh = float(sys.argv[3]) if len(sys.argv) > 3 else 0
    ref_ld = (line.split() for line in f)
    ld = collections.defaultdict(list)
    for s in snps:
        ld[s].append((s, 1))
    for k, v, r, d in ref_ld:
        if float(r) >= thresh and k != v and k in snps and v in snps:
            if snps[k] < snps[v]:
                k, v = v, k
            ld[k].append((v, r))

for tag, score in sorted(snps.items(), key=operator.itemgetter(1), reverse=True):
    if tag in ld:
        for tagged, r in ld.pop(tag):
            if tagged in ld:
                print(tag, tagged, r)
                ld.pop(tagged)
