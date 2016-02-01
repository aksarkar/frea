"""Print minimal set of tags that covers a set of SNPs

Usage: python tag.py TAGS SNPS

TAGS and SNPS are lists of IDs (one per line)

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import sys

with open(sys.argv[1]) as f:
    tags = set(line.strip() for line in f)

with open(sys.argv[2]) as f:
    snps = set(line.strip() for line in f)

with gzip.open(sys.argv[3], 'rt') as f:
    seen = set()
    ld = (line.split() for line in f)
    for a, b, r, d in ld:
        if a in tags and b in snps:
            tags.remove(a)
            seen.add(a)
        elif b in tags and a in snps:
            tags.remove(b)
            seen.add(b)
    for s in seen:
        print(s)
