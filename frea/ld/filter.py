"""Filter SNPs in LD with specified tags

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import sys

with open(sys.argv[1]) as f:
    a = set(line.strip() for line in f)

with gzip.open(sys.argv[2], 'rt') as f:
    b = set(line.split()[3] for line in f)

snps = a | b

with gzip.open(sys.argv[3], 'rt') as f:
    data = (line.split() for line in f)
    for a_, b_, r, d in data:
        if float(r) > 0.8 and a_ in snps and b_ in snps:
            print(a_, b_, r, d)
