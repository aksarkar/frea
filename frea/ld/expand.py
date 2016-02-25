"""Get all SNPs in LD with a set of tags

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import sys

data = (line.split() for line in sys.stdin)
tags, rest = set(), dict()
for row in data:
    if row[-1] == '1':
        tags.add(row[3])
        print(*row, sep='\t')
    else:
        rest[row[3]] = row

with gzip.open(sys.argv[1], 'rt') as f:
    data = (line.split() for line in sys.stdin)
    for a, b, r, d in data:
        if b in tags and a in rest:
            a, b = b, a
        if a in tags and b in rest:
            print(*rest[b], sep='\t')
