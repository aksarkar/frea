"""Convert OXSTATS legend file to ImpG map

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import sys

with gzip.open(sys.argv[1], 'rt') as f:
    next(f)
    data = (line.split() for line in f)
    print("snp", "pos", "ref", "alt", sep='\t')
    for row in data:
        if row[4] == "SNP" and float(row[13]) >= .01:
            print(*row[:4], sep='\t')
