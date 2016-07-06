"""Entry points for ImpG pipeline

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import gzip
import operator
import sys

from .algorithms import hashjoin
from .formats import parse, legend_format

def oxstats_legend_to_impg_map(filename, min_maf=0.01):
    """Convert OXSTATS legend file to ImpG map file

    Write one file per chromosome, then use GenMaps (included with ImpG) to
    produce files for each imputation window.

    ImpG currently only supports SNPs, so we need to filter on the information
    provided in the reference. We additionally impose a MAF filter to reduce
    the computational burden.

    """
    with gzip.open(filename, 'rt') as f:
        next(f)
        data = (line.split() for line in f)
        print("snp", "pos", "ref", "alt", sep='\t')
        for row in data:
            if row[4] == "SNP" and float(row[13]) >= .01:
                print(*row[:4], sep='\t')

def oxstats_haps_to_impg_haps(haplotypes, legend):
    """Convert OXSTATS haplotype file to ImpG format

    Write one file per chromosome, then use split_impg_haps to produce files
    for each imputation window

    """
    with gzip.open(haplotypes, 'rt') as f, gzip.open(legend, 'rt') as g:
        haps = (line.split() for line in f)
        info = (line.split() for line in g)
        next(info)
        for s, h in zip(info, haps):
            print(s[0], ''.join(chr(ord(x) + 1) for x in h))

def split_impg_haps(haplotypes, map_file, sample_file):
    """Split ImpG formatted haplotypes for one imputation window

    Imputation windows are defined by the map_file. Hash-join on rsids to get
    the correct subset of reference haplotypes.

    Assume haplotypes are gzip'ed (this is done externally in a shell
    pipeline).

    """
    with open(sample_file) as f:
        data = (line.split() for line in f)
        print('SNP', ' '.join('{} {}'.format(row[0], row[0]) for row in data))
    with gzip.open(haplotypes, 'rt') as f, open(map_file) as g:
        haps = (line.split() for line in f)
        snps = (line.split() for line in g)
        next(snps)
        for s, h in hashjoin(snps, haps, key1=operator.itemgetter(0)):
            print(h)
