"""Entry points for ImpG pipeline

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import argparse
import gzip
import operator
import sys

from .algorithms import hashjoin
from .formats import oxstats_haplotypes

def oxstats_legend_to_impg_map():
    """Convert OXSTATS legend file to ImpG map file

    Write one file per chromosome, then use GenMaps (included with ImpG) to
    produce files for each imputation window.

    ImpG currently only supports SNPs, so we need to filter on the information
    provided in the reference. We additionally impose a MAF filter to reduce
    the computational burden.

    """
    parser = argparse.ArgumentParser(description='Convert OXSTATS legend to ImpG map file')
    parser.add_argument('file', help='OXSTATS legend file')
    parser.add_argument('--min-maf', type=float, default=0.01, help='Minimum MAF in EUR samples')
    args = parser.parse_args()
    with gzip.open(args.file, 'rt') as f:
        next(f)
        data = (line.split() for line in f)
        print("snp", "pos", "ref", "alt", sep='\t')
        for row in data:
            if row[4] == "SNP" and float(row[13]) >= args.min_maf:
                print(*row[:4], sep='\t')

def oxstats_haps_to_impg_haps():
    """Convert OXSTATS haplotype file to ImpG format

    Write one file per chromosome, then use split_impg_haps to produce files
    for each imputation window

    """
    parser = argparse.ArgumentParser(description='Convert OXSTATS haplotypes to ImpG haps file')
    parser.add_argument('haps_file', help='OXSTATS haplotypes file')
    parser.add_argument('legend_file', help='OXSTATS legend file')
    parser.add_argument('sample_file', help='OXSTATS sample file')
    args = parser.parse_args()
    with oxstats_haplotypes(args.sample_file, args.legend_file, args.haps_file, group='EUR') as data:
        for s, h in data:
            print(s[0], ''.join(str(1 + x) for x in h))

def split_impg_haps():
    """Split ImpG formatted haplotypes for one imputation window

    Imputation windows are defined by the map_file. Hash-join on rsids to get
    the correct subset of reference haplotypes.

    Assume haplotypes are gzip'ed (this is done externally in a shell
    pipeline).

    """
    parser = argparse.ArgumentParser(description='Split ImpG haplotypes into windows with overlapping buffers')
    parser.add_argument('haps_file', help='OXSTATS haplotypes file')
    parser.add_argument('map_file', help='ImpG map file')
    parser.add_argument('sample_file', help='OXSTATS sample file')
    args = parser.parse_args()
    with open(args.sample_file) as f:
        data = (line.split() for line in f)
        print('SNP', ' '.join('{} {}'.format(row[0], row[0]) for row in data))
    with gzip.open(args.haps_file, 'rt') as f, open(args.map_file) as g:
        haps = (line.split() for line in f)
        snps = (line.split() for line in g)
        next(snps)
        for s, h in hashjoin(snps, haps, key1=operator.itemgetter(0)):
            print(h)
