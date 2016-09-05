""" Shard a BED file by chromosome for pruning

Author: Abhishek Sarkar

"""
import argparse
import gzip
import itertools
import operator
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--ref', action='store_true',
                    help='Shard reference LD information')
parser.add_argument('file', help='File to split by chromosome')
args = parser.parse_args()

if args.ref:
    key = lambda x: x[0].split('|')[1]
else:
    key=operator.itemgetter(0)

with gzip.open(args.file, 'rt') as f:
    data = (line.split() for line in f)
    for k, g in itertools.groupby(data, key=key):
        with gzip.open('{}.txt.gz'.format(k), 'wt') as out:
            for row in g:
                print(*row, file=out)
