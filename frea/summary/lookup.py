""" Join SNP lists to 1KG reference information

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import functools
import gzip
import itertools
import operator
import math
import re
import sys

from ..algorithms import join
from ..formats import *

def get_pouyak_name(chr_, name, pos, a0, a1):
    if len(a0) <= len(a1):
        start = pos
        end = pos
        if len(a0) == len(a1):
            delta = a1
        else:  # insertion
            delta = a1[len(a0):]
    else:  # deletion
        start = pos + len(a1)
        end = pos + len(a0) - len(a1)
        delta = ""
    return '|'.join([name, chr_, str(end), str(end), delta])

def get_plink_alleles(a0, a1):
    if len(a0) > 1:
        a0 = 'I{}'.format(len(a0))
        a1 = 'D'
    elif len(a1) > 1:
        a0 = 'D'
        a1 = 'I{}'.format(len(a1))
    return a0, a1

def is_alignable(a0, a1, b0, b1):
    return (b0, b1) == (a0, a1) or (b1, b0) == (a0, a1)

def output(bed_entry, ref_entry, key=operator.itemgetter(0, -1), output_fn=output_ucsc_bed,
           score_fn=None, file=sys.stdout):
    if score_fn is None:
        score_fn = lambda x: -math.log10(x)
    chr_, score = key(bed_entry)
    name, pos, a0, a1, _ = ref_entry
    print(output_fn(chr_, score_fn(score), name, pos, a0, a1), file=file)

def lookup(parsed_input, output_fn, input_sort_key=operator.itemgetter(0),
           input_join_key=operator.itemgetter(1), ref_join_key=operator.itemgetter(1)):
    """Lookup input SNPs in 1KG

    parsed_input - iterable of parsed entries
    output_fn - function to output hits
    input_sort_key - index of chromosome column
    input_join_key - key to lookup in 1KG
    ref_join_key - key to match in 1KG
    """
    for k, g in itertools.groupby(parsed_input, key=input_sort_key):
        with gzip.open('/broad/compbio/aksarkar/data/1kg/ALL.{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(k), 'rt') as f:
            next(f)
            ref = parse(legend_format, (line.split() for line in f))
            keep = (x for x in ref if x[-1] == 'SNP')
            for pair in join(g, keep, input_join_key, ref_join_key):
                output_fn(*pair)
