"""Draw master regulator-cofactor network

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import contextlib
import collections
import itertools
import operator
import math
import re
import sys

from .enrichments import generate_long_form

phenotype_patt = re.compile('^.*-')
tf_patt = re.compile('_[0-9]*$')

colors = {
    'AD': '#e41a1c',
    'BIP': '#377eb8',
    'CAD': '#4daf4a',
    'CD': '#984ea3',
    'RA': '#ff7f00',
    'SCZ': '#ffff33',
    'T1D': '#a65628',
    'T2D': '#f781bf',
}

chromosomes = collections.OrderedDict(
    [('chr1', 249250621),
     ('chr2', 243199373),
     ('chr3', 198022430),
     ('chr4', 191154276),
     ('chr5', 180915260),
     ('chr6', 171115067),
     ('chr7', 159138663),
     ('chrX', 155270560),
     ('chr8', 146364022),
     ('chr9', 141213431),
     ('chr10', 135534747),
     ('chr11', 135006516),
     ('chr12', 133851895),
     ('chr13', 115169878),
     ('chr14', 107349540),
     ('chr15', 102531392),
     ('chr16', 90354753),
     ('chr17', 81195210),
     ('chr18', 78077248),
     ('chr20', 63025520),
     ('chrY', 59373566),
     ('chr19', 59128983),
     ('chr22', 51304566),
     ('chr21', 48129895),
    ])

@contextlib.contextmanager
def block(name):
    print(name, '{')
    yield
    print('}')

def output_gv():
    key = operator.itemgetter(0)
    data = sorted(generate_long_form(), key=key)
    tfs = set()
    co_tfs = set()
    edges = []
    for k, g in itertools.groupby(data, key=key):
        phenotype = phenotype_patt.sub('', k).upper()
        for record in g:
            _, _, enrichment, cofactor = record
            tf = tf_patt.sub('', enrichment.motif)
            tfs.add(tf)
            co_tf = tf_patt.sub('', cofactor.motif)
            co_tfs.add(tf)
            edges.append('"{}" -- "{}" [color="{}"];'.format(tf, co_tf, colors[phenotype]))
    with block('graph {}'.format(phenotype)):
        print('start=0;')
        print('overlap=0;')
        print('splines=1;')
        print('fontname="Nimbus Sans L Regular";')
        print('node [fontsize=16,shape="plaintext"];')
        for i, tf in enumerate(tfs):
            theta = 2 * math.pi * i / len(tfs)
            r = 4
            print('"{}" [pos="{:.3f},{:.3f}",shape=box]'.format(tf, r * (1 + math.cos(theta)), r * (1 + math.sin(theta))))
        for co_tf in co_tfs:
            print('"{}";'.format(co_tf))
        for edge in edges:
            print(edge)

if __name__ == '__main__':
    output_gv()
