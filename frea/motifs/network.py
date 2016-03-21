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
    'AD':	'#1b9e7780',
    'BIP':	'#d95f0280',
    'CAD':	'#7570b380',
    'CD':	'#e7298a80',
    'RA':	'#66a61e80',
    'SCZ':	'#e6ab0280',
    'T1D':	'#a6761d80',
    'T2D':	'#66666680',
}

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
            if tf != co_tf:
                edges.append('"{}" -- "{}" [color="{}"];'.format(tf, co_tf, colors[phenotype]))
    with block('graph {}'.format(phenotype)):
        print('start=0;')
        print('overlap="prism";')
        print('size="8,5"')
        print('splines="curved";')
        print('outputorder="edgesfirst"')
        print('node [fontname="Helvetica",fontsize=24,shape="plaintext"];')
        for i, tf in enumerate(tfs):
            print('"{}" [label=<<b>{}</b>>,fontsize=48]'.format(tf, tf))
        for co_tf in co_tfs:
            if co_tf not in tfs:
                print('"{}";'.format(co_tf))
        for edge in edges:
            print(edge)

if __name__ == '__main__':
    output_gv()
