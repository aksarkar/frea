"""Summarize motif enrichments.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import glob
import itertools
import operator
import sys

Enrichment = collections.namedtuple('Enrichment', ['cluster', 'motif', 'odds_ratio', 'p'])
Cofactor = collections.namedtuple('Cofactor', ['motif', 'score'])

def parse_motif_matches(filename):
    """Return a dict mapping master regulators (string) to cofactors (Cofactor)"""
    result = collections.defaultdict(list)
    with open(filename) as f:
        data = (line.split() for line in f)
        key = operator.itemgetter(0, 1)
        for k, g in itertools.groupby(sorted(data, key=key), key=key):
            master_regulator, cofactor = k
            result[master_regulator].append(Cofactor(cofactor, max(x[6] for x in g)))
    return result

def parse_motif_enrichments(filename):
    """Yield Enrichment namedtuples

    Take the motif with best p-value as cluster representative. This assumes
    input data is already appropriately sorted.

    """
    with open(filename) as f:
        data = (line.split() for line in f)
        for k, g in itertools.groupby(data, key=operator.itemgetter(0)):
            yield Enrichment(*next(iter(g)))

def generate_long_form():
    """Generate long-form tuples for use in plotting"""
    for filename in glob.glob('*.summary'):
        phenotype, cluster, *_ = filename.split('.')
        cofactors = parse_motif_matches('{}.{}.matches'.format(phenotype, cluster))
        for enrichment in parse_motif_enrichments(filename):
            for cofactor in cofactors[enrichment.motif]:
                yield phenotype, cluster, enrichment, cofactor

def output_long_form(output_key):
    """Output long-form tuples for use in plotting"""
    key = operator.itemgetter(0, 1)
    data = sorted(generate_long_form(), key=key)
    for k, g in itertools.groupby(data, key=key):
        for record in g:
            phenotype, cluster, enrichment, cofactor = record
            print(phenotype, cluster, *output_key(enrichment, cofactor))
    
def output_enrichments():
    """Output (phenotype, cluster, motif, odds, p) records"""
    output_long_form(lambda enrichment, _: enrichment)

def output_cofactors():
    """Output (phenotype, cluster, master regulator, cofactor, score) records"""
    def output_key(enrichment, cofactor):
        return enrichment.motif, cofactor.motif, cofactor.score
    output_long_form(output_key)
