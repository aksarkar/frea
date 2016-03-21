"""Simulate phenotypes with different genetic architectures.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import gzip
import random
import pdb
import sys

import numpy

def load_annotations(file):
    data = (line.split() for line in f)
    result = collections.defaultdict(list)
    for row in data:
        result[row[-1]].append(row[:-1])
    return result

def load_dosages(file):
    data = itertools.dropwhile(lambda x: x.startswith('#'), file)
    for line in data:
        yield line.split()

def simulate(samples, dosages, annotations, causal_proportion=None, pve=None,
             num_causal=1000):
    if causal_proportion is None and pve is None:
        raise ValueError('At least one of causal_proportion and pve must be set')
    if causal_proportion is None:
        num_causal_per_annot = {k: num_causal // len(pve) for k in pve}
    else:
        num_causal_per_annot = {k: int(num_causal * v) for k, v in
                                causal_proportion}
    if pve is None:
        pve = {k: .5 / len(num_causal_per_annot) for k in num_causal_per_annot}
    N = numpy.random.normal
    causal_vars = {k: random.sample(annotations[k], num_causal_per_annot[k])
                   for k in num_causal_per_annot}
    causal_effects = {tuple(snp): N(scale=numpy.sqrt(pve[k] / num_causal_per_annot[k]))
                      for k in causal_vars for snp in causal_vars[k]}
    genetic_values = numpy.zeros(len(samples))
    for row in dosages:
        snp, dose = row[:5], row[9:]
        snp = tuple(snp)
        if snp in causal_effects:
            genetic_values += numpy.array(dose, dtype='float') * causal_effects[snp]
    noise_scale = numpy.sqrt(numpy.var(genetic_values) * (1 / sum(pve.values()) - 1))
    return genetic_values + N(scale=noise_scale, size=len(samples))

if __name__ == '__main__':
    random.seed(sys.argv[1])
    with open(sys.argv[2]) as f:
        samples = [line.split() for line in f][2:]
    with gzip.open(sys.argv[3], 'rt') as f:
        annotations = load_annotations(f)
    with gzip.open(sys.argv[4], 'rt') as f:
        dosages = load_dosages(f)
        phenotype = simulate(samples, dosages, annotations, pve={'77': .4, '66': .1})
        for p in phenotype:
            print(p)
