""" Merge simulation replicates for easy plotting

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import operator
import sys

types = [int, str, str, str, float, float]

def parse(filename):
    with open(filename) as f:
        for line in f:
            yield [g(x) for g, x in zip(types, line.split(','))]

def advance(replicates):
    while True:
        yield [next(r) for r in replicates]

def quantile(p, sequence):
    return sequence[int(p * len(sequence))]

if __name__ == '__main__':
    replicates = [parse(f) for f in sys.argv[1:]]
    for items in advance(replicates):
        total, _, celltype, feature, *_ = items[0]
        expected = sum([item[-1] for item in items]) / len(items)
        ordered = sorted([item[4] for item in items])
        quantiles = quantile(.05, ordered), quantile(.5, ordered), quantile(.95, ordered)
        for phenotype, q in zip(["lower", "median", "upper"], quantiles):
            print(total, phenotype, celltype, feature, q, expected, sep=',')
