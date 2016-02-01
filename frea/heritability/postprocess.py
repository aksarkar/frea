import itertools
import operator
import glob

from ..algorithms import benjamini_hochberg

def parse_gcta(file):
    data = (line.split() for line in file)
    pve, se, p = 0, 0, 1
    for row in data:
        if row[0] == "V(G1)/Vp_L":
            pve, se = [float(x) for x in row[1:]]
        if row[0] == "Pval":
            p = float(row[1])
    return pve, se, p

def count_snps(file):
    result = 0
    for line in file:
        result += 1
    return result

def handle(file):
    pheno, cluster, *_ = file.split('.')
    with open(file) as f:
        pve, se, p = parse_gcta(f)
    with open('{}.{}.tags'.format(pheno, cluster)) as f:
        n = count_snps(f)
    return pheno, int(cluster), pve, se, p, n

if __name__ == '__main__':
    key = operator.itemgetter(4)
    result = sorted([handle(f) for f in glob.glob('*.hsq')], key=key)
    keep = sorted(benjamini_hochberg(result, key=key), key=operator.itemgetter(0, 1))
    for row in keep:
        print(*row)
