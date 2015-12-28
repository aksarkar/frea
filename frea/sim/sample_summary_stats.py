""" Sample null z-scores

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import gzip 
import itertools
import operator
import sys

import numpy.linalg
import numpy.random
import scipy.stats

N = scipy.stats.norm.rvs

def sample_conditional(g, covariance, g_buff, z_buff, ntrials, ridge_lambda):
    new_covariance = numpy.dot(g_buff, g) / g.shape[0]
    schur_complement = new_covariance.dot(numpy.linalg.pinv(covariance))
    conditional_mean = schur_complement.dot(z_buff)
    conditional_variance = 1 - schur_complement.dot(new_covariance)
    assert conditional_variance > 0
    new_covariance = numpy.vstack((numpy.hstack((covariance, new_covariance[:,numpy.newaxis])),
                                   numpy.append(new_covariance, 1 + ridge_lambda)))
    return new_covariance, N(loc=conditional_mean, scale=conditional_variance, size=ntrials)

def sample(data, block_size=None, ntrials=1000, ridge_lambda=1e-4):
    name, g = next(standardized)
    if block_size is None:
        block_size = len(g)
    z = N(size=ntrials)
    yield name, z
    g_buff = [g]
    z_buff = [z]
    covariance = (1 + ridge_lambda) * numpy.eye(1)
    for name, g in standardized:
        if len(g_buff) >= block_size:
            g_buff.pop(0)
            z_buff.pop(0)
            covariance = covariance[1:,1:]
        covariance, z = sample_conditional(g, covariance, g_buff, z_buff, ntrials, ridge_lambda)
        g_buff.append(g)
        z_buff.append(z)
        yield name, z

if __name__ == '__main__':
    numpy.random.seed(0)
    data = (line.split() for line in sys.stdin)
    snps = {row[0]: row[1] for row in data}
    with gzip.open(sys.argv[1], 'rt') as f:
        data = (line.split('\t') for line in f)
        parsed = (('|'.join(row[:5]), [int(x) for x in row[5].split(';')])
                  for row in data)
        keep = (row for row in parsed if row[0] in snps)
        filter_monomorphic = ((name, haps) for name, haps in keep 
                              if .01 < haps.count(0) / len(haps) < .99)
        genotypes = ((name, numpy.sum(numpy.array(haps).reshape(-1, 2), axis=1))
                     for name, haps in filter_monomorphic)
        standardized = ((name, (g - numpy.mean(g)) / numpy.std(g)) for name, g in genotypes)
        for snpid, x in sample(data):
            print(snpid, *x)
