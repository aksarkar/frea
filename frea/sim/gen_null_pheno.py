""" Generate a null phenotype

Expects OXSTATS sample on stdin. Writes OXSTATS sample on stdout

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import sys

import numpy.random
import scipy.stats

numpy.random.seed(0)
N = scipy.stats.norm()
data = (line.split() for line in sys.stdin)
h1 = next(data)
print(*h1)
h2 = next(data)
h2[-1] = 'P'
print(*h2)
for row in data:
    row[-1] = N.rvs()
    print(*row)
