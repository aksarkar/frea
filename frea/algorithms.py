"""Algorithms needed throughout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import operator

def join(seq1, seq2, key1, key2=None):
    """Yield pairs of elements s1, s2 in seq1, seq2 such that key1(s1) == key2(s2).

    This implementation performs a sort-merge join of two sequences. It
    requires the keys for both sequences to be comparable and for each sequence
    to be sorted on its key.

    seq1, seq2 - iterables
    key1 - Key for sequence 1
    key2 - Key for sequence 2. If None, assumed to be equal to key1

    """
    if key2 is None:
        key2 = key1
    seq1_buckets = itertools.groupby(seq1, key1)
    seq2_buckets = itertools.groupby(seq2, key2)
    k1, g1 = next(seq1_buckets)
    k2, g2 = next(seq2_buckets)
    while True:
        if k1 == k2:
            for pair in itertools.product(g1, g2):
                yield pair
            k1, g1 = next(seq1_buckets)
            k2, g2 = next(seq2_buckets)
        elif k1 < k2:
            k1, g1 = next(seq1_buckets)
        else:
            k2, g2 = next(seq2_buckets)

def hashjoin(seq1, seq2, key1=operator.itemgetter(3), key2=operator.itemgetter(1)):
    """Yield pairs of elements s1, s2 in seq1, seq2 such that key1(s1) == key2(s2)

    This implementation performs a hash-join of two sequences. It requires that
    both keys be hashable. It also always builds a hash on seq1 since it may
    not be possible to identify which sequence is smaller.

    seq1, seq2 - iterables
    key1, key2 - join key for each sequence

    """
    lookup = {key1(s): s for s in seq1}
    for s in seq2:
        if key2(s) in lookup:
            yield lookup[key2(s)], s

def moments(xs):
    """Return the mean and variance of xs.

    This implementation performs a single pass through the elements of xs,
    which means it works for generators and also has favorable numerical
    stability.

    The algorithm is given in TAOCP, Vol. 2.

    """
    running_mean = next(iter(xs))
    running_variance = 0
    for i, x in enumerate(xs):
        new_mean = running_mean + (x - running_mean) / (i + 2)
        running_variance += (x - running_mean) * (x - new_mean)
        running_mean = new_mean
    running_variance /= ntrials
    return running_mean, running_variance

def benjamini_hochberg(data, key=operator.itemgetter(0), fdr=.05, n_tests=None):
    """Returns a generator which yields rows of data filtered by p-value to achieve
target FDR.

    This method implements the Benjamini-Hochberg stepwise procedure.

    """
    if n_tests is None:
        n_tests = len(data)
    threshold = None
    for i, row in enumerate(data):
        p = key(row)
        if threshold is None and p > fdr * (i + 1) / n_tests:
            threshold = p
        if threshold is None or p <= threshold:
            yield row
