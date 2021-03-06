"""Algorithms needed throughout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import heapq
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

def hashjoin(seq1, seq2, key1, key2=None):
    """Yield pairs of elements s1, s2 in seq1, seq2 such that key1(s1) == key2(s2)

    This implementation performs a hash-join of two sequences. It requires that
    both keys be hashable. It also always builds a hash on seq1 since it may
    not be possible to identify which sequence is smaller.

    By default, the sequences are joined on the same key.

    seq1, seq2 - iterables
    key1, key2 - join key for each sequence

    """
    if key2 is None:
        key2 = key1
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
    n = 1
    for x in xs:
        n += 1
        new_mean = running_mean + (x - running_mean) / n
        running_variance += (x - running_mean) * (x - new_mean)
        running_mean = new_mean
    running_variance /= n
    return running_mean, running_variance

def benjamini_hochberg(data, key=operator.itemgetter(0), fdr=.05, n_tests=None):
    """Returns a generator which yields rows of data filtered by p-value to achieve
target FDR.

    This method implements the Benjamini-Hochberg stepwise procedure.

    """
    if n_tests is None:
        n_tests = len(data)
    threshold = None
    for i, row in enumerate(sorted(data, key=key)):
        p = key(row)
        if threshold is None and p > fdr * (i + 1) / n_tests:
            threshold = p
        if threshold is None or p <= threshold:
            yield row

def merge(*iterables, key):
    """Merge sorted iterables according to key.

    Backport this from Python 3.5

    https://hg.python.org/cpython/rev/f5521f5dec4a

    """
    h = []
    h_append = h.append
    _heapify = heapq.heapify
    _heappop = heapq.heappop
    _heapreplace = heapq.heapreplace
    direction = 1
    for order, it in enumerate(map(iter, iterables)):
        try:
            next = it.__next__
            value = next()
            h_append([key(value), order * direction, value, next])
        except StopIteration:
            pass
    _heapify(h)
    while len(h) > 1:
        try:
            while True:
                key_value, order, value, next = s = h[0]
                yield value
                value = next()
                s[0] = key(value)
                s[2] = value
                _heapreplace(h, s)
        except StopIteration:
            _heappop(h)
    if h:
        key_value, order, value,  next = h[0]
        yield value
        yield from next.__self__

def kwise(iterable, k):
    it = iter(iterable)
    return zip(*[it for _ in range(k)])
