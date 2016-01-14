import collections
import itertools
import operator
import sys

seed = sys.argv[1]
data = (line.split() for line in sys.stdin)
for k, g in itertools.groupby(data, key=operator.itemgetter(8)):
    offsets = collections.Counter(int(x[1]) - int(x[6]) for x in g)
    for o, count in offsets.items():
        print(seed, k, o, count)
