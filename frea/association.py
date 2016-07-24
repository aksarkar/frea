import sys

from frea.simulation import *
from frea.formats import *

with oxstats_genotypes(sample_file=sys.argv[1], gen_file=sys.argv[2]) as (_, _, samples, data):
    y = numpy.array([float(row[-1]) for row in samples])
    y -= y.mean()
    y /= y.std()
    for k, g in itertools.groupby(enumerate(data), key=lambda x: x[0] // 1e3):
        block = [x for i, x in g]
        x = numpy.array([oxstats_gen_to_dosage(x_j[5:]) for x_j in block]).T
        x -= x.mean(axis=0)
        b, se, logp = compute_marginal_stats(x, y)
        for (_, name, pos, a0, a1, *_), b_j, s, p in zip(block, b, se, logp):
            print(name, pos, a0, a1, b_j, s, p)
