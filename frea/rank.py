"""Estimate robustness of GWAS p-value ranks to sample size.

Usage: python -m frea.rank SUMMARY_FILE

Compute correlation between held-out cohorts and the remaining in large
meta-analyses where per-cohort summary statistics are available.

- Stahl, E. A., Raychaudhuri, S., Remmers, E. F., Xie, G., Eyre, S., Thomson,
  B. P., Li, Y., Kurreeman, F. A. S., Zhernakova, A., Hinks, A., Guiducci, C.,
  Chen, R., Alfredsson, L., Amos, C. I., Ardlie, K. G., Barton, A., Bowes, J.,
  Brouwer, E., Burtt, N. P., Catanese, J. J., Coblyn, J., Coenen, M. J. H.,
  Costenbader, K. H., Criswell, L. A., Crusius, J. B. A., Cui, J., de Bakker,
  P. I. W., De Jager, P. L., Ding, B., Emery, P., Flynn, E., Harrison, P.,
  Hocking, L. J., Huizinga, T. W. J., Kastner, D. L., Ke, X., Lee, A. T., Liu,
  X., Martin, P., Morgan, A. W., Padyukov, L., Posthumus, M. D., Radstake, T.
  R. D. J., Reid, D. M., Seielstad, M., Seldin, M. F., Shadick, N. A., Steer,
  S., Tak, P. P., Thomson, W., van der Helm-van Mil, A. H. M., van der
  Horst-Bruinsma, I. E., van der Schoot, C. E., van Riel, P. L. C. M.,
  Weinblatt, M. E., Wilson, A. G., Wolbink, G. J., Wordsworth, B. P., Wijmenga,
  C., Karlson, E. W., Toes, R. E. M., de Vries, N., Begovich, A. B.,
  Worthington, J., Siminovitch, K. A., Gregersen, P. K., Klareskog, L., Plenge,
  R. M., .. (2010). Genome-wide association study meta-analysis identifies
  seven new rheumatoid arthritis risk loci. Nat Genet, 42(6), 508–514.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import csv
import gzip
import math
import itertools
import operator
import sys

def pearson(cohort):
    """Coroutine implementation of Pearson correlation

    https://groups.google.com/d/msg/sci.stat.math/HMrw3pWSvQs/-QP3BehUV1MJ

    """
    n = 0
    mean_x = 0
    mean_y = 0
    sum_x_sq = 0
    sum_y_sq = 0
    sum_xy = 0
    def output():
        print(n, cohort, sum_xy / math.sqrt(sum_x_sq * sum_y_sq))
    cutoffs = set(x * y for x, y in itertools.product([1, 2, 3, 4, 5], [10 ** i for i in range(3, 7)]))
    try:
        while True:
            x, y = yield
            n += 1
            if n in cutoffs:
                output()
            dev_x = x - mean_x
            dev_y = y - mean_y
            mean_x += dev_x / n
            mean_y += dev_y / n
            sum_x_sq += dev_x * (x - mean_x)
            sum_y_sq += dev_y * (y - mean_y)
            sum_xy += dev_x * (y - mean_y)
    except GeneratorExit:
        output()

def hold_out_cor(zscores, cohorts):
    targets = [pearson(cohort) for cohort in cohorts] + [pearson("Overall")]
    for t in targets:
        next(t)
    for row in zscores:
        sample_weighted = sum(x * w for x, w in zip(row, weights))
        for z, w, t in zip(row, weights, targets):
            t.send((z, (sample_weighted - z * w) / (total_weight - w)))
        targets[-1].send((row[-1], sample_weighted / total_weight))
    for t in targets:
        t.close()

cohorts = ["WTCCC", "NARAC1", "NARAC2", "EIRA", "CANADA", "BRASS"]
samples = [(1525, 10608),  
           (867, 1041),    
           (902, 4510),    
           (1173, 1089),   
           (489, 1472),    
           (483, 1449)]    
total_samples = sum(sum(cohort) for cohort in samples)
weights = [math.sqrt(sum(cohort) / total_samples) for cohort in samples]
total_weight = sum(weights)

with gzip.open(sys.argv[1], "rt") as f:
    data = (line.split() for line in f)
    next(data)  # skip header
    filter_chr6 = (row for row in data if row[1] != '6')
    projection = (row[11:17] + row[-4:-3] for row in filter_chr6)
    filter_missing = (row for row in projection if 'NA' not in row)
    zscores = sorted(([float(x) for x in row] for row in filter_missing), key=operator.itemgetter(-1))
    hold_out_cor(zscores, cohorts)
