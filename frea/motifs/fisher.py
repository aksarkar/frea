"""Fisher's exact test for motif enrichments

Usage: python fisher.py FG_COUNTS BG_COUNTS FG_TOTAL BG_TOTAL FDR

FG_COUNTS and BG_COUNTS are space-separated (motif, count) pairs, one per
line. FG_TOTAL and BG_TOTAL are integers.

Prints motif, odds ratio, p-value, and '*' if significant at specified FDR on
stdout

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import operator
import sys
import scipy.stats

def load(filename):
    with open(filename) as f:
        return {k: int(v) for k, v in (line.split() for line in f)}

def count(filename):
    result = 0
    with open(filename) as f:
        for line in f:
            result += 1
    return result

if __name__ == "__main__":
    fg_motif_counts, bg_motif_counts = [load(x) for x in sys.argv[1:3]]
    fg_total_count = count(sys.argv[3])
    bg_total_count = count(sys.argv[4])
    fdr = float(sys.argv[5])
    T = scipy.stats.fisher_exact
    enrichment = [(k,) + T([[fg_motif_counts[k], fg_total_count - fg_motif_counts[k]],
                            [bg_motif_counts[k], bg_total_count - bg_motif_counts[k]]],
                           alternative='greater')
                  for k in fg_motif_counts]
    for i, (motif, odds, p, *_) in enumerate(sorted(enrichment, key=operator.itemgetter(2))):
        if p < fdr * (i + 1) / len(enrichment):
            print(motif, odds, p, '*')
        else:
            print(motif, odds, p)
