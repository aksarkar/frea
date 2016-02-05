"""Parsers for common formats

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

ucsc_bed_format = [str, int, int, str, float]
impg_format = [str, int, str, str, float, float]
legend_format = [str, int, str, str, str]
plink_format = [lambda x: 'chr{}'.format(x), str, int, int, str, str]

def get_pouyak_name(chr_, name, pos, a0, a1):
    if len(a0) <= len(a1):
        start = pos
        end = pos
        if len(a0) == len(a1):
            delta = a1
        else:  # insertion
            delta = a1[len(a0):]
    else:  # deletion
        start = pos + len(a1)
        end = pos + len(a0) - len(a1)
        delta = ""
    return '|'.join([name, chr_, str(end), str(end), delta])

def get_plink_alleles(a0, a1):
    if len(a0) > 1:
        a0 = 'I{}'.format(len(a0))
        a1 = 'D'
    elif len(a1) > 1:
        a0 = 'D'
        a1 = 'I{}'.format(len(a1))
    return a0, a1

def parse(types, iterable):
    for row in iterable:
        yield [f(x) for f, x in zip(types, row)]

def output_ucsc_bed(chr_, score, name, pos, a0, a1):
    return '\t'.join([chr_, str(pos - 1), str(pos + len(a0) - 1), get_pouyak_name(chr_, name, pos, a0, a1), str(score)])

def output_impg_zscores(chr_, score, name, pos, a0, a1):
    return ' '.join([name, str(pos), a0, a1, str(score)])

def output_plink_bim(chr_, score, name, pos, a0, a1):
    return '\t'.join([chr_, get_pouyak_name(chr_, name, pos, a0, a1), "0", str(pos), a0, a1])

_isf = scipy.stats.chi2(1).isf

def zscore(p, odds):
    z = math.sqrt(_isf(float(p)))
    if float(odds) < 1:
        z *= -1
    return z

def logp(p):
    return -math.log10(float(p))

def chromosome(n):
    return 'chr{}'.format(n)

def summary_stats(filename, input_key, skip_header=True):
    """Parse summary statistics

    input_key - function that returns (chromosome, position, id, z)
    """
    autosomes = set('chr{}'.format(x) for x in range(1, 23))
    with gzip.open(filename, 'rt') as f:
        if skip_header:
            next(f)
        data = (line.split() for line in f)
        for row in data:
            parsed = input_key(row)
            if parsed[0] in autosomes:
                yield parsed

