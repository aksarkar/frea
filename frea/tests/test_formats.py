import gzip
import os

import pytest

from frea.formats import *

@pytest.fixture
def set_working_dir():
    os.chdir('/broad/compbio/aksarkar/projects/roadmap/data/1kg/')
    sample_file = 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
    legend_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'
    haps_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'
    return sample_file, legend_file, haps_file

def test_parse_oxstats_haps(set_working_dir):
    sample_file, legend_file, haps_file = set_working_dir
    with open(sample_file) as f:
        samples = [line.split() for line in f][1:]
    with gzip.open(legend_file, 'rt') as f, gzip.open(haps_file, 'rt') as g:
        result = next(parse_oxstats_haps(samples, f, g))
    assert len(result[1]) == 2 * len([x for x in samples if x[2] == 'EUR'])

def test_parse_oxstats_haps_amr(set_working_dir):
    sample_file, legend_file, haps_file = set_working_dir
    with open(sample_file) as f:
        samples = [line.split() for line in f][1:]
    with gzip.open(legend_file, 'rt') as f, gzip.open(haps_file, 'rt') as g:
        result = next(parse_oxstats_haps(samples, f, g, group='AMR'))
    assert len(result[1]) == 2 * len([x for x in samples if x[2] == 'AMR'])

def test_oxstats_haplotypes(set_working_dir):
    sample_file, legend_file, haps_file = set_working_dir
    with oxstats_haplotypes(sample_file, legend_file, haps_file) as data:
        next(data)
