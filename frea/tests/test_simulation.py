import gzip
import itertools
import os
import random

import pytest

from frea.algorithms import kwise
from frea.formats import *
from frea.simulation import *

@pytest.fixture
def set_working_dir():
    os.chdir('/broad/compbio/aksarkar/projects/roadmap/data/1kg/')
    sample_file = 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
    legend_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'
    haps_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'
    random.seed(0)
    return sample_file, legend_file, haps_file

def test_sample_chromosome_mosaic(set_working_dir):
    sample_file, legend_file, haps_file = set_working_dir
    ancestors, variants, effects = next(sample_chromosome_mosaic(10, sample_file, legend_file, haps_file))
    assert len(ancestors) == 20

def test_reconstruct_chromosome_mosaic(set_working_dir):
    sample_file, legend_file, haps_file = set_working_dir
    with oxstats_haplotypes(sample_file, legend_file, haps_file) as data:
        expected_haps = [h[:20] for l, h in itertools.islice(data, 0, 10)]
        expected_genotypes = [[sum(pair) for pair in kwise(h, 2)] for h in expected_haps]
    ancestors = list(range(20))
    reconstuction = reconstruct_chromosome_mosaic(ancestors, sample_file, legend_file, haps_file)
    genotypes = list(itertools.islice(reconstuction, 0, 10))
    assert genotypes == expected_genotypes
