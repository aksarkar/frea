import gzip
import itertools
import os
import random

import pytest

from frea.simulation import *

@pytest.fixture
def set_working_dir():
    os.chdir('/broad/hptmp/aksarkar/dummy')
    sample_file = 'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
    legend_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'
    haps_file = 'ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz'
    hotspot_file = 'decode-chr1.txt.gz'
    return sample_file, legend_file, haps_file, hotspot_file

def test_blocks_data(set_working_dir):
    sample_file, legend_file, haps_file, hotspot_file = set_working_dir
    iterable = blocks(hotspot_file, sample_file=sample_file,
                      legend_file=legend_file, haps_file=haps_file)
    block, rate = next(iterable)
    block = list(block)
    import pdb; pdb.set_trace()
    assert len(block) == 21069

def test_blocks_final_block(set_working_dir):
    sample_file, legend_file, haps_file, hotspot_file = set_working_dir
    l = 0
    for block, rate in blocks(hotspot_file, sample_file=sample_file,
                              legend_file=legend_file, haps_file=haps_file):
        l = len(list(block))
    assert l == 10136

def test_blocks_rate(set_working_dir):
    sample_file, legend_file, haps_file, hotspot_file = set_working_dir
    block, rate = next(blocks(hotspot_file, sample_file=sample_file,
                              legend_file=legend_file, haps_file=haps_file))
    assert numpy.isclose(rate, 12.227070)

def test_sample_events(set_working_dir):
    sample_file, legend_file, haps_file, hotspot_file = set_working_dir
    y, events = sample_events(0, 100, hotspot_file, sample_file=sample_file,
                              legend_file=legend_file, haps_file=haps_file)
    assert y.shape == (100,)
