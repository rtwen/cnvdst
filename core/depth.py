#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''tanbowen'''

import os,sys,re
import numpy as np
import pandas as pd
import pysam
from . import region


def bin_depth(bamf,sregion,mapq=20):
    #sregion = region.shift_window(bedregion,binSize,shiftSize,minSize)
    print("shift window done")
    bin_depth = region_coverage(bamf,sregion,mapq)
    return bin_depth


def chrom_mean_bin_depth(region_depth):
    chrom_mean_depth = region_depth['depth'].groupby(region_depth['chrom']).mean()
    return chrom_mean_depth


def chrom_median_bin_depth(region_depth):
    chrom_median_depth = region_depth['depth'].groupby(region_depth['chrom']).mean()
    return chrom_median_depth


def region_coverage(bamf,tregion,min_mapq=20):

    bamfile = pysam.AlignmentFile(bamf, 'rb')
    region_depth = tregion.copy(deep=True)
    cov =[]
    for index,row in tregion.iterrows():
        depth = depth_count(bamfile,row[0],row[1],row[2],min_mapq)
        cov.append(depth)
    region_depth['depth'] = cov
    region_depth['depth'] = region_depth['depth'].astype('float')
    #bamfile.close()
    return region_depth


def depth_count(bamf,chrom,start,end,min_mapq=20):
    depth = bamf.count_coverage(chrom,start,end,quality_threshold = min_mapq)
    rlen = end - start
    allbase = sum(depth[0] + depth[1] + depth[2] + depth[3])
    mean_depth = round(allbase/rlen,2)
   # print(mean_depth)
    return mean_depth


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("python3 {0} in.bed".format(sys.argv[0]))
        sys.exit(-1)

