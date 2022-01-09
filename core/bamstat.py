#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''author: tanbowen '''
import os,sys,re
import numpy as np
import pandas as pd
import pysam
from . import region

def bam_index_stat(bamf):

    indexStat = False

    if os.path.isfile(bamf + '.bai') or os.path.isfile(bamf.replace('.bam','.bai')):
        indexStat = True

    return indexStat


def getPloidy(data,gender):
    ploidy = dict()
    autosomes_dict = region.autosomes_dict()
    chr_list = region.get_chrom_list(data)
    for c in chr_list:
        if autosomes_dict.get(c):
            ploidy[c]  =2
        elif c.replace('chr','') == 'Y':
            if gender == 'M':
                ploidy[c] = 1
            else:
                ploidy[c] = 0
        elif c.replace('chr','') == 'X':
            if gender == 'M':
                ploidy[c] = 1
            else:
                ploidy[c] = 2
        else:
            ploidy[c] = 1

    return ploidy




def get_ploidy(region_depth,chr_list):

    gender = guess_gender(region_depth)
    autosomes_dict = region.autosomes_dict()
    ploidy = dict()

    for c in chr_list:

        if autosomes_dict.get(c):
            ploidy[c] = 2

        elif c.replace('chr','') == 'Y':
            if gender == 'M':
                ploidy[c] = 1
            else:
                ploidy[c] = 0
        elif c.replace('chr','') == 'X':

            if gender == 'M':
                ploidy[c] = 1
            else:
                ploidy[c] = 2
        else:
            ploidy[c] = 1

    return ploidy


def guess_gender(regin_depth):

    chr_dict = region.chr_dict()
    chr_list = list(chr_dict.keys())
    chr_depth = region_detph[regin_depth.chrom.isin(chr_list)]
    chrom_mean_depth = chrom_mean_bin_depth(chr_depth)
    clists = region.get_chrom_list(chrom_mean_depth)
    cdict = dict(zip(clist,[1 for i in range(len(clists))]))

    chrX_depth = chrY_depth = 0
    chrX_tag = chrY_tag = 1
    gender = 'U'

    if cdict.get('chrX'):
        chrX_depth = chrom_mean_depth['chrX']
    elif cdict.get('X'):
        chrX_depth = chrom_mean_depth['X']
    else:
        chrX_tag = 0

    if cdict.get('chrY'):
        chrY_depth = chrom_mean_depth['chrY']
    elif cdict.get('Y'):
        chrY_depth = chrom_mean_depth['Y']
    else:
        chrY_tag = 0

    if chrX_tag and chrY_tag :
        ratio = chrX_depth /float(chrY_depth + 0.00001)
        if ratio < 0.2 :
            gender = 'F'
        else:
            gender = 'M'
    elif chrX_tag:
        chrX_z_score = (chrX_depth - chrom_mean_depth.mean()) / chrom_mean_depth.std()
        if chrX_z_score < -3 :
            gender = 'M'
        else:
            gender = 'F'
    elif chrY_tag:
        ratio = chrY_depth / (chrom_mean_depth.median() + 0.00001)
        if raito < 0.2 :
            gender = 'F'
        else:
            gender = 'M'

    return gender


def chrom_mean_bin_depth(region_depth):
    chrom_mean_depth = region_depth['depth'].groupby(region_depth['chrom']).mean()
    return chrom_mean_depth


def chrom_median_bin_depth(region_depth):
    chrom_median_depth = region_depth['depth'].groupby(region_depth['chrom']).mean()
    return chrom_median_depth


def region_coverage(bamf,region,min_mapq=20):

    bamfile = pysam.AlignmentFile(bamf, 'rb')
    region_depth = region.copy(deep=True)
    cov =[]
    for index,row in region.iterrows():
        depth = bamfile.count_coverage(row[0],row[1],row[2],quailty_threshold = min_mapq)
        cov.append(depth)
    region_depth['depth'] = cov
    region_depth['depth'] = region_depth['depth'].astype('float')
    bamfile.close()

    return region_depth


def depth_count(bamf,chrom,start,end,min_mapq=20):
    depth = bamf.count_coverage(chrom,start,end,quality_threshold = min_mapq)
    rlen = end - start
    allbase = sum(depth[0])
    mean_depth = round(allbase/rlen,2)
    return mean_depth


if __name__ =='__main__':

    if len(sys.argv) < 5:
        print("python3 {0} inbam chrom start end".format(sys.argv[0]))
        sys.exit(-1)

    bamfile = pysam.AlignmentFile(sys.argv[1], 'rb')
    mean_depth = depth_count(bamfile,sys.argv[2],int(sys.argv[3]),int(sys.argv[4]))
    print(mean_depth)


