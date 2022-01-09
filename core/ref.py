#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' author tanbowen '''

import os,re,sys
import pyfaidx
import numpy as np
import pandas as pd


## 获得region的GC含量
def region2gc(ref,region):

    gc_frac = []
    with pyfaidx.Fasta(ref, as_raw=True) as fa_file :
        for index,row in region.iterrows():
             subseq = fa_file[row[0]][row[1]:row[2]]
             frac_gc,frac_at = calculate_gc(subseq)
             gc_frac.append(frac_gc)

    region['gc'] = gc_frac
    region['gc'] = region['gc'].astype('float')



## 获得序列GC含量
def calculate_gc(subseq):

    total = 0
    cnt_at = subseq.count('a') + subseq.count('A') + subseq.count('t') +subseq.count('T')
    cnt_gc = subseq.count('g') + subseq.count('G') + subseq.count('c') + subseq.count('C')

    total = float(cnt_gc + cnt_at)
    if total  == 0:
        return 0.0, 0.0
    frac_gc = round(cnt_gc / total,2)
    frac_at = round(cnt_at / total,2)

    return frac_gc, frac_at


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("python3 {0} ref bed".format(sys.argv[0]))
        sys.exit(-1)

    region = region.getBedRegion(sys.argv[2])
    #region = region.loadBedRegion(sys.argv[2])
    region2gc(sys.argv[1],region)
    print(region)
