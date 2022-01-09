#!/usr/bin/env python3
# -*- coding:utf-8 -*-
''' author tanbowen '''

import os,sys,re
import numpy as np
import pandas as pd
from itertools import takewhile

# 读取 bed 文件，获得 target 区间


def getBedRegion(bedf):
    region = loadBedRegion(bedf)
    region = region_sort(region)
    region = region_merge(region)
    return region


def loadBedRegion(bedf):

    region = pd.DataFrame(columns=['chrom', 'start', 'end'])
    row_id = 0
    with open(bedf) as f:
        for line in iter(f):
            line = line.strip()
            if line.startswith('#'):
                continue
            infos = re.split(r'\s+', line)
            if len(infos) < 3:
                print(len(infos))
                print("the input bed fomart error")
                sys.exit(-1)

            region.loc[row_id] = [infos[0], infos[1],infos[2]]
            row_id += 1
    region[['start','end']] = region[['start','end']].astype('int')
    return region


def chr_dict():

    chr_list = [i for i in range(1,23)]
    chr_list.extend(['X','Y'])
    list2 = []
    for i in chr_list:
        list2.append("chr{0}".format(i))
    chr_list.extend(list2)
    chr_dict = dict(zip(chr_list,[1 for i in range(len(chr_list))]))

    return chr_dict

def autosomes_dict():
    chr_list = ["{0}".format(i) for i in range(1,23)]
    list2 = []
    for i in chr_list:
        list2.append("chr{0}".format(i))
    chr_list.extend(list2)
    chr_dict = dict(zip(chr_list,[1 for i in range(len(chr_list))]))
    return chr_dict


def region_sort(region):

    chrom_list = sort_chrom(region['chrom'])
    region = (region.assign(_sort_key_=chrom_list)
            .sort_values(by=['chrom', 'start', 'end'],kind='mergesort')
            .drop('_sort_key_', axis=1)
            .reset_index(drop=True))

    return region


def shift_window(region,win=80,shift=50,min_win=30):

    win_region = pd.DataFrame(columns = ['chrom','start','end'])
    m = 0

    for index,row in region.iterrows():

        for i in range(row[1],row[2],shift):
            j = i + win
            if j > row[2]:
                j = row[2]
            bin_len = j - i
            if bin_len < min_win:
                next
            win_region.loc[m] = [row[0],i,j]
            m += 1

    win_region[['start','end']] = win_region[['start','end']].astype('int')
    return win_region


def region_merge(region):

    region = region_sort(region)
    mregion = pd.DataFrame(columns=['chrom', 'start', 'end'])
    row_id = 0
    chrom = ''
    start = -1
    end = -1
    for index,row in region.iterrows():
        if index == 0:
            chrom,start,end  = row[0:3]
        else:
            if chrom == row[0]:
                if  end >= row[1]:
                    if end < row[2]:
                        end = row[2]
                else:
                    mregion.loc[row_id] = [chrom,start,end]
                    row_id += 1
                    chrom,start,end = row[0:3]
            else:
                mregion.loc[row_id] = [chrom,start,end]
                row_id +=1
                chrom,start,end = row[0:3]
    mregion.loc[row_id] = [chrom,start,end]
    mregion[['start','end']] = mregion[['start','end']].astype('int')

    return mregion


def sort_chrom(chroms):
    sorted_chrom = sorted(chroms,key=lambda chrom : sort_chrom_key(chrom))
    return sorted_chrom


def get_chrom_list(region):
    clist = list(set(region['chrom']))
    clist = sort_chrom(clist)
    return clist

def get_chrom_dict(region):
    clist = get_chrom_list(region)
    cdict = dict(zip(clist,[1*len(clist)]))
    return cdict



def is_sex_chrom(chrom):
    if chrom == 'chrX' or chrom == 'chrY' or chrom =='X' or chrom == 'Y' :
        return True
    return False

def is_autosomes(chrom):
    chrom =str(chrom).replace('chr','')
    autosomes = autosomes_dict()
    if autosomes.get("{0}".format(chrom)):
        return True
    return False



def is_chrY(chrom):
    if chrom == 'chrY' or chrom == 'Y':
        return True
    return False


def is_chrX(chrom):
    if chrom =='chrX' or chrom == 'X':
        return True
    return False



# sort_chrom_key copy from skgenome chromsort

def sort_chrom_key(label):
    """Create a sorting key from chromosome label.

    Sort by integers first, then letters or strings. The prefix "chr"
    (case-insensitive), if present, is stripped automatically for sorting.

    E.g. chr1 < chr2 < chr10 < chrX < chrY < chrM
    """
    # Strip "chr" prefix
    chrom = (label[3:] if label.lower().startswith('chr')
             else label)
    if chrom in ('X', 'Y'):
        key = (1000, chrom)
    else:
        # Separate numeric and special chromosomes
        nums = ''.join(takewhile(str.isdigit, chrom))
        chars = chrom[len(nums):]
        nums = int(nums) if nums else 0
        if not chars:
            key = (nums, '')
        elif len(chars) == 1:
            key = (2000 + nums, chars)
        else:
            key = (3000 + nums, chars)
    return key


if __name__ =='__main__':

    chrom = [1,2,3,4,5,6,7,8,9,10,'chrX','chrX','chrY','chrY','chr2',1,2,3,4,5,6,7,8,9,10,'chrX','chrX','chrY','chrY','chr2','chr3','chrY','chrX']
    for i in chrom:
        print(i,is_autosomes(i))


    if len(sys.argv) < 2:
        print("python3 {0} in.bed".format(sys.argv[0]))
        sys.exit(-1)

    region = loadBedRegion(sys.argv[1])
    chrom = region['chrom']
    get_chrom_list(region)
    region = region_sort(region)
    region = region_merge(region)
    #region['start'] = region['start'].astype('int')
    hh = region['start'].groupby(region['chrom']).mean()
    dict1 = chr_dict()
    print(hh.mean())
    print(hh.std())

