#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
@author:tanbowen
"""

import os,sys,re
import numpy as np
import pandas as pd
import subprocess as sp
from multiprocessing import Pool
import matplotlib.colors as colors
#from plotnine import *
import core
from core import bamstat
import core.depth
import core.lowess
import core.ref
import core.region
import core.gc_corrent
import core.baseline
import core.hmm
import core.cnv
import core.anno


from optparse import OptionParser


def parseCommand():
    usage = "usage: ./" + sys.argv[0] + " -i sample.list -o . -s all -c 3"
    version = "%prog 1.0"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i", "--input", dest="input", help="sample list file")
    parser.add_option("-r", "--ref", dest="ref", help="genome ref file")
    parser.add_option("-t", "--target", dest="target", help="target bed file")
    parser.add_option("-b", "--baseLine", dest="baseLine", help="baseLine file")
    parser.add_option("-o", "--outdir", dest="outdir", default=".", help="output file")
    parser.add_option("-c", "--control", dest="control", help="control sample list")
    parser.add_option("-w", "--cpu-nbr", dest="cpuNbr", help="cpu nbr used for processing files, default set to 1",type="int",default=1)

    return parser.parse_args()


def main():

    (options, args) = parseCommand()
    if options.input == None:
        print("see -h for help")
        sys.exit(-1)

    if options.outdir == None:
        options.outdir ='.'
    outdir = os.path.abspath(options.outdir)

    if not os.path.exists(outdir):
        os.mkdir(outdir)


    caseSampleInfo = getSampleInfo(options.input)
    controlSampleInfo = dict()
    if options.control !=None:
        controlSampleInfo = getSampleInfo(options.control)


    # load target bed 
    targetBed = core.region.getBedRegion(options.target)
    #core.ref.region2gc(options.ref,targetBed)
    chromList = core.region.get_chrom_list(targetBed)

    saDepth = dict()
    cnDepth = dict()
    fixDepth = dict()
    baseline_ratio = pd.DataFrame()
    i = 0

    binBed = core.region.shift_window(targetBed,win=150,shift=100,min_win=100)
    core.ref.region2gc(options.ref,binBed)

    ### 获得target 区间深度 并进行GC 修正

    for sample in caseSampleInfo:

        sdepth = pd.DataFrame()
        chrJobs = []
        frames = []
        chrBinBed = []
        mainPool = Pool(6)
        for c in range(0,len(chromList)):
            chrBinBed.append(binBed[binBed["chrom"]==chromList[c]])
            chrJobs.append(mainPool.apply_async(core.depth.bin_depth, args=(caseSampleInfo[sample]["bamfile"],chrBinBed[c],20)))
        mainPool.close()
        mainPool.join()
        for b in range(0,len(chrJobs)):
            frames.append(chrJobs[b].get())
        sdepth = pd.concat(frames)

        print("mainPool test")
        sdepth = core.region.region_sort(sdepth)
        sdepth.reindex()
        saDepth[sample] = sdepth
        saDepth[sample].to_csv("{0}/{1}.raw.depth.csv".format(outdir,sample),sep="\t")


def getSampleInfo(file):
    sampleInfo = dict()

    with open(file,'r') as inf:
        for line in iter(inf):
            line = line.strip()
            info = re.split(r'\s+',line)
            if not sampleInfo.get(info[0]):
                sampleInfo.update({info[0]:{"gender":info[1],"bamfile":info[2]}})
            else:
                print ("ERROR This sample {0} appeared more than twice.\n".format(info))
                sys.exit(-1)

    return sampleInfo



if __name__ == '__main__':

    main()
