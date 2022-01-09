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
    parser.add_option("-e", "--refgene", dest="refgene", help="ucsc refgene file longest transcript file",type="str",default="/db/GATK/longest_transcript.sort.txt.gz")
    parser.add_option("-m", "--sample", dest="sample", help="which sample need to analysis,default all",default="all")
    parser.add_option("-s", "--steps",dest="stepTag",
    help="select step to run, f for fqFilter, b for bwaMem,and bam merge, m for bamMarkdup, r for bamRealignment and bamRecalibration, v for snv, j for jointVariant, t for TMB, c for cnvkit, i for MSI, q for bamQC,s for manta ,u for mergeRsult  default set to all",default="all")

    return parser.parse_args()


def main():

    (options, args) = parseCommand()
    if options.input == None:
        print("see -h for help")
        sys.exit(-1)

    if options.outdir == None:
        options.outdir ='.'
    outdir = os.path.abspath(options.outdir)


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

    binBed = core.region.shift_window(targetBed,win=100,shift=80,min_win=60)
    core.ref.region2gc(options.ref,binBed)

    ### 获得target 区间深度 并进行GC 修正

    for sample in caseSampleInfo:

        #sdepth = core.depth.bin_depth(caseSampleInfo[sample]["bamfile"],binBed,mapq=20)
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

        sdepth = core.region.region_sort(sdepth)
        sdepth.reindex()
        gc_loe = core.lowess.gc_depth_lowess(sdepth)
        core.gc_corrent.gc_corrent(gc_loe,sdepth,100)
        core.baseline.getRatio(sdepth,'M')
        saDepth[sample] = sdepth


    for sample in controlSampleInfo:

        mainPool = Pool(6)
        chrJobs = []
        frames = []
        chrBinBed = []
        sdepth = pd.DataFrame()

        for c in range(0,len(chromList)):
            chrBinBed.append(binBed[binBed["chrom"]==chromList[c]])
            chrJobs.append(mainPool.apply_async(core.depth.bin_depth, args=(controlSampleInfo[sample]["bamfile"],chrBinBed[c],20)))
        mainPool.close()
        mainPool.join()

        for b in range(0,len(chrJobs)):
            frames.append(chrJobs[b].get())
        sdepth = pd.concat(frames)

        gc_loe = core.lowess.gc_depth_lowess(sdepth)
        core.gc_corrent.gc_corrent(gc_loe,sdepth,100)
        #core.baseline.getControlRatio(sdepth,'M')
        core.baseline.getControlRatio(sdepth,'M')
        sdepth = core.region.region_sort(sdepth)
        sdepth.reindex()
        cnDepth[sample] = sdepth
        cnDepth[sample].to_csv("{0}.csv".format(sample),sep="\t")


        if i == 0:
            baseline_ratio["chrom"] = sdepth["chrom"]
            baseline_ratio["start"] = sdepth["start"]
            baseline_ratio["end"] = sdepth["end"]
            baseline_ratio["gc"] = sdepth["gc"]
        baseline_ratio[sample] = sdepth["ratio"]
        i += 1

    core.baseline.getBaseline(baseline_ratio)
    finalbaseline = core.baseline.baselineFilter(baseline_ratio)
    finalbaseline.to_csv("baseline.csv",sep="\t")


    for sample in caseSampleInfo:

        mainPool = Pool(6)
        chromFixDepth = []
        cnvJobs = []
        cn = []
        fdepth = core.baseline.baseline_fix(saDepth[sample],finalbaseline)
        fdepth.to_csv("{0}.depth.csv".format(sample),sep="\t")
        fixDepth[sample] = fdepth
        cns = core.hmm.fit_hmm(fixDepth[sample],0.1,0.01)
        fixDepth[sample]["cn"] = cns
        kk = 0

        """
        for chrom in core.region.get_chrom_list(fixDepth[sample]):

            cfd = fdepth[fdepth["chrom"] == chrom]
            cfd.to_csv("{0}.csv".format(chrom))
            chromFixDepth.append(fdepth[fdepth["chrom"] == chrom])
            #np.random.seed(1)
            cnvJobs.append(mainPool.apply_async(core.hmm.fit_hmm, args=(chromFixDepth[kk],0.1,0.01)))
            #rst = core.hmm.fit_hmm(cfd,0.1,0.01)
            kk+=1
            #cn.extend(rst)

        mainPool.close()
        mainPool.join()

        for j in range(0,len(cnvJobs)):
            x = cnvJobs[j].get()
            cn.extend(x)

        fixDepth[sample]["cn"] =cn
        """
        rst_cnvs = core.cnv.mergeCnv(fixDepth[sample],'M')
        core.anno.anno(rst_cnvs,options.refgene)
        genes = rst_cnvs["gene"]
        rst_cnvs.drop('gene',axis=1,inplace = True)
        rst_cnvs.insert(4,'gene',genes)
        fixDepth[sample].to_csv("{0}.csv".format(sample),sep="\t")
        rst_cnvs.to_csv("{0}.cnv.csv".format(sample),sep="\t")

    return


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
