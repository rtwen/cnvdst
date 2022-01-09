#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os,re,sys
import numpy as np
import pandas as pd
import pysam
from optparse import OptionParser


def parseCommand():
    usage = "usage: ./" + sys.argv[0] + "-i cnv.csv -o out.csv -r refgene.txt "
    version = "%prog 1.0"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-i", "--input", dest="input", help="input file")
    parser.add_option("-o", "--output", dest="output", help="output file")
    parser.add_option("-r", "--refgene", dest="refgene", help="ucsc refgene file longest transcript file",default="/tanbowen/db/GATK/longest_transcript.sort.txt.gz")

    return parser.parse_args()




def anno(cnvs,refgene):
    refGeneTbx = pysam.TabixFile(refgene)
    geneList  = []

    for index,row in cnvs.iterrows():
        genes = getGeneAnno(refGeneTbx,row["chrom"],row["start"],row["end"])
        geneList.append(genes)
    cnvs["gene"] = geneList

    refGeneTbx.close()

    return

def getGeneAnno(refTbx,chrom,start,end):

    start =int(start)
    end = int(end)
    glist = []
    gdict = dict()

    for row in refTbx.fetch(chrom,start,end):
        sflag = 0
        ratio = 0
        aflag  = 0
        exon = '-'
        cstart = 0
        cend = 0
        info = row.strip().split('\t')
        estart = info[9].split(',')
        eend = info[10].split(',')
        enum = int(info[8])
        if start > int(info[5]) or end < int(info[4]):
            exon = '-'
            ratio = 0
        elif start <=int(info[4]) and end >=int(info[5]):
            exon = "Exon1-" + str(enum)
            ratio = 1
            aflag = 1
        else:

            sflag = 0

            for i in range(0,len(estart)-1):
                if start <= int(estart[i]) and sflag == 0:
                    sflag = 1
                    cstart = i

                if sflag == 1 and end <= int(eend[i]):
                    cend = i
                    break
            lend = len(eend) - 2
            if end > int(eend[lend]) :
                cend = lend

            if re.search(r'\+',info[3],re.I):
                cstart += 1
                cend += 1
            elif re.search(r'\-',info[3],re.I):
                ctmp = cstart
                cstart = enum - cend
                cend = enum - ctmp
            exon = "Exon{0}-{1}".format(cstart,cend)
            ratio = round((cend-cstart +1)/enum,2)

        if not gdict.get(info[0]):
            if aflag == 1:
                glist.append(info[0])
            else:
                glist.append("{0}:{1}".format(info[0],exon))
            gdict[info[0]] = 1
    genes = '-'
    if len(glist) >=1:
        genes = ",".join(glist)

    return genes


if __name__ == '__main__':

    (options, args) = parseCommand()
    if options.input == None or options.output == None or options.refgene == None:
        print("see -h for help\n")
        sys.exit(-1)
    cnvs = pd.read_csv(options.input,sep="\t")
    anno(cnvs,options.refgene)
    cnvs.to_csv(options.output,sep="\t")




