#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os,re,sys
import numpy as np
import pandas as pd


def gc_corrent(gcLowe,depth,win):

    k = 10000
    depth["gc"] = depth["gc"] * k
    meanDepth = depth["depth"].mean()

    dratio = mdepth_lowessDepth_Ratio(gcLowe,meanDepth)

    for index,row in depth.iterrows():
        row["depth"] = row["depth"] / dratio[int(row["gc"])]

    depth["gc"] = depth["gc"] / k

    return


def mdepth_lowessDepth_Ratio(gcLowe,meanDepth):

    dratio = dict()
    k = 10000
    gc = gcLowe["gc"] * k
    minGC = int(pd.DataFrame(gc).min()[0])
    maxGC = int(pd.DataFrame(gc).max()[0])

    for index,row in gcLowe.iterrows():
        i = int(row["gc"] * k)
        if row["depth"] <=0:
            dratio[i] = 1
        else:
            dratio[i] = meanDepth / row["depth"]

    for i in range(0,minGC):
        dratio[i] = dratio[minGC]

    for i in range(maxGC+1,k+1):
        dratio[i] = dratio[maxGC]

    for i in range(minGC+1,maxGC):

        if not dratio.get(i):
            lpos = i - 1
            rpos = i + 1
            lratio = dratio[lpos]
            while not dratio.get(rpos):
                rpos += 1
            rratio = dratio[rpos]
            dratio[i] = (lratio + (rpos- lpos + 0.0) * rratio) / (rpos - lpos + 1.0);

    return dratio



if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("python3 {0} infdepth gcLowess output\n".format(sys.argv[0]))
        sys.exit()

    depthf =sys.argv[1]
    lof = sys.argv[2]
    depth = pd.read_csv(depthf,sep="\t",names=["gc","depth"])
    loe = pd.read_csv(depthf,sep="\t",names=["gc","depth"])
    win = 200
    gc_corrent(loe,depth,win)
    depth.to_csv("gc_corrent.depth.csv",sep="\t")






