#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os,re,sys
import numpy as np
import statsmodels.api as sm
import pandas as pd
import matplotlib.colors as colors
from plotnine import *


def gc_depth_lowess(data):
    lowess = sm.nonparametric.lowess
    #result = lowess(data["depth"], data["gc"],frac=0.25, it=3, delta=0.0)
    result = lowess(data["depth"], data["gc"],frac=0.25,delta=0.0)
    rst = pd.DataFrame(result)
    rst.columns = ["gc","depth"]
    #rst.mean
    return rst


def drawGCDepth(data,lowe,outFile):
    p1 = (ggplot(data) + geom_point(aes(x="gc",y="depth"),color=colors.cnames["blue"]) +
            geom_line(aes(x=lowe["gc"],y=lowe["depth"]),color=colors.cnames["black"],size=2))

    p1.save(outFile)



if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("python3 {0} inputfile outPng".format(sys.argv[0]))
        sys.exit(-1)
    infile = sys.argv[1]
    outPng = sys.argv[2]
    data = pd.read_csv(infile,sep="\t",names=["gc","depth"])
    lowe = gc_depth_lowess(data)
    #lowe[lowe.depth < 2000].loc[] = 0
    #print(lowe)
    lowe.to_csv("lowess.csv",sep="\t")
    gcmean = lowe.mean()


    #print(gcmean["gc"])
    #drawGCDepth(data,lowe,outPng)
