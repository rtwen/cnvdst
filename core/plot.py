#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.colors as colors
from plotnine import *
import numpy as np
import pandas as pd
import sys,re,os


def chromPlot(data,outFile):

    binLen = data.shape[0]
    bins = [ x+1 for x in range(0,binLen) ]
    data["bin"] = bins

    p1 = (ggplot(data) + geom_point(aes(x="bin",y="fixRatio"),color=colors.cnames["blue"]) +
             lims(y=(0, np.max(data["fixRatio"]))) + lims(x=(0,binLen)))

    p1.save(outFile)



def cnvPlot(cnvs,datas,gender):
    return



if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("python3 {0} in.csv out.png".format(sys.argv[0]))
        sys.exit(-1)

    inf = sys.argv[1]
    outf = sys.argv[2]
    datas = pd.read_csv(inf,sep="\t")
    chromPlot(datas,outf)












