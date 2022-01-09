#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os,re,sys
import numpy as np
import pandas as pd
from . import bamstat



def mergeCnv(data,gender):

    ploidy = bamstat.getPloidy(data,gender)
    isin = 0
    i = 0
    copyn = -1
    chrom = ''
    cnv = []

    for index,row in data.iterrows():

        if row["cn"] == ploidy[row["chrom"]]  and isin == 0:
            continue

        binLen = row["end"] - row["start"]

        if isin == 1:

            if row["cn"] != copyn or row["chrom"] !=chrom:
                i += 1
                chrom = row["chrom"]
                if row["cn"] == ploidy[row["chrom"]] :
                    isin = 0
                else:
                    isin = 1
                    cnv.append([row["chrom"],row["start"],row["end"],1,binLen,row["depth"] * binLen,
                        row["ratio"]*binLen, row["MeanRatio"]*binLen,row["fixRatio"]*binLen,row["cn"]])

                    copyn = row["cn"]
            else:
                cnv[i][2] = row["end"]
                cnv[i][3] += 1
                cnv[i][4] += binLen
                cnv[i][5] += row["depth"] * binLen
                cnv[i][6] += row["ratio"] * binLen
                cnv[i][7] += row["MeanRatio"] * binLen
                cnv[i][8] += row["fixRatio"] * binLen
        else :

            if row["cn"] != ploidy[row["chrom"]] :
                isin = 1
                chrom = row["chrom"]
                cnv.append([row["chrom"],row["start"],row["end"],1,binLen,row["depth"]*binLen,
                    row["ratio"]*binLen,row["MeanRatio"]*binLen,row["fixRatio"]*binLen,row["cn"]])
                copyn = row["cn"]

    for i in range(len(cnv)):

        cnv[i][5] = round(cnv[i][5]/cnv[i][4],2)
        cnv[i][6] = round(cnv[i][6]/cnv[i][4],2)
        cnv[i][7] = round(cnv[i][7]/cnv[i][4],2)
        cnv[i][8] = round(cnv[i][8]/cnv[i][4],2)

    result = pd.DataFrame(cnv,columns = ['chrom', 'start', 'end',"binNum","binLen","depth","ratio","ControlMeanRatio","fixRatio","cn"])

    return result


if __name__ == '__main__':
    if len(sys.argv) < 4 :
        print("python3 {0} in.cn.csv out.cnv.csv gender ".format(sys.argv[0]))
        sys.exit(-1)
    inf = sys.argv[1]
    outf = sys.argv[2]
    gender = sys.argv[3]
    data = pd.read_csv(inf,sep="\t")
    result =mergeCnv(data,gender)
    result.to_csv(outf,sep="\t")

