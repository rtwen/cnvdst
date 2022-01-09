#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,re,sys
import numpy as np
import statsmodels.api as sm
import pandas as pd
import matplotlib.colors as colors
from . import region
from plotnine import *
from scipy.stats import norm

CV_MAX = 0.15


def baseline_fix(data,baseline):

    baseline1 = baseline[["chrom","start","end","MeanRatio"]]
    chrdict = region.get_chrom_dict(data)
    chrdict1  = region.get_chrom_dict(baseline)
    fix_data = pd.DataFrame()

    if chrdict.get("chrY") or chrdict.get("Y"):

        fix_chrY = pd.DataFrame()

        if chrdict1.get("chrY") or chrdict1.get("Y"):
            baseline_chrY = baseline1[(baseline1["chrom"] =='chrY') | (baseline1["chrom"] =='Y')]
            #baseline_chrY = baseline1[[region.is_chrY(baseline1["chrom"][i]) for i in range(0,len(baseline1["chrom"]))]]
            data_chrY = data[(data["chrom"] =='chrY') | (data["chrom"] =='Y')]
            fix_chrY = ratioFilter(pd.merge(data_chrY,baseline_chrY,on = ['chrom','start','end']))
            fix_chrY["fixRatio"] = round(fix_chrY["ratio"] / fix_chrY["MeanRatio"],4)
        else:
            fix_chrY = data_chrY
            fix_chrY["MeanRatio"] = -1
            fix_chrY["fixRatio"] = fix_chrY["ratio"]
        if len(region.get_chrom_list(data)) >1:
            baseline_other = baseline1[(baseline1["chrom"] !='chrY') & (baseline1["chrom"] !='Y')]
            #baseline_other = baseline1[[not region.is_chrY(baseline1["chrom"][i]) for i in range(0,len(baseline1["chrom"]))]]
            data_other = data[(data["chrom"] !='chrY') & (data["chrom"] !='Y')]
            fix_data = ratioFilter(pd.merge(data_other,baseline_other,on = ['chrom','start','end']))
            fix_data["fixRatio"] = round(fix_data["ratio"] / fix_data["MeanRatio"],4)

        fix_data.append(fix_chrY)
    else:
        fix_data = ratioFilter(pd.merge(data,baseline1,on = ['chrom','start','end']))
        fix_data["fixRatio"] = round(fix_data["ratio"] / fix_data["MeanRatio"],4)


    return fix_data


def getBaseline(data):

    cvList = []
    mratioList = []
    stdList =[]
    # data : chrom start end GC ratio1 ratio2 ratio3...
    for index,row in data.iterrows():
        ratios = row[4:]
        ratios_clean_mad = ratios[mad_filter(ratios)]
        cv = getcv(ratios_clean_mad)
        meanRatio = round(np.mean(ratios_clean_mad),4)
        std = round(np.std(ratios_clean_mad),4)
        cvList.append(cv)
        mratioList.append(meanRatio)
        stdList.append(std)

    data["MeanRatio"] = mratioList
    data["std"] = stdList
    data["CV"] = cvList


def baselineFilter(baseline):
    fbsline = baseline[baseline["CV"] < CV_MAX]
    fbsline.reset_index()
    return fbsline


def ratioFilter(data,thread=0.05):
    data1 = data[data["MeanRatio"] >= thread]
    return data1

def getRatio(depth,gender):

    meanDepth = 0
    index = [region.is_autosomes(row["chrom"]) for i,row in depth.iterrows()]
    depth["index"] = index
    de = depth[index]
    if len(de["chrom"]) < 1:
        meanDepth = depth["depth"].mean()
    else:
        meanDepth = de["depth"].mean()
    meanDepth = round(depth["depth"].mean(),4)

    depth["ratio"] = round(depth["depth"] / meanDepth,4)

    return


def getControlRatio(depth,gender):

    meanDepth = 0
    index = [region.is_autosomes(row["chrom"]) for i,row in depth.iterrows()]
    de = depth[index]
    if len(de["chrom"]) < 1:
        meanDepth = depth["depth"].mean()
    else:
        meanDepth = de["depth"].mean()

    depth["ratio"] = depth["depth"] / meanDepth
    index_sex_chrom = [region.is_sex_chrom(depth["chrom"][i])  for i in range(0,len(depth["chrom"]))]
    index_chrY_chrom = [region.is_chrY(depth["chrom"][i])  for i in range(0,len(depth["chrom"]))]

    if gender == 'M':
        depth.loc[index_sex_chrom,'ratio'] *= 2
        #ratios[[region.is_sex_chrom(depth["chrom"][i])  for i in range(0,len(depth["chrom"]))]] *=2
    elif gender == 'F':
        depth.loc[index_chrY_chrom,'ratio'] = 1
        #ratios[[region.is_chrY(depth["chrom"][i])  for i in range(0,len(depth["chrom"]))]] =1

    """
    chrom_list = region.get_chrom_list(depth)
    for chrom in chrom_list:
        index = depth[depth.chrom == chrom].index.tolist()
        if gender == 'M' :
            if region.is_sex_chrom(chrom):
                depth["ratio"][index] = depth["dpeth"][index] * 2 / meanDepth
            else :
                depth["ratio"][index] = depth["depth"][index] / meanDepth
        elif gender == 'F':
            if not region.is_chrY(chrom):
                depth["ratio"][index] = depth["depth"][index] / meanDepth
    """

    return



def getcv(data):

    if len(data) <=1:
        return 0
    if type(data) is list:
        data = np.asarray(data)
    med = np.mean(data) + 0.001
    sd = np.std(data)
    cv = round(sd / med,4)

    return cv


def mad_filter(data,thresh=3.5):

    if len(data)<=1:
        return [True]

    if type(data) is list:
        data = np.asarray(data)

    med = np.median(data)
    #print(med)
    abs_dev = np.absolute(data - med + 0.00001)
    med_abs_dev = np.median(abs_dev)
    score = norm.ppf(0.75) * abs_dev / med_abs_dev
    return score<thresh


if __name__ == '__main__':
    main()
