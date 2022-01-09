#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os,re,sys
import numpy as np
import pandas as pd


def getBaseLine(data):
    row_id = 0
    mm = []
    stds= []
    for index,row in data.iterrows():
        ratios = row[3:]
        print(list(ratios))
        meanRatio = np.mean(ratios)
        std = np.std(ratios)
        mm.append(meanRatio)
        stds.append(std)
    data["std"] = stds
    data["mean"] = mm


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("python3 t.py inf")
        sys.exit(-1)

    data = pd.read_csv(sys.argv[1],sep="\s+",header=None)
    getBaseLine(data)
    print(data)
    data[['std','mean']] = data[['std','mean']].astype('float')
    data2 = data[data["std"] < 0.15]
    print(data2)
