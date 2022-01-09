#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os,re,sys
import numpy as np
import pandas as pd

from hmmlearn import hmm



def fit_hmm(
        data, # normalised coverage array 
        variance, # variance per copy 
        variance_fixed,  # variance for the zero copy number state 
        min_copy_number=0,  # minimum copy number to consider in the model
        max_copy_number=6,  # maximum copy number to consider in the model 
        n_iter=1000000,  # number of iterations to perform when fitting the model
        params='st',  # parameters that can be changed through fitting 
        init_params='',  # parameters that are initialised from the data
        ):

    n_states = max_copy_number - min_copy_number
    means = []
    covars = []


    startprob_prior = []

    if n_states < 6:
        print("Error the states must >=6 ")
        sys.exit(-1)

    if n_states == 6:
        startprob_prior.extend([0.01,0.03,0.9,0.03,0.02,0.01])
        means.extend([[0], [0.5], [1.0], [1.5] ,[2.0], [2.5]])
        covars.extend([[0.15], [0.15], [0.2], [0.25] ,[0.25], [0.25]])
    else:
        startprob_prior.extend([0.01,0.03,0.9,0.03,0.02,0.01])
        startprob_prior.extend([0.00] * (n_states -6))
        means = means.extend([[0], [0.5], [1.0], [1.5], [2.0], [2.5] ])
        covars.extend([[0.15], [0.15], [0.2], [0.25] ,[0.25], [0.25]])
        for i in range(7,n_states):
            means.append([2.5 + (i-6)/2])
            covars.extend([0.25])

    means = np.array(means)

    #strtprob = 
    # construct the transition matrix
    transmat = np.zeros((n_states, n_states))
    transition_probability = 1 / (n_states + 1)
    transmat[:] = transition_probability
    transmat[np.diag_indices(n_states)] = 1-((n_states-1)*transition_probability)


    # construct means and covariance
    #means = np.array([[n] for n in range(min_copy_number, max_copy_number)])
    #covars = np.array([[variance*n + variance_fixed] for n in range(min_copy_number, max_copy_number)])

    # setup HMM 
    """
    model = hmm.GaussianHMM(n_states,
                        covariance_type='diag',
                        startprob_prior=startprob_prior,
                        n_iter=n_iter,
                        transmat_prior=transmat,
                        params=params,
                        init_params=init_params)
    """

    np.random.seed(2)
    model = hmm.GaussianHMM(n_states,
            covariance_type='diag',
            n_iter=n_iter,
            params=params,
            init_params=init_params)

    model.startprob_ = startprob_prior
    model.transmat_ = transmat

    model.means_ = means
    model.covars_ = covars

    # fit HMM
    obs = np.column_stack([data["fixRatio"]])
    model.fit(obs)
    print(model.covars_)
    print(model.means_)
    #print(model.startprob_)
    #print(model.transmat_)

    # predict hidden states
    h = model.predict(obs)
    return h

if __name__ =='__main__':


    if len(sys.argv) < 2:
        print("python3 {0} in.csv out.csv".format(sys.argv[0]))
        sys.exit(-1)
    inf = sys.argv[1]
    outf = sys.argv[2]
    datas = pd.read_csv(inf,sep="\t")
    re = fit_hmm(datas,0.1,0.01)
    datas["cn"] = re
    datas.to_csv(outf,sep="\t")



