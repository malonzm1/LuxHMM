#!/usr/bin/env python

import pickle
import numpy
import numpy as np
import numpy.random
import scipy.stats
import os

def get_stan_input(bsTot,bsC,D,bsEff,bsBEff,seqErr):
        sigmaB2 = 15
        alpha = 5
        beta = 5
        n_replicates = D.shape[0]
        n_predictors = D.shape[1]

        # data and init dictionaries
        data = {'bsBEff': bsBEff, 'n_replicates':n_replicates, 'n_predictors':n_predictors,
        'bsC': bsC, 'bsTot': bsTot,
        'X': D, 'sigmaB2': sigmaB2,
        'alpha': alpha, 'beta': beta,
        'bsEff': bsEff, 'seqErr': seqErr}

        return data

def savagedickey(locus,P):
        sigma2_B = 5
        beta = np.loadtxt("output.csv", delimiter=',', skiprows=38, usecols=(4,))
        density = scipy.stats.kde.gaussian_kde(np.transpose(beta),bw_method='scott')
        numerator = scipy.stats.multivariate_normal.pdf([0],mean=[0],cov=np.eye(1)*sigma2_B)
        denominator = density.evaluate([0])[0]
        bf=numerator/denominator
        return bf


